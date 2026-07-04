# Pink Pigeon Immune-Gene Annotation Pipeline

Identify and annotate immune-related genes in the Pink Pigeon
(*Nesoenas mayeri*) genome by integrating three independent lines of evidence
and assigning each gene a confidence tier.

- **InterProScan** - protein-domain based immune GO-term filtering
- **KofamScan** - KEGG KO annotation restricted to immune pathways
- **OrthoFinder** - orthology to immune genes in three reference species
  (Mouse, Chicken, Zebra Finch), plus orthogroup enrichment statistics
- **Tiering** - integrate the three sources, enrich functional annotation,
  predict gene symbols, and produce summary figures

The final output is a gene-level table with an evidence **Tier** (1-3), predicted
symbols, orthology relationships, and orthogroup immune-enrichment statistics.

> **Built for:** this pipeline was developed for the Pink Pigeon genome
> assembly **bNesMay2.1 (GCA_963082525.1)** and its annotation
> **`Nesoenas_mayeri-GCA_963082525.1-2023_10-genes.gff3`** (ENA/INSDC).
> The chromosome-ID map in `config/chromosome_id_map.tsv` is specific to this
> assembly.
>
> It can be adapted to **other species or a new genome annotation** with minor
> changes: point the paths in `config/config.yaml` at the new genome FASTA,
> GFF3 and protein FASTA, update `config/chromosome_id_map.tsv` for the new
> assembly (or skip stage 00 if the IDs already match), and adjust the reference
> species set in `species:` together with the corresponding BioMart / curated
> inputs. The analysis logic itself is species-agnostic.

---

## 1. Requirements

**Python stages** (this repo): create the conda environment once.

```bash
conda env create -f envs/environment.yml
conda activate pigeon_immune
```

**Heavy annotation tools** are loaded as HPC modules inside the SLURM scripts,
not installed by the environment above:

| Tool | Version | Used by |
|------|---------|---------|
| InterProScan | 5.73-104.0 | `01_interproscan/interproscan_run.sh` |
| seqkit | (module) | FASTA splitting for InterProScan |
| KofamScan (`exec_annotation`) | conda env `pigeon_immune` | `02_kofamscan/run_kofam.sh` |
| OrthoFinder | 3.0.1b1 | `03_orthofinder/run_orthofind.sh` |

---

## 2. Configuration

All data paths live in [`config/config.yaml`](config/config.yaml). Edit the
absolute paths under `reference:` and `external_data:` once for your
environment. Every Python script reads its defaults from this file via
`--config config/config.yaml`; any value can be overridden with a CLI flag
(e.g. `--input`, `--output`).

The SLURM scripts keep a small **USER CONFIG** block at the top for HPC-specific
values (partition, memory, email, module versions). Keep those in sync with the
config file by hand.

---

## 3. External data to prepare manually

These files are **not** produced by the pipeline and must be obtained before
running (paths set under `external_data:` in the config):

| File | Source |
|------|--------|
| `go_terms_immune_system_process.txt` | AmiGO: all children of GO:0002376 (immune system process) plus "immune" keyword hits in MF/CC, comma-separated |
| `processed_chicken_biomart.csv` | Ensembl BioMart (Chicken immune genes from the Avian Immune Database) |
| `processed_zebrafinch_biomart.csv` | Ensembl BioMart (Zebra Finch immune genes) |
| `processed_mouse_biomart.csv` | Ensembl BioMart (Mouse gene IDs from the curated list) |
| `ImmuneGeneFunction_20240520.csv` | Published curated mouse immune-gene list (Category1, Subcategory, UniProt_function) |
| `gene_info.csv` | Conserved avian immune-gene functional database (unpublished) |
| `Nesoenas_genes_by_name.tsv` | Pink Pigeon native gene-symbol table from the genome annotation (Symbol `<TAB>` Gene ID) |

The reference genome FASTA, GFF3 and the two protein FASTAs come from the genome
annotation project (see `reference:` in the config).

---

## 4. Running the pipeline

Run stage by stage. The heavy tools are SLURM array/batch jobs that take hours
to days; submit them, wait, then run the local Python steps. The driver
[`run_pipeline.sh`](run_pipeline.sh) documents the exact order and commands.

### Stage 00 - preprocessing
```bash
bash 00_preprocess/gffid_change.sh -g <raw.gff3> -m config/chromosome_id_map.tsv \
     -r NM_genes.renamed.gff3 -a NM_genes.final.gff3
# optional QC (inspect the output by hand):
bash 00_preprocess/seqid_list.sh -f <genome.fasta> -g <raw.gff3>
```

### Stage 01 - InterProScan  (HPC)
```bash
# 1. MANUAL: split the protein FASTA first
module load seqkit && seqkit split -p 8 -O ipr_chunks CDS_protein_seq.fa
# 2. array job (one InterProScan run per chunk)
sbatch 01_interproscan/interproscan_run.sh
# 3. merge the chunk outputs
sbatch --array=9 01_interproscan/interproscan_run.sh
# 4. filter to immune GO terms (local)
python 01_interproscan/immunewash.py --config config/config.yaml
```

### Stage 02 - KofamScan  (HPC)
```bash
python 02_kofamscan/prepare_db.py --config config/config.yaml   # download KEGG JSON
sbatch 02_kofamscan/run_kofam.sh                                # HPC
python 02_kofamscan/extract_map.py     --config config/config.yaml
python 02_kofamscan/filter_kofam.py    --config config/config.yaml
python 02_kofamscan/map_to_gene.py     --config config/config.yaml
python 02_kofamscan/add_ko_metadata.py --config config/config.yaml
```

### Stage 03 - OrthoFinder  (HPC)
```bash
python 03_orthofinder/data_merge_new.py --config config/config.yaml
sbatch 03_orthofinder/run_orthofind.sh                          # HPC
python 03_orthofinder/final_filtering_v2.py       --config config/config.yaml
python 03_orthofinder/compute_og_stats.py         --config config/config.yaml
python 03_orthofinder/compute_orthology_stats.py  --config config/config.yaml
```

### Stage 04 - tiering  (local)
```bash
python 04_tiering/merge_immune_annotations.py           --config config/config.yaml
python 04_tiering/add_gene_info_3.py                    --config config/config.yaml
python 04_tiering/process_ppg_genes.py                  --config config/config.yaml
python 04_tiering/merge_og_stats_to_final_gene_table.py --config config/config.yaml
python 04_tiering/merge_orthology_to_final_gene_table.py --config config/config.yaml
python 04_tiering/visualization.py                      --config config/config.yaml
```

Final table: `PinkPigeon_Immune_Predict_Result_Final_with_OG_and_Orthology.csv`
plus 7 PNG figures in the configured `figure_dir`.

---

## 5. Manual steps and checkpoints

These are **not** automated and need human action or review:

1. **External data downloads** (section 3) - BioMart, AmiGO, curated list,
   `gene_info.csv`, PPG native symbols.
2. **KEGG JSON download** (`prepare_db.py`) - needs internet. HPC compute nodes
   are often offline; download the JSON elsewhere and copy it in if it fails.
3. **FASTA splitting** before the InterProScan array job (`seqkit split`).
4. **ID consistency check** (`seqid_list.sh`) - a diagnostic; inspect
   `ids_compare.csv` by eye before trusting the annotation.
5. **Symbol comparison review** - in `process_ppg_genes.py` the
   `Symbol_Comparison` categories `Paralog_Likely` and `Mismatch` flag cases
   that should be curated by hand.

---

## 6. Repository layout

```
immune-annotation-pipeline/
├── README.md
├── run_pipeline.sh              driver documenting run order
├── config/
│   ├── config.yaml              all paths and parameters
│   └── chromosome_id_map.tsv    simple ID -> INSDC accession
├── envs/environment.yml         conda environment
├── 00_preprocess/              chromosome-ID normalization + QC
├── 01_interproscan/            InterProScan run + immune GO filtering
├── 02_kofamscan/               KofamScan run + KEGG immune filtering
├── 03_orthofinder/             OrthoFinder + orthogroup statistics
├── 04_tiering/                 evidence integration + symbols + figures
├── docs/PIPELINE_OVERVIEW.md   detailed per-script I/O reference
└── extras/                     side scripts not part of the core pipeline
```

See [`docs/PIPELINE_OVERVIEW.md`](docs/PIPELINE_OVERVIEW.md) for a per-script
description of inputs, outputs and logic.
