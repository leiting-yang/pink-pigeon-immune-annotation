# Immune-Gene Annotation Pipeline

Identify and annotate immune-related genes in a genome by integrating three
independent lines of evidence and assigning each gene a confidence tier. The
pipeline was originally developed for the Pink Pigeon (*Nesoenas mayeri*), but
its analysis logic is species-agnostic and driven entirely by
`config/config.yaml`, so it can be reused for another species without editing
code.

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
> The analysis logic is **species-agnostic**: the target and reference species,
> their BioMart inputs and the symbol-prediction rules are all driven by the
> `species:` block in `config/config.yaml`. See
> [section 7 (Reusing for another species)](#7-reusing-the-pipeline-for-another-species).

---

## 1. Requirements

**Python stages** (this repo): create the conda environment once.

```bash
conda env create -f envs/environment.yml
conda activate immune_annotation
```

**Heavy annotation tools.** InterProScan, seqkit and OrthoFinder are loaded as
HPC modules inside the SLURM scripts. KofamScan is installed by the conda
environment above (it is rarely available as a module) but still needs its HMM
database downloaded once (see "KofamScan database" below).

| Tool | Version | Provided by | Used by |
|------|---------|-------------|---------|
| InterProScan | 5.73-104.0 | HPC module | `01_interproscan/interproscan_run.sh` |
| seqkit | (module) | HPC module | FASTA splitting for InterProScan |
| KofamScan (`exec_annotation`) | (bioconda) | conda env `immune_annotation` | `02_kofamscan/run_kofam.sh` |
| OrthoFinder | 3.0.1b1 | HPC module | `03_orthofinder/run_orthofind.sh` |

Each SLURM script `module load`s its own prerequisites (e.g. `interproscan_run.sh`
loads `openjdk` before `interproscan`). To check every dependency resolves on
your cluster before running, use `module spider <tool>` and the helper
`test/check_deps.sh`.

**KofamScan database.** KofamScan needs the KO HMM profiles + `ko_list`, which are
NOT installed by conda. Download them once and point `KOFAM_DB_CONFIG`
(`config/cluster.sh`) at a `config.yml` that references them:

```bash
mkdir -p <db_dir> && cd <db_dir>
wget https://www.genome.jp/ftp/db/kofam/ko_list.gz && gunzip ko_list.gz
wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz && tar xzf profiles.tar.gz
cat > config.yml <<EOF
profile: <db_dir>/profiles
ko_list: <db_dir>/ko_list
EOF
```

---

## 2. Configuration

All data paths live in [`config/config.yaml`](config/config.yaml). Edit the
absolute paths under `reference:` and `external_data:`, and the `species:` block
(target + reference species with their BioMart files), once for your
environment. Every Python script reads its defaults from this file via
`--config config/config.yaml`; any value can be overridden with a CLI flag
(e.g. `--input`, `--output`).

The three SLURM scripts (`interproscan_run.sh`, `run_kofam.sh`,
`run_orthofind.sh`) cannot parse YAML, so their editable values (protein FASTA
name, work dirs, conda/OrthoFinder paths) live in
[`config/cluster.sh`](config/cluster.sh), which they `source`. Edit that one
file instead of the scripts. Submit the SLURM jobs from the repo root (they find
`config/cluster.sh` via `$SLURM_SUBMIT_DIR`), or set `CLUSTER_CONFIG` to its
absolute path. The only thing left inside the SLURM scripts is the `#SBATCH`
header (partition, memory, email), which SLURM requires to be literal.

So in total you edit **two** config files - `config/config.yaml` (Python stages)
and `config/cluster.sh` (SLURM scripts) - and none of the analysis scripts.

---

## 3. Input files to prepare

Nothing in this section is produced by the pipeline; you must obtain every file
below **before** running and point the config at it. The examples in the tables
are the Pink Pigeon filenames shipped in the default `config/config.yaml` —
replace them with your own.

### 3.1 Complete input inventory

**A. Target species (from your genome annotation project)** — set under
`reference:` and `external_data:`:

| File (example) | Config key | What it is |
|------|------------|------------|
| `GCA_..._genomic.fasta` | `reference.genome_fasta` | Genome nucleotide FASTA. Sequence IDs are whatever your assembly uses (INSDC accessions, RefSeq, or simple `1/2/Z`). |
| `Nesoenas_...-genes.gff3` | `reference.gff3_raw` | Gene-model annotation GFF3 for the same assembly. |
| `PinkPigeon.faa` | `reference.protein_fasta` | Target protein FASTA. **Produced** by stage 00 `extract_proteins.sh` (which also cleans gffread's `.`/`*` characters so InterProScan/OrthoFinder/KofamScan accept it); you can also supply your own. Feeds InterProScan and KofamScan. |
| `PinkPigeon_genes_by_name.tsv` | `external_data.native_symbols` | The annotation's own gene symbols: 2 columns, `Symbol <TAB> GeneID`. **Produced by default** by stage 00 `extract_native_symbols.sh` from the GFF3; or supply your own (an empty file is fine if you have none). Used to compare predicted vs. native symbols. |

**B. Reference species (one set per entry in `species.reference`)** — used by
OrthoFinder and the reference-symbol lookup:

| File (example) | Config key | Source |
|------|------------|--------|
| `<Species>.fa` protein FASTA | *(OrthoFinder input dir, see 3.2)* | That species' proteome. **The FASTA basename must equal the species `name`** in the config (e.g. `Chicken.fa`), because OrthoFinder derives its column names from the filename. |
| `processed_<species>_biomart.csv` | `species.reference[].biomart` | Ensembl BioMart export of that species' immune genes (gene ID, symbol, description, GO terms). Required for every reference species. |
| `ImmuneGeneFunction_20240520.csv` | `species.reference[].curated_list` | **Optional** published/curated functional list for one species (columns `GeneStableID`, `Category1`, `Subcategory`, `UniProt_function`). Only the Pink Pigeon setup attaches one, to Mouse. |

**C. Immune reference sets / databases** — set under `external_data:`:

| File (example) | Config key | Source |
|------|------------|--------|
| `go_terms_immune_system_process.txt` | `external_data.go_terms_immune` | AmiGO export: all children of GO:0002376 plus "immune" MF/CC hits, as a single comma-separated list. Defines the InterProScan immune filter. |
| `gene_info.csv` | `external_data.gene_info` | **Optional** conserved immune-gene functional database (avian-oriented, keyed by `entrezgene_accession`). Set `species.gene_info_enrichment: false` to skip it for non-avian work. |

Small copies of the reference tables in **B** and **C** (the BioMart/immune-gene
CSVs and the GO-term list) are bundled in [`reference_data/`](reference_data/);
see [`reference_data/CITATIONS.md`](reference_data/CITATIONS.md) for their
provenance and for how to download the large reference proteome FASTAs (Ensembl).

**D. Downloaded once, not in the repo:**
- **KofamScan HMM database** (profiles + `ko_list`) — see "KofamScan database" in
  section 1. Point `KOFAM_DB_CONFIG` (`config/cluster.sh`) at its `config.yml`.
- **KEGG `ko00001.json`** hierarchy — fetched by `02_kofamscan/prepare_db.py`. If
  the compute node has no internet, download it elsewhere and drop it into the
  KofamScan `processing_dir`.

### 3.2 Where to put the files (directory organization)

The shipped config uses **absolute paths for every input**, so inputs can live
anywhere (typically the workspace root) and each stage finds them regardless of
its `work_dir`. Outputs are written into **per-stage folders** under
`<workspace>/outputs/`:

```
workspace/
├── <genome>.fasta  <genome>.gff3  <Target>.faa            (target inputs)
├── go_terms_immune_system_process.txt  gene_info.csv       (external data)
├── processed_*_biomart.csv  ImmuneGeneFunction_*.csv       (reference tables)
├── <Target>_genes_by_name.tsv                              (native symbols)
├── ortho_input/   <Target>.faa Mouse.faa Chicken.faa ZebraFinch.faa
│   └── OrthoFinder/Results_<name>/   (OrthoFinder output; -n makes the name fixed)
└── outputs/
    ├── 01_interproscan/   ipr_chunks/  ipr_out/  interproscan_immune_results.csv
    ├── 02_kofamscan/       result_kofam_detail.txt  step*.txt  final_kofam_*.tsv
    ├── 03_orthofinder/     master_lookup_table.csv  Final_Filtered_List*  OG_stats_*
    └── 04_tiering/         Immune_*.csv  figures/
```

Notes:
- Keep the per-stage `*_WORK_DIR` values in `config/cluster.sh` in sync with the
  matching `work_dir` in `config/config.yaml`, so the SLURM raw outputs
  (`ipr_out/`, `result_kofam_detail.txt`) land where the Python stages look.
- The **OrthoFinder input directory** (`ORTHO_INPUT_DIR` / `orthofinder.input_dir`)
  must contain **exactly** the target + reference protein FASTAs, each named
  `<SpeciesName>.fa`/`.faa` to match the config. Nothing else.
- Create the output folders once before running:
  `mkdir -p <workspace>/outputs/{01_interproscan,02_kofamscan,03_orthofinder,04_tiering}`
  (the SLURM scripts also `mkdir -p` their own).

---

## 4. Running the pipeline

Run stage by stage. The heavy tools are SLURM array/batch jobs that take hours
to days; submit them, wait, then run the local Python steps. The driver
[`run_pipeline.sh`](run_pipeline.sh) documents the exact order and commands.

### Stage 00 - preprocessing
```bash
# extract the protein FASTA from the genome + GFF3 (input to stages 01 and 02).
# extract_proteins.sh also cleans gffread's '.'/'*' so all downstream tools accept it.
bash 00_preprocess/extract_proteins.sh -g <genome.fasta> -a <raw.gff3> -o PinkPigeon.faa

# native gene symbols from the GFF3 (ON by default; skip if you have your own,
# or if your GFF has no symbols). Writes the 2-column native-symbol table.
bash 00_preprocess/extract_native_symbols.sh -g <raw.gff3> -o PinkPigeon_genes_by_name.tsv

# OPTIONAL: normalize chromosome IDs. Only needed when the genome FASTA and the
# GFF3 use different sequence IDs (e.g. an Ensembl GFF with simple IDs). Most
# NCBI RefSeq assemblies already match, so you can skip this and use the raw
# GFF3 directly. Run the QC below first to check.
bash 00_preprocess/gffid_change.sh -g <raw.gff3> -m config/chromosome_id_map.tsv \
     -r NM_genes.renamed.gff3 -a NM_genes.final.gff3

# QC (inspect ids_compare.csv by hand): are the FASTA and GFF3 sequence IDs consistent?
bash 00_preprocess/seqid_list.sh -f <genome.fasta> -g <raw.gff3>
```

### Stage 01 - InterProScan  (HPC)
```bash
# 1. MANUAL: split the protein FASTA into the stage work dir's ipr_chunks/
#    (must be <IPR_WORK_DIR>/ipr_chunks, e.g. outputs/01_interproscan/ipr_chunks)
module load seqkit
seqkit split -p 8 -O <workspace>/outputs/01_interproscan/ipr_chunks <workspace>/PinkPigeon.faa
# 2. array job (one InterProScan run per chunk). Pass CLUSTER_CONFIG so the job
#    finds config/cluster.sh (or submit from the repo root).
sbatch --export=ALL,CLUSTER_CONFIG=$PWD/config/cluster.sh 01_interproscan/interproscan_run.sh
# 3. merge the chunk outputs
sbatch --array=9 --export=ALL,CLUSTER_CONFIG=$PWD/config/cluster.sh 01_interproscan/interproscan_run.sh
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
python 04_tiering/process_symbols.py                    --config config/config.yaml
python 04_tiering/merge_og_stats_to_final_gene_table.py --config config/config.yaml
python 04_tiering/merge_orthology_to_final_gene_table.py --config config/config.yaml
python 04_tiering/visualization.py                      --config config/config.yaml
```

Final table: `Immune_Predict_Result_Final_with_OG_and_Orthology.csv`
plus 7 PNG figures in the configured `figure_dir`. Every column of the final
table is documented in [`docs/OUTPUT_COLUMNS.md`](docs/OUTPUT_COLUMNS.md).

---

## 5. Manual steps and checkpoints

These are **not** automated and need human action or review:

1. **External data downloads** (section 3) - BioMart, AmiGO, curated list,
   `gene_info.csv`, target native symbols.
2. **KEGG JSON download** (`prepare_db.py`) - needs internet. HPC compute nodes
   are often offline; download the JSON elsewhere and copy it in if it fails.
3. **FASTA splitting** before the InterProScan array job (`seqkit split`).
4. **ID consistency check** (`seqid_list.sh`) - a diagnostic; inspect
   `ids_compare.csv` by eye before trusting the annotation.
5. **Symbol comparison review** - in `process_symbols.py` the
   `Symbol_Comparison` categories `Paralog_Likely` and `Mismatch` flag cases
   that should be curated by hand.

---

## 6. Repository layout

```
immune-annotation-pipeline/
├── README.md
├── CHANGELOG.md                 bug fixes / improvements (see for known issues)
├── run_pipeline.sh              driver documenting run order
├── pipeline_common.py          shared species/column resolution helpers
├── config/
│   ├── config.yaml              all paths, species and parameters (Python stages)
│   ├── cluster.sh               shell-side values for the SLURM scripts
│   └── chromosome_id_map.tsv    simple ID -> INSDC accession
├── reference_data/             bundled small reference tables + CITATIONS.md
├── envs/environment.yml         conda environment (incl. kofamscan)
├── 00_preprocess/              protein extraction, native symbols, chr-ID QC
├── 01_interproscan/            InterProScan run + immune GO filtering
├── 02_kofamscan/               KofamScan run + KEGG immune filtering
├── 03_orthofinder/             OrthoFinder + orthogroup statistics
├── 04_tiering/                 evidence integration + symbols + figures
├── test/                       dependency check + example per-species config
│   ├── check_deps.sh           verify every external dependency resolves
│   ├── config.yaml  cluster.sh example config for a second-species run
└── docs/
    ├── PIPELINE_OVERVIEW.md    detailed per-script I/O reference
    └── OUTPUT_COLUMNS.md       column dictionary for the final table
```

See [`docs/PIPELINE_OVERVIEW.md`](docs/PIPELINE_OVERVIEW.md) for a per-script
description of inputs, outputs and logic, and
[`docs/OUTPUT_COLUMNS.md`](docs/OUTPUT_COLUMNS.md) for what every column in the
final table means.

---

## 7. Reusing the pipeline for another species

The pipeline is species-agnostic. All species-specific behaviour is driven by
the `species:` block in `config/config.yaml`, so no code needs editing.

To adapt it:

1. **Set the target and reference species** in `species:`:
   ```yaml
   species:
     target: MyBird
     reference:
       - name: Chicken
         biomart: chicken_biomart.csv
         curated_list: chicken_curated.csv   # optional
       - name: ZebraFinch
         biomart: zebrafinch_biomart.csv
   ```
2. **Name the OrthoFinder input FASTAs to match** these names exactly. The
   species columns in `Orthogroups.tsv` come from the FASTA basenames, so the
   target FASTA must be `MyBird.fa`, the references `Chicken.fa`,
   `ZebraFinch.fa`, etc. This is the one hard requirement.
3. **Provide a BioMart CSV** for each reference species (same column layout as
   the Pink Pigeon inputs). A `curated_list` is optional per species.
4. **Tune symbol prediction** (optional): `consensus_group` is the set of
   species whose intersection is treated as high-confidence (use `[]` to
   disable); `prediction_priority` is the single-species fallback order.
5. **Point `reference:` paths** at the new genome FASTA, GFF3 and protein FASTAs.
6. **Chromosome renaming is optional:** stage 00 `gffid_change.sh` is only needed
   when the genome FASTA and GFF3 use different sequence IDs. Skip it (and the
   `chromosome_id_map.tsv`) if they already match, e.g. most NCBI RefSeq
   assemblies.
7. **Set `orthofinder.orthologues_dir`** to your OrthoFinder output for the
   target; the final folder is named `Orthologues_<target>` (defaults to that if
   left unset).
8. **Non-avian work:** set `gene_info_enrichment: false` to skip the
   avian-oriented `gene_info.csv` enrichment.

Everything downstream (columns, orthology summaries, figures) follows the
configured species automatically.

---

## 8. Source-specific blocks you may need to review

The `species:` block covers naming and species logic, but a few steps make
assumptions about **ID conventions and external-database formats** that are not
config-driven. None of these block a first run; review them if your data comes
from a different source (e.g. NCBI RefSeq instead of Ensembl/ENA) or a different
tool version. They are listed here rather than generalized because the right
value depends on your specific inputs.

1. **Transcript/gene ID-prefix cleaning.** Several scripts strip
   Ensembl/ENA-style prefixes (`transcript:`, `mRNA:`, `gene:`, `transcript_`,
   `gene_`) so the transcript and gene IDs coming from the GFF3, InterProScan,
   KofamScan and OrthoFinder line up on the same join key. Files:
   `02_kofamscan/extract_map.py`, `02_kofamscan/map_to_gene.py`,
   `04_tiering/merge_immune_annotations.py`,
   `04_tiering/merge_og_stats_to_final_gene_table.py`,
   `04_tiering/merge_orthology_to_final_gene_table.py`,
   `03_orthofinder/compute_orthology_stats.py`. If your annotation uses a
   different convention (e.g. NCBI `rna-`, `gene-`, or already-bare IDs), adjust
   these `.replace(...)` calls so the transcript-to-gene join keys match. If your
   IDs are already consistent across tools, the replacements are harmless no-ops.

2. **Ensembl protein-version suffix.** `03_orthofinder/final_filtering_v2.py`
   strips a trailing `.N` version from reference protein IDs (e.g. `ENSMUSP...9`).
   Only relevant to Ensembl-style versioned IDs; harmless otherwise.

3. **Chromosome-ID map (optional feature).** `config/chromosome_id_map.tsv` and
   `00_preprocess/gffid_change.sh` exist only to reconcile a genome FASTA and a
   GFF3 that use *different* sequence IDs (e.g. Ensembl `1/2/Z` vs INSDC
   accessions). Most NCBI RefSeq assemblies already match, so you can **skip
   stage 00 renaming entirely** and point `gff3_raw` at the raw GFF3. The map is
   assembly-specific and is not required by the rest of the pipeline. Run
   `00_preprocess/seqid_list.sh` first if you are unsure whether the IDs match.

4. **Immune GO-term whitelist tokens.** `01_interproscan/immunewash.py` removes
   InterProScan source tags such as `(InterPro,PANTHER)` from the GO column
   before matching. If your InterProScan version formats the GO column
   differently, check that token list.

5. **KEGG immune-pathway selection.** `02_kofamscan/filter_kofam.py` finds immune
   KOs by walking the KEGG hierarchy through the literal node names
   `09150 Organismal Systems` and `09151 Immune system`. If KEGG renames those
   nodes in a future `ko00001.json`, update the strings. (The score cutoff uses
   KofamScan's own per-KO adaptive threshold, so it needs no tuning.)

6. **KEGG grouped-symbol expansion.** `04_tiering/process_symbols.py`
   `expand_kegg_symbol` splits compact KO symbols like `GRK4_5_6` into
   `GRK4/GRK5/GRK6`. This is a heuristic for KEGG's naming style; adjust it if
   your KO symbols follow a different pattern.
