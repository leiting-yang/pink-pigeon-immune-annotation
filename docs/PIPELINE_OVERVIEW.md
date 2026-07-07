# Pipeline Overview

Per-script reference for the immune-gene annotation pipeline. Filenames and
species below are from the reference Pink Pigeon (*Nesoenas mayeri*) setup —
target species with Mouse (*Mus musculus*), Chicken (*Gallus gallus*) and Zebra
Finch (*Taeniopygia guttata*) as references — but the flow is the same for any
species configured in `config/config.yaml`.

> This overview reflects the **complete** flow, including the orthogroup
> enrichment-statistics layer (`compute_og_stats.py`,
> `merge_og_stats_to_final_gene_table.py`) that connects OrthoFinder filtering to
> the tiering stage. Earlier drafts of the overview omitted this layer.

---

## Full flow

```
Stage 00  extract_proteins.sh ─► PinkPigeon.faa (protein FASTA, feeds 01 + 02)
          gffid_change.sh     ─► renamed / autosome GFF3
          seqid_list.sh       ─► ids_compare.csv (manual QC)

Stage 01  interproscan_run.sh ─► interproscan_merged.tsv
          immunewash.py       ─► interproscan_immune_results.csv

Stage 02  prepare_db.py     ─► ko00001.json
          run_kofam.sh      ─► result_kofam_detail.txt
          extract_map.py    ─► transcript_to_gene.txt
          filter_kofam.py   ─► step1/2/3 *.txt
          map_to_gene.py    ─► final_kofam_gene_level.tsv
          add_ko_metadata.py─► final_kofam_annotated.tsv

Stage 03  data_merge_new.py         ─► master_lookup_table.csv
          run_orthofind.sh          ─► Orthogroups.tsv, GeneCount, Duplications, Orthologues/
          final_filtering_v2.py     ─► Final_Filtered_List.csv
          compute_og_stats.py       ─► ..._with_OG_stats.csv (+ OG_stats_summary.tsv)
          compute_orthology_stats.py─► ..._with_OG_stats_and_orthology.csv

Stage 04  merge_immune_annotations.py           ─► Immune_Gene_Master_List.csv
          add_gene_info_3.py                    ─► Immune_Annotated_Final.csv
          process_symbols.py                    ─► Immune_Predict_Result_Final.csv
          merge_og_stats_to_final_gene_table.py ─► ..._Final_with_OG_stats.csv
          merge_orthology_to_final_gene_table.py─► ..._Final_with_OG_and_Orthology.csv
          visualization.py                      ─► 7 PNG figures
```

---

## External data

| File | Source | Used by |
|------|--------|---------|
| `Nesoenas_mayeri-...-genes.gff3` | ENA/INSDC (GCA_963082525.1) | preprocess, tiering |
| `GCA_963082525.1_..._genomic.fasta` | ENA/INSDC | protein extraction, seqid QC |
| `PinkPigeon.faa` | extract_proteins.sh (gffread) | InterProScan / KofamScan |
| `go_terms_immune_system_process.txt` | AmiGO (GO:0002376 subtree + "immune" MF/CC) | immunewash |
| `processed_{chicken,zebrafinch,mouse}_biomart.csv` | Ensembl BioMart | data_merge_new |
| `ImmuneGeneFunction_20240520.csv` | published curated mouse list | data_merge_new |
| `gene_info.csv` | unpublished avian immune-gene database | add_gene_info_3 |
| `Nesoenas_genes_by_name.tsv` | genome annotation (native symbols) | add_gene_info_3 |

---

## Stage 00 - preprocessing

### `extract_proteins.sh`
Translate the genome + GFF3 into the target-species protein FASTA with gffread
(`gffread -y PinkPigeon.faa -g <genome.fasta> <gff3>`). This FASTA is the input
to both InterProScan and KofamScan.
- In: genome FASTA, raw GFF3
- Out: `PinkPigeon.faa`

### `gffid_change.sh`  *(optional)*
Replace simple chromosome IDs with INSDC accessions using
`config/chromosome_id_map.tsv`; also emit an autosome-only GFF3 (Z/W removed).
Sex chromosomes are identified by their **original** IDs before renaming. Only
needed when the genome FASTA and GFF3 use different sequence IDs; skip it (and
the map) when they already match, e.g. most NCBI RefSeq assemblies. The core
pipeline uses the raw GFF3 regardless.
- In: raw GFF3, chromosome map TSV
- Out: `NM_genes.renamed.gff3` (all), `NM_genes.final.gff3` (autosomes)

### `seqid_list.sh`  *(manual diagnostic)*
Compare sequence IDs between the genome FASTA and the GFF3.
- In: genome FASTA, GFF3 → Out: `ids_compare.csv` (inspect by hand)

---

## Stage 01 - InterProScan

### `interproscan_run.sh`  *(SLURM array)*
InterProScan (`-f tsv -dp -goterms`) per FASTA chunk, then merge. The `seqkit`
split is a manual prerequisite; the merge runs as a separate array index.
- In: `PinkPigeon.faa` (split into chunks)
- Out: `ipr_out/interproscan_merged.tsv`

### `immunewash.py`
Keep proteins whose GO annotations hit the immune whitelist.
- In: merged InterProScan TSV, immune GO-term whitelist
- Out: `interproscan_immune_results.csv` (+ `immune_go_term` column)

---

## Stage 02 - KofamScan

### `prepare_db.py`
Download the KEGG `ko00001` hierarchy JSON (manual fallback if node offline).
- Out: `ko00001.json`

### `run_kofam.sh`  *(SLURM)*
KofamScan `exec_annotation -f detail-tsv`.
- In: `PinkPigeon.faa`, KofamScan DB config → Out: `result_kofam_detail.txt`

### `extract_map.py`
Transcript → gene map from GFF3 mRNA lines.
- Out: `transcript_to_gene.txt`

### `filter_kofam.py`
Parse detail format; keep score > threshold; keep immune-pathway KOs;
deduplicate per transcript by lowest E-value.
- In: `result_kofam_detail.txt`, `ko00001.json`
- Out: `step1_filtered_score.txt`, `step2_immune_only.txt`, `step3_final_deduplicated.txt`

### `map_to_gene.py`
Lift transcript-level KO to gene level; dedup per gene by lowest E-value.
- Out: `final_kofam_gene_level.tsv`

### `add_ko_metadata.py`
Add `KO_Gene_Symbol` and `KO_Gene_Name` from the KEGG hierarchy.
- Out: `final_kofam_annotated.tsv`

---

## Stage 03 - OrthoFinder

### `data_merge_new.py`
Build the reference immune-protein lookup: clean text, join mouse BioMart with
the curated list, tag `Immune_Source`, concatenate all three species.
- Out: `reference_mouse_final.csv`, `master_lookup_table.csv`

### `run_orthofind.sh`  *(SLURM)*
OrthoFinder on the four-species protein set.
- Out: `Orthogroups.tsv`, `Orthogroups.GeneCount.tsv`, `Duplications.tsv`,
  `Orthologues/`

### `final_filtering_v2.py`
Exist-to-keep: keep all Pink Pigeon proteins in any orthogroup containing a
reference immune protein.
- In: `Orthogroups.tsv`, `master_lookup_table.csv`
- Out: `Final_Filtered_List.csv`

### `compute_og_stats.py`
Per candidate orthogroup: immune fraction/enrichment vs a genome-wide
background, one-sided Fisher exact test (hypergeometric, no scipy) + BH FDR,
normalized Shannon entropy, gene-copy variance/SD, duplication density.
- In: `Orthogroups.tsv`, `Orthogroups.GeneCount.tsv`, `Duplications.tsv`,
  `master_lookup_table.csv`, `Final_Filtered_List.csv`
- Out: `OG_stats_summary.tsv`, `Final_Filtered_List_with_OG_stats.csv`

### `compute_orthology_stats.py`
Classify each transcript's orthology type (1:1, 1:N, N:1, N:M) per species and
combine with reference immune status.
- In: `{target}__v__{species}.tsv`, `master_lookup_table.csv`,
  `..._with_OG_stats.csv`
- Out: `..._with_OG_stats_and_orthology.csv`

---

## Stage 04 - tiering

### `merge_immune_annotations.py`
Gene-level outer join of the three evidence sources; assign Tier
(3 sources = Tier 1, 2 = Tier 2, 1 = Tier 3).
- In: `interproscan_immune_results.csv`, `Final_Filtered_List.csv`,
  `final_kofam_annotated.tsv`, GFF3
- Out: `Immune_Gene_Master_List.csv`

### `add_gene_info_3.py`
Enrich with `gene_info.csv` via an all-in symbol strategy (target native,
Chicken, ZebraFinch, Mouse, Kofam symbols).
- Out: `Immune_Annotated_Final.csv`, `Unmapped_Gene_Info_Symbols.txt`

### `process_symbols.py`
Predict a gene symbol (avian consensus > Chicken > ZebraFinch > Mouse > Kofam)
and compare against the native symbol (Match, Alternative_Match, Group_Match,
Broad_Match, Paralog_Likely, Mismatch, Novel_Annotation, ...).
- Out: `Immune_Predict_Result_Final.csv`

### `merge_og_stats_to_final_gene_table.py`
Aggregate OG statistics to gene level; choose a representative OG (Primary_OG by
smallest FDR, highest enrichment, highest ref_total_genes).
- In: `..._with_OG_stats.csv`, `Immune_Predict_Result_Final.csv`, GFF3
- Out: `Immune_Predict_Result_Final_with_OG_stats.csv`,
  `Unmapped_OG_ProteinIDs.txt`

### `merge_orthology_to_final_gene_table.py`
Aggregate transcript-level orthology to gene level; add `Orth_Sub_Tier` (0-3
distinct species giving high-weight evidence: 1:1 immune, or 1:N all-immune).
- In: `..._with_OG_stats_and_orthology.csv`, `..._Final_with_OG_stats.csv`, GFF3
- Out: `Immune_Predict_Result_Final_with_OG_and_Orthology.csv`

### `visualization.py`
Seven 300-dpi PNG figures: evidence Venn, Tier distribution, top-10 immune
pathways, mouse category distribution, symbol-comparison pie, prediction-source
pie, ambiguity (7a/7b).
- In: final annotated table → Out: PNGs in `figure_dir`

---

## Final-table columns

Every column of `Immune_Predict_Result_Final_with_OG_and_Orthology.csv` is
documented in [`OUTPUT_COLUMNS.md`](OUTPUT_COLUMNS.md).
