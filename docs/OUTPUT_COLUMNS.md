# Output column dictionary

Column reference for the final table
`Immune_Predict_Result_Final_with_OG_and_Orthology.csv` (one row per gene). The
columns are grouped below by the stage that adds them; in the file itself the
groups are interleaved (OG-stats columns are inserted just after the tier
columns, orthology columns after the OG-stats block), but the meaning of each
column is the same regardless of position.

Two families of columns depend on your configuration and inputs:

- `<Species>_*` columns are generated once **per reference species** listed in
  `species.reference` (e.g. `Chicken_GeneSymbols`, `Mouse_Description`). The
  examples below use the reference Pink Pigeon species (Mouse, Chicken,
  ZebraFinch).
- The extra functional columns copied from `gene_info.csv` depend entirely on
  the columns present in that file, so they are not enumerated here.

---

## 1. Gene identity and confidence tier
*(added by `merge_immune_annotations.py`)*

| Column | Meaning |
|--------|---------|
| `GeneID` | Gene identifier from the GFF3 (the row key). |
| `Tier` | Confidence tier from the number of independent evidence sources. **Tier 1** = all three (InterProScan + OrthoFinder + KofamScan), **Tier 2** = any two, **Tier 3** = a single source. |
| `Evidence_Sources` | Which sources support the gene, `;`-joined: `Interproscan`, `Orthofinder`, `Kofam`. |
| `Evidence_Count` | Number of supporting sources (1–3); drives `Tier`. |

## 2. Evidence detail from the three sources
*(added by `merge_immune_annotations.py`, aggregated to gene level)*

| Column | Source | Meaning |
|--------|--------|---------|
| `signature_description` | InterProScan | Member-database signature descriptions matched on the gene's proteins. |
| `interpro_annotations_description` | InterProScan | InterPro entry descriptions. |
| `Total_Score` | OrthoFinder | Score carried over from the OrthoFinder immune filtering (`Final_Filtered_List.csv`). |
| `Filter_Reason` | OrthoFinder | Why the gene passed the OrthoFinder immune filter. |
| `<Species>_GeneSymbols` | OrthoFinder / BioMart | Immune-gene symbols of that reference species' orthologs. |
| `<Species>_Description` | OrthoFinder / BioMart | Gene descriptions for those orthologs. |
| `<Species>_GO_Terms` | OrthoFinder / BioMart | GO terms for those orthologs. |
| `<Species>_Category1`, `<Species>_Subcategory`, `<Species>_UniProt_Function`, `<Species>_Immune_Source` | Curated list | Present **only** for a reference species that has a `curated_list` in the config (Mouse in the reference setup). Functional category, subcategory, UniProt function text, and whether the annotation came from the curated list or BioMart GO. |
| `KO_Gene_Symbol` | KofamScan | KEGG KO gene symbol. |
| `KO_Gene_Name` | KofamScan | KEGG KO gene name. |
| `KO` | KofamScan | Matched KEGG Orthology (KO) ID. |
| `Immune_Pathway` | KofamScan | KEGG Immune-system pathway(s) the KO belongs to. |
| `KO_definition` | KofamScan | KO definition text. |

## 3. Native symbol and functional enrichment
*(added by `add_gene_info_3.py`)*

| Column | Meaning |
|--------|---------|
| `Native_Symbol` | The gene symbol from the target genome's own annotation (from the `native_symbols` table). Used as the ground truth for symbol comparison. |
| `Matched_Gene_Symbol` | Which candidate symbol(s) matched an entry in `gene_info.csv`. Empty when enrichment is skipped or no match. |
| *(gene_info columns)* | Extra functional annotation copied from `gene_info.csv`; the exact columns depend on that file. Absent when `gene_info_enrichment: false`. |

## 4. Predicted symbol and comparison
*(added by `process_symbols.py`)*

| Column | Meaning |
|--------|---------|
| `All_Predict_Symbols` | Union of every reference-species and KO symbol available for the gene, `/`-joined (the pool the prediction draws from). |
| `Predicted_Symbol` | Final predicted symbol(s). Multiple candidates are `/`-joined. Prediction order: consensus of `consensus_group` > each species in `prediction_priority` > KofamScan KO. |
| `Predict_Sources` | Which rule produced the prediction: `Consensus`, a reference species name, or `Kofam`; an ambiguous pick carries a `(Ambiguous)` suffix. |
| `Ambiguity_Flag` | `True` when more than one candidate symbol survived (i.e. the prediction is ambiguous). |
| `Symbol_Comparison` | Category comparing `Predicted_Symbol` against `Native_Symbol`: `Match` (exact), `Alternative_Match` (matches another candidate in the pool), `Group_Match` (matches an expanded KEGG group symbol), `Broad_Match` (prefix overlap), `Paralog_Likely` (high string similarity), `Mismatch`, `Novel_Annotation` (predicted but no native symbol), `Not_Predictable` (native symbol but nothing predicted), `No_Info` (neither). |

## 5. Orthogroup enrichment statistics
*(added by `merge_og_stats_to_final_gene_table.py`; computed by `compute_og_stats.py`)*

A gene can belong to more than one orthogroup (OG). A single **representative
OG** (`Primary_OG`) is chosen per gene by smallest FDR, then highest enrichment,
then most reference genes. `Primary_*` columns describe that OG; `max_*` / `min_*`
columns summarize across all of the gene's OGs.

| Column | Meaning |
|--------|---------|
| `OG_List` | All orthogroups the gene's proteins fall into, `;`-joined. |
| `OG_Count` | Number of distinct orthogroups. |
| `Primary_OG` | The representative orthogroup ID. |
| `OG_Filter_Reason_List` | Filter reasons across the gene's OGs. |
| `Primary_immune_enrichment` | Immune fraction of the primary OG divided by the genome-wide background immune fraction (>1 means immune-enriched). |
| `Primary_FDR` | Benjamini-Hochberg FDR of a one-sided Fisher exact test for immune enrichment of the primary OG. |
| `Primary_immune_fraction` | Fraction of the reference-species genes in the primary OG that are immune. |
| `Primary_ref_total_genes` | Total reference-species genes in the primary OG. |
| `Primary_dup_count` | Number of gene-duplication events OrthoFinder reported for the OG. |
| `Primary_dup_density` | `dup_count` divided by total genes in the OG. |
| `Primary_copy_variance` | Variance of gene-copy number across all species in the OG (a paralog-expansion signal). |
| `Primary_immune_entropy` | Normalized Shannon entropy of immune-gene counts across reference species (how evenly the immune signal is spread). |
| `Primary_immune_<Species>` | Immune-gene count contributed by each reference species in the primary OG. |
| `max_immune_enrichment` | Highest `immune_enrichment` across all of the gene's OGs. |
| `min_FDR` | Smallest FDR across all of the gene's OGs. |
| `max_dup_density` | Highest duplication density across the gene's OGs. |
| `max_copy_variance` | Highest copy-number variance across the gene's OGs. |

## 6. Orthology relationships
*(added by `merge_orthology_to_final_gene_table.py`; computed by `compute_orthology_stats.py`)*

Per gene, the best relationship is kept per reference species (type priority
`1:1 > 1:N > N:1 > N:M`, counts maxed across transcripts).

| Column | Meaning |
|--------|---------|
| `Gene_Transcript_Count` | Number of transcripts aggregated for this gene. |
| `Orth_<Species>_type` | Orthology relationship of the gene to that reference species: `1:1`, `1:N`, `N:1`, or `N:M`. |
| `Orth_<Species>_ref_total` | Number of that species' orthologs in the best-supporting orthogroup. |
| `Orth_<Species>_ref_immune` | How many of those orthologs are annotated immune. |
| `Orth_<Species>_all_immune` | `True` when every ortholog of that species is immune. |
| `Orth_<Species>_rows` | Number of supporting orthogroup rows for that species. |
| `Orth_1to1_immune_count` | Number of reference species where the gene has a 1:1 ortholog that is immune. |
| `Orth_1to1_immune_species` | Which species those are, `;`-joined. |
| `Orth_1toN_allimmune_count` | Number of reference species where the gene has a 1:N relationship and **all** orthologs are immune. |
| `Orth_1toN_allimmune_species` | Which species those are, `;`-joined. |
| `Orth_Sub_Tier` | Count of distinct reference species giving high-weight orthology evidence (a 1:1 immune ortholog, or a 1:N all-immune relationship). Range 0–3 for three references; blank when the gene has no OrthoFinder evidence. |
