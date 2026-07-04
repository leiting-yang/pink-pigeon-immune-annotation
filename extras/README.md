# extras/

Scripts kept for reference but **not part of the core annotation pipeline**.
They are preserved as-is (original code, not refactored).

## `immune_bed_export/`

A side branch that exports immune genes to BED format for a separate
genomic-coordinate analysis (window / ROH overlap, etc.). It is **not** an input
to the tiering stage.

- `extract_immune_mrna.sh` - pull the mRNA lines of immune proteins from the
  renamed GFF3, producing `immune_genes.gff` and a transcript-gene id table.
- `immunebed.py` - convert `immune_genes.gff` to a 6-column BED, keeping the
  longest transcript per gene.

Run order (if you need the BED): `extract_immune_mrna.sh` then `immunebed.py`.
Both consume the output of stage 01 (`interproscan_immune_results.csv`) and the
renamed GFF3 from stage 00.

## `test_pipeline.py`

A standalone self-test of the KofamScan filtering logic (threshold filter →
immune-KO filter → E-value deduplication) on tiny mock data. Useful for
sanity-checking `02_kofamscan/filter_kofam.py`, not for production runs.
