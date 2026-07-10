#!/usr/bin/env bash
# =============================================================================
# extract_native_symbols.sh   (OPTIONAL preprocessing step, ON by default)
# -----------------------------------------------------------------------------
# Build the target-species "native symbol" table straight from the annotation
# GFF3, so downstream symbol comparison (04_tiering/process_symbols.py) has the
# annotation's own gene symbols to compare the predicted ones against.
#
# Output: 2-column TSV  <GeneSymbol> <TAB> <GeneID>  matching what
# external_data.native_symbols expects. The GeneID is the gene feature's ID
# attribute (e.g. NCBI RefSeq "gene-LOC..." / "gene-ACYX44_..."), which is the
# same gene key the pipeline uses everywhere else, so the symbols attach correctly.
#
# The symbol is taken from the GFF attributes, preferring `gene=` then `Name=`.
# This matches NCBI RefSeq and Ensembl GFF3 conventions. Uncharacterized genes
# keep their placeholder (e.g. LOCxxxxx), which is the honest native annotation.
#
# WHEN TO SKIP: if you already have a curated native symbol table, or your GFF
# has no gene symbols, just point external_data.native_symbols at your own file
# (or an empty file) and do not run this script.
#
# Usage:
#   extract_native_symbols.sh -g <annotation.gff3> -o <Species_genes_by_name.tsv>
# =============================================================================
set -euo pipefail

GFF=""
OUT=""

usage() { grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

while getopts "g:o:h" opt; do
  case "$opt" in
    g) GFF="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    h|*) usage ;;
  esac
done

[[ -n "$GFF" ]] || { echo "ERROR: -g <annotation.gff3> is required" >&2; usage; }
[[ -n "$OUT" ]] || { echo "ERROR: -o <output.tsv> is required" >&2; usage; }
[[ -s "$GFF" ]] || { echo "ERROR: GFF3 not found: $GFF" >&2; exit 1; }

# gene / pseudogene feature lines carry ID=<gene id> and gene=/Name=<symbol>.
awk -F'\t' '$3=="gene" || $3=="pseudogene" {
  id=""; sym="";
  n=split($9, a, ";");
  for (i=1; i<=n; i++) {
    if      (a[i] ~ /^ID=/)                id  = substr(a[i], 4);
    else if (a[i] ~ /^gene=/)              sym = substr(a[i], 6);
    else if (a[i] ~ /^Name=/ && sym == "") sym = substr(a[i], 6);
  }
  if (id != "" && sym != "") print sym "\t" id;
}' "$GFF" > "$OUT"

echo "Done. Native symbol table written to: $OUT"
echo "Rows: $(wc -l < "$OUT")"
