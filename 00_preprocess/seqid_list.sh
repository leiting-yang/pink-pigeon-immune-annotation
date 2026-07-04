#!/usr/bin/env bash
# =============================================================================
# seqid_list.sh  (MANUAL DIAGNOSTIC - not part of the automated flow)
# -----------------------------------------------------------------------------
# Compare the sequence IDs in the genome FASTA against the sequence IDs in the
# GFF3 to spot ID mismatches before running any annotation. The output is meant
# to be inspected by a human; nothing downstream consumes it.
#
# Usage:
#   seqid_list.sh -f <genome.fasta> -g <annotation.gff3> [-o ids_compare.csv]
# =============================================================================
set -euo pipefail

FASTA=""
GFF=""
OUT="ids_compare.csv"

usage() { grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

while getopts "f:g:o:h" opt; do
  case "$opt" in
    f) FASTA="$OPTARG" ;;
    g) GFF="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    h|*) usage ;;
  esac
done

[[ -s "$FASTA" ]] || { echo "ERROR: FASTA not found: $FASTA" >&2; exit 1; }
[[ -s "$GFF" ]]   || { echo "ERROR: GFF not found: $GFF" >&2; exit 1; }

tmp_fasta="$(mktemp)"
tmp_gff="$(mktemp)"
trap 'rm -f "$tmp_fasta" "$tmp_gff"' EXIT

# Sequence IDs from the FASTA headers (first whitespace-delimited token)
grep "^>" "$FASTA" | cut -d " " -f1 | sed 's/^>//' | sort -u > "$tmp_fasta"

# Distinct sequence IDs from GFF column 1 (skip comments)
grep -v "^#" "$GFF" | cut -f1 | sort -u > "$tmp_gff"

{
  echo "FASTA_ID,GFF3_ID"
  paste -d "," "$tmp_fasta" "$tmp_gff"
} > "$OUT"

echo "Wrote comparison table: $OUT"
echo "Inspect it manually to confirm FASTA and GFF3 sequence IDs are consistent."
