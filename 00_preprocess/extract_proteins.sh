#!/usr/bin/env bash
# =============================================================================
# extract_proteins.sh
# -----------------------------------------------------------------------------
# Extract the target-species protein FASTA from the genome FASTA and GFF3 with
# gffread. The output (default PinkPigeon.faa) is the input to BOTH InterProScan
# (stage 01) and KofamScan (stage 02).
#
# Usage:
#   extract_proteins.sh -g <genome.fasta> -a <annotation.gff3> [-o PinkPigeon.faa]
# =============================================================================
set -euo pipefail

GENOME=""
GFF=""
OUT="PinkPigeon.faa"

usage() { grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

while getopts "g:a:o:h" opt; do
  case "$opt" in
    g) GENOME="$OPTARG" ;;
    a) GFF="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    h|*) usage ;;
  esac
done

[[ -n "$GENOME" ]] || { echo "ERROR: -g <genome.fasta> is required" >&2; usage; }
[[ -n "$GFF" ]]    || { echo "ERROR: -a <annotation.gff3> is required" >&2; usage; }
[[ -s "$GENOME" ]] || { echo "ERROR: genome FASTA not found: $GENOME" >&2; exit 1; }
[[ -s "$GFF" ]]    || { echo "ERROR: GFF3 not found: $GFF" >&2; exit 1; }

module load gffread

# -y : write the protein translation of each transcript
# -g : the genome FASTA to pull sequences from
gffread -y "$OUT" -g "$GENOME" "$GFF"

echo "Done. Protein FASTA written to: $OUT"
echo "Sequences: $(grep -c '^>' "$OUT")"
