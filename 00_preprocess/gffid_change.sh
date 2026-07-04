#!/usr/bin/env bash
# =============================================================================
# gffid_change.sh
# -----------------------------------------------------------------------------
# Normalize chromosome IDs in the raw GFF3: replace simple IDs (1, 2, Z, W, ...)
# with INSDC accessions (OY720056.1, ...) using a mapping table. Produces two
# files: one with every chromosome renamed, and one with only autosomes (the
# Z and W sex chromosomes removed).
#
# Usage:
#   gffid_change.sh -g <raw.gff3> -m <chromosome_id_map.tsv> \
#                   -r <renamed.gff3> -a <autosomes.gff3>
#
# Defaults match config/config.yaml (preprocess section). The mapping table is
# a 2-column TSV: <simple_id> <TAB> <insdc_accession>, '#' lines are comments.
# =============================================================================
set -euo pipefail

# ----- defaults (override with flags) ----------------------------------------
GFF=""
MAP="config/chromosome_id_map.tsv"
OUT_RENAMED="NM_genes.renamed.gff3"
OUT_AUTOSOMES="NM_genes.final.gff3"

usage() { grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

while getopts "g:m:r:a:h" opt; do
  case "$opt" in
    g) GFF="$OPTARG" ;;
    m) MAP="$OPTARG" ;;
    r) OUT_RENAMED="$OPTARG" ;;
    a) OUT_AUTOSOMES="$OPTARG" ;;
    h|*) usage ;;
  esac
done

[[ -n "$GFF" ]]      || { echo "ERROR: -g <raw.gff3> is required" >&2; usage; }
[[ -s "$GFF" ]]      || { echo "ERROR: GFF not found: $GFF" >&2; exit 1; }
[[ -s "$MAP" ]]      || { echo "ERROR: chromosome map not found: $MAP" >&2; exit 1; }

# Sex chromosomes are identified by their ORIGINAL simple IDs before renaming.
SEX_CHROMS="Z W"

awk -v sexlist="$SEX_CHROMS" '
BEGIN {
    FS = OFS = "\t"
    n = split(sexlist, tmp, " ")
    for (i = 1; i <= n; i++) sex[tmp[i]] = 1
}
# Load the mapping table (first file). Skip comment/blank lines.
NR == FNR {
    if ($0 ~ /^#/ || $0 == "") next
    map[$1] = $2
    next
}
# Header/comment lines go to both outputs unchanged.
/^#/ { print > out_renamed; print > out_autosomes; next }
{
    chr = $1
    is_sex = (chr in sex)
    if (chr in map) $1 = map[chr]
    print > out_renamed
    if (!is_sex) print > out_autosomes
}
' out_renamed="$OUT_RENAMED" out_autosomes="$OUT_AUTOSOMES" "$MAP" "$GFF"

echo "Done."
echo "  - all chromosomes (renamed): $OUT_RENAMED"
echo "  - autosomes only:            $OUT_AUTOSOMES"
