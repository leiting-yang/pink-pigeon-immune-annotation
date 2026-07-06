#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh
# -----------------------------------------------------------------------------
# Reference driver documenting the full run order of the immune-annotation
# pipeline. It runs the local Python "glue" stages directly, and STOPS at each
# heavy HPC step (InterProScan / KofamScan / OrthoFinder), which must be
# submitted with sbatch and allowed to finish before continuing.
#
# It is deliberately conservative: run it stage by stage rather than end to end
# the first time. Every Python step also accepts --config, so you can run any
# single command by hand.
#
# Usage:
#   ./run_pipeline.sh <stage>
#   stages: preprocess | interproscan | kofamscan | orthofinder | tiering | all
# =============================================================================
set -euo pipefail

CONFIG="config/config.yaml"
PY="python"

banner() { echo; echo "=============== $* ==============="; }
manual() { echo; echo ">>> MANUAL / HPC STEP: $*"; echo ">>> Submit it, wait for completion, then continue."; }

# -----------------------------------------------------------------------------
stage_preprocess() {
  banner "Stage 00 - preprocess"
  bash 00_preprocess/extract_proteins.sh -g "<genome.fasta>" -a "<raw.gff3>" -o PinkPigeon.faa
  # OPTIONAL: only needed if the genome FASTA and GFF3 use different chromosome
  # IDs (e.g. Ensembl simple IDs). Skip it if they already match (most NCBI
  # RefSeq assemblies). The core pipeline uses the raw GFF3 regardless.
  # bash 00_preprocess/gffid_change.sh -g "<raw.gff3>" -m config/chromosome_id_map.tsv
  echo "(optional QC) bash 00_preprocess/seqid_list.sh -f <genome.fasta> -g <raw.gff3>"
}

stage_interproscan() {
  banner "Stage 01 - InterProScan"
  manual "seqkit split the protein FASTA, then: sbatch 01_interproscan/interproscan_run.sh"
  manual "merge chunks: sbatch --array=9 01_interproscan/interproscan_run.sh"
  $PY 01_interproscan/immunewash.py --config "$CONFIG"
}

stage_kofamscan() {
  banner "Stage 02 - KofamScan"
  $PY 02_kofamscan/prepare_db.py --config "$CONFIG"
  manual "sbatch 02_kofamscan/run_kofam.sh"
  $PY 02_kofamscan/extract_map.py     --config "$CONFIG"
  $PY 02_kofamscan/filter_kofam.py    --config "$CONFIG"
  $PY 02_kofamscan/map_to_gene.py     --config "$CONFIG"
  $PY 02_kofamscan/add_ko_metadata.py --config "$CONFIG"
}

stage_orthofinder() {
  banner "Stage 03 - OrthoFinder"
  $PY 03_orthofinder/data_merge_new.py --config "$CONFIG"
  manual "sbatch 03_orthofinder/run_orthofind.sh"
  $PY 03_orthofinder/final_filtering_v2.py       --config "$CONFIG"
  $PY 03_orthofinder/compute_og_stats.py         --config "$CONFIG"
  $PY 03_orthofinder/compute_orthology_stats.py  --config "$CONFIG"
}

stage_tiering() {
  banner "Stage 04 - tiering"
  $PY 04_tiering/merge_immune_annotations.py          --config "$CONFIG"
  $PY 04_tiering/add_gene_info_3.py                   --config "$CONFIG"
  $PY 04_tiering/process_ppg_genes.py                 --config "$CONFIG"
  $PY 04_tiering/merge_og_stats_to_final_gene_table.py --config "$CONFIG"
  $PY 04_tiering/merge_orthology_to_final_gene_table.py --config "$CONFIG"
  $PY 04_tiering/visualization.py                     --config "$CONFIG"
}

case "${1:-}" in
  preprocess)   stage_preprocess ;;
  interproscan) stage_interproscan ;;
  kofamscan)    stage_kofamscan ;;
  orthofinder)  stage_orthofinder ;;
  tiering)      stage_tiering ;;
  all)
    stage_preprocess; stage_interproscan; stage_kofamscan
    stage_orthofinder; stage_tiering ;;
  *)
    echo "Usage: $0 {preprocess|interproscan|kofamscan|orthofinder|tiering|all}"
    exit 1 ;;
esac
