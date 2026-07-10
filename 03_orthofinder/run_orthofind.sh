#!/usr/bin/env bash
#SBATCH --job-name=orthofinder
#SBATCH --output=ortho.%J.out
#SBATCH --error=ortho.%J.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00
#SBATCH --partition=cpuqueue
#SBATCH --mail-type=fail,end
#SBATCH --mail-user=your.email@example.com
# =============================================================================
# run_orthofind.sh
# -----------------------------------------------------------------------------
# Run OrthoFinder on the four-species protein FASTA set (PinkPigeon, Mouse,
# Chicken, ZebraFinch). The input directory must contain exactly those four
# protein FASTA files and nothing else.
#
# Usage: sbatch run_orthofind.sh
# =============================================================================
set -euo pipefail

# ----- Shared cluster config (edit config/cluster.sh, not this script) -------
CLUSTER_CONFIG="${CLUSTER_CONFIG:-${SLURM_SUBMIT_DIR:-.}/config/cluster.sh}"
[ -f "$CLUSTER_CONFIG" ] || { echo "ERROR: cluster config not found: $CLUSTER_CONFIG (submit from the repo root or set CLUSTER_CONFIG)" >&2; exit 1; }
source "$CLUSTER_CONFIG"
INPUT_DIR="$ORTHO_INPUT_DIR"
# -----------------------------------------------------------------------------

module load orthofinder/3.0.1b1

# Deterministic results folder: -n makes OrthoFinder write to
# Results_${ORTHO_RUN_NAME} instead of the default date-stamped Results_<date>,
# so the config paths never need editing between runs. Default the name to
# "run" if cluster.sh does not set ORTHO_RUN_NAME.
RUN_NAME="${ORTHO_RUN_NAME:-run}"

# Remove a prior run of the same name first; otherwise OrthoFinder appends a
# _1/_2 suffix and the folder name stops being deterministic.
rm -rf "${INPUT_DIR}/OrthoFinder/Results_${RUN_NAME}"

echo "Running OrthoFinder on: ${INPUT_DIR} (results: Results_${RUN_NAME})"
echo "Threads: ${SLURM_CPUS_PER_TASK}"

# -f: input folder (4 species FASTAs)
# -n: fixed results-folder name suffix
# -t: parallel sequence-search threads
# -a: parallel analysis threads
orthofinder \
    -f "${INPUT_DIR}" \
    -n "${RUN_NAME}" \
    -t "${SLURM_CPUS_PER_TASK}" \
    -a "${SLURM_CPUS_PER_TASK}"

echo "Done. Results are under ${INPUT_DIR}/OrthoFinder/Results_${RUN_NAME}"
