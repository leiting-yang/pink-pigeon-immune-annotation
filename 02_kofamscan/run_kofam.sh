#!/usr/bin/env bash
#SBATCH --job-name=kofam_run
#SBATCH --output=kofam_run.%J.out
#SBATCH --error=kofam_run.%J.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --partition=cpuqueue
#SBATCH --mail-type=fail,end
#SBATCH --mail-user=your.email@example.com
# =============================================================================
# run_kofam.sh
# -----------------------------------------------------------------------------
# Run KofamScan (HMMER-based KEGG KO assignment) on the Pink Pigeon protein
# FASTA, emitting the detail-tsv format (keeps the score/threshold columns that
# filter_kofam.py needs).
#
# Usage: sbatch run_kofam.sh
# Edit the USER CONFIG block to match your environment / config.yaml.
# =============================================================================
# Do not use `set -u` here: it breaks conda activation.
set -e
set -o pipefail

# ----- Shared cluster config (edit config/cluster.sh, not this script) -------
CLUSTER_CONFIG="${CLUSTER_CONFIG:-${SLURM_SUBMIT_DIR:-.}/config/cluster.sh}"
[ -f "$CLUSTER_CONFIG" ] || { echo "ERROR: cluster config not found: $CLUSTER_CONFIG (submit from the repo root or set CLUSTER_CONFIG)" >&2; exit 1; }
source "$CLUSTER_CONFIG"

# ----- Script-local settings -------------------------------------------------
WORK_DIR="$KOFAM_WORK_DIR"
INPUT_FASTA="${WORK_DIR}/${PROTEIN_FASTA}"
CONFIG_FILE="$KOFAM_DB_CONFIG"
OUTPUT_FILE="result_kofam_detail.txt"
# -----------------------------------------------------------------------------

# 1. Environment
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

# 2. Work directory
cd "${WORK_DIR}"

# 3. Remove stale output
[ -f "${OUTPUT_FILE}" ] && { echo "Removing old output: ${OUTPUT_FILE}"; rm -f "${OUTPUT_FILE}"; }

# 4. Private temp dir
export TMPDIR="${WORK_DIR}/tmp_kofam_$$"
mkdir -p "${TMPDIR}"

# 5. Environment check
echo "=== Environment Check ==="
echo "Date:            $(date)"
echo "Host:            $(hostname)"
echo "Work dir:        $(pwd)"
echo "Input FASTA:     ${INPUT_FASTA}"
echo "Sequences:       $(grep -c '^>' "${INPUT_FASTA}")"
echo "CPUs:            ${SLURM_CPUS_PER_TASK}"
echo "Config:          ${CONFIG_FILE}"
echo "Temp dir:        ${TMPDIR}"
echo "========================="

# 6. Run KofamScan (detail-tsv keeps thresholds and E-values)
echo "Starting KofamScan at $(date) with ${SLURM_CPUS_PER_TASK} threads..."
exec_annotation \
    -f detail-tsv \
    -c "${CONFIG_FILE}" \
    --cpu "${SLURM_CPUS_PER_TASK}" \
    --tmp-dir "${TMPDIR}" \
    "${INPUT_FASTA}" \
    -o "${OUTPUT_FILE}" 2>&1 | tee kofam_run.log

echo "Finished at $(date). Results in ${OUTPUT_FILE}"

# 7. Cleanup
rm -rf "${TMPDIR}"
echo "Done."
