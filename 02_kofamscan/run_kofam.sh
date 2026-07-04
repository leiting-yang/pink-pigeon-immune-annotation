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

# ----- USER CONFIG (keep in sync with config/config.yaml) --------------------
CONDA_SH="/maps/usermodules/shared_software/shared_envmodules/conda/conda-25.1.1/etc/profile.d/conda.sh"
CONDA_ENV="pigeon_immune"
WORK_DIR="/path/to/your/workspace"
INPUT_FASTA="${WORK_DIR}/PinkPigeon.faa"
CONFIG_FILE="${WORK_DIR}/database/kofam_db/config.yml"
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
