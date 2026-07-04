#!/usr/bin/env bash
#SBATCH --job-name=iprscan
#SBATCH --output=iprscan.%A_%a.out
#SBATCH --error=iprscan.%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=96:00:00
#SBATCH --partition=cpuqueue
#SBATCH --array=1-8%8
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@example.com
# =============================================================================
# interproscan_run.sh
# -----------------------------------------------------------------------------
# Run InterProScan on the Pink Pigeon protein FASTA as a SLURM array job.
#
# WORKFLOW (three steps, run in order):
#   1. MANUAL, once, before submitting the array:
#        module load seqkit
#        seqkit split -p ${NCHUNK} -O ${CHUNK_DIR} ${INFA}
#   2. Array job (this script, array 1-N): one InterProScan run per chunk.
#        sbatch interproscan_run.sh
#   3. Merge (this script with a single index past the array, e.g. task 9):
#        sbatch --array=9 interproscan_run.sh
#
# Adjust NCHUNK and --array above together (they must match the chunk count).
# =============================================================================
set -euo pipefail

# ----- USER CONFIG (keep in sync with config/config.yaml) --------------------
INFA="CDS_protein_seq.fa"          # Pink Pigeon CDS protein FASTA
NCHUNK=8                            # number of chunks (must match --array 1-N)
CHUNK_DIR="$PWD/ipr_chunks"
OUTDIR="$PWD/ipr_out"
TMPDIR="${SLURM_TMPDIR:-$PWD/ipr_tmp}"
MERGE_TASK_ID=9                    # array index reserved for the merge step
# -----------------------------------------------------------------------------

mkdir -p "$CHUNK_DIR" "$OUTDIR" "$TMPDIR"

module load openjdk
module load interproscan/5.73-104.0
module load seqkit

# ----- Merge step ------------------------------------------------------------
if [[ "${SLURM_ARRAY_TASK_ID}" -eq "${MERGE_TASK_ID}" ]]; then
  echo "[INFO] Merging per-chunk InterProScan TSVs..."
  # InterProScan TSV output has no header, so plain concatenation is correct.
  cat "${OUTDIR}"/*.tsv > "${OUTDIR}/interproscan_merged.tsv"
  echo "[DONE] Merged into ${OUTDIR}/interproscan_merged.tsv"
  exit 0
fi

# ----- Scan step (one chunk per array task) ----------------------------------
if [[ "${SLURM_ARRAY_TASK_ID}" -ge 1 && "${SLURM_ARRAY_TASK_ID}" -le "${NCHUNK}" ]]; then
  CHUNK=$(ls "${CHUNK_DIR}"/* | sed -n "${SLURM_ARRAY_TASK_ID}p")
  [ -s "${CHUNK}" ] || { echo "ERROR: chunk not found: ${CHUNK}" >&2; exit 1; }

  interproscan.sh \
    -i "${CHUNK}" \
    -f tsv \
    -dp \
    -goterms \
    --output-dir "${OUTDIR}" \
    --tempdir "${TMPDIR}"

  echo "[ARRAY ${SLURM_ARRAY_TASK_ID}] done."
  exit 0
fi

echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} is out of range (1-${NCHUNK} scan, ${MERGE_TASK_ID} merge)." >&2
exit 1
