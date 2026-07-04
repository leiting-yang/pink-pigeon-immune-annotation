#!/usr/bin/env bash
#SBATCH --job-name=ortho_pigeon
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

# ----- USER CONFIG (keep in sync with config/config.yaml) --------------------
INPUT_DIR="/path/to/your/workspace/orthofinder/ortho_input"
# -----------------------------------------------------------------------------

module load orthofinder/3.0.1b1

echo "Running OrthoFinder on: ${INPUT_DIR}"
echo "Threads: ${SLURM_CPUS_PER_TASK}"

# -f: input folder (4 species FASTAs)
# -t: parallel sequence-search threads
# -a: parallel analysis threads
orthofinder \
    -f "${INPUT_DIR}" \
    -t "${SLURM_CPUS_PER_TASK}" \
    -a "${SLURM_CPUS_PER_TASK}"

echo "Done. Results are under ${INPUT_DIR}/OrthoFinder/Results_*"
