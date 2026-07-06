# =============================================================================
# TEST cluster config - Mauritius kestrel (bFalPun1.1)
# -----------------------------------------------------------------------------
# Used by the SLURM scripts for the kestrel test. Because the scripts default to
# config/cluster.sh, point them here by exporting CLUSTER_CONFIG before sbatch:
#
#   cd /maps/projects/echo/people/blv970/immune_test/pink-pigeon-immune-annotation
#   export CLUSTER_CONFIG="$PWD/config/cluster_kestrel.sh"
#   sbatch 01_interproscan/interproscan_run.sh
#
# (sbatch inherits your shell environment, so the export carries into the job.)
# =============================================================================

# Target-species protein FASTA (produced by extract_proteins.sh in the workspace)
PROTEIN_FASTA="Kestrel.faa"

# --- Conda (run_kofam.sh) ----------------------------------------------------
CONDA_SH="/maps/usermodules/shared_software/shared_envmodules/conda/conda-25.1.1/etc/profile.d/conda.sh"
CONDA_ENV="immune_annotation"

# --- KofamScan (run_kofam.sh) ------------------------------------------------
KOFAM_WORK_DIR="/maps/projects/echo/people/blv970/immune_test/pink-pigeon-immune-annotation"
# KofamScan HMM database config.yml (species-independent - reuse your existing one).
# >>> EDIT: point this at the KofamScan DB config you already have. <<<
KOFAM_DB_CONFIG="/path/to/your/kofam_db/config.yml"

# --- OrthoFinder (run_orthofind.sh) ------------------------------------------
# Must contain exactly 4 protein FASTAs named to match config species:
#   Kestrel.faa, Mouse.faa, Chicken.faa, ZebraFinch.faa
ORTHO_INPUT_DIR="/maps/projects/echo/people/blv970/immune_test/pink-pigeon-immune-annotation/ortho_input"
