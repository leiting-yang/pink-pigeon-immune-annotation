# =============================================================================
# test/cluster.sh - shell-side config for the HPC SLURM scripts (Kestrel test)
# -----------------------------------------------------------------------------
# The three SLURM wrapper scripts source THIS file. Because it lives in test/
# (not config/), submit the jobs with CLUSTER_CONFIG pointing here, e.g.:
#
#   export CLUSTER_CONFIG=/maps/.../immune-annotation-pipeline/test/cluster.sh
#   sbatch --export=ALL,CLUSTER_CONFIG="$CLUSTER_CONFIG" \
#          /maps/.../immune-annotation-pipeline/01_interproscan/interproscan_run.sh
# =============================================================================

# Workspace root: where the INPUT files live (genome, GFF, protein FASTA, etc.).
WORKSPACE="/maps/projects/echo/people/blv970/immune_test"

# Target-species protein FASTA (produced by 00_preprocess/extract_proteins.sh).
# Absolute path so both InterProScan and KofamScan find it regardless of their
# per-stage output directories.
PROTEIN_FASTA_PATH="${WORKSPACE}/Kestrel.faa"

# --- Per-stage output directories (keep in sync with work_dir in config.yaml) -
IPR_WORK_DIR="${WORKSPACE}/outputs/01_interproscan"     # ipr_chunks/, ipr_out/
KOFAM_WORK_DIR="${WORKSPACE}/outputs/02_kofamscan"      # result_kofam_detail.txt etc.

# --- Conda (used by run_kofam.sh) --------------------------------------------
CONDA_SH="/maps/usermodules/shared_software/shared_envmodules/conda/conda-25.1.1/etc/profile.d/conda.sh"
CONDA_ENV="immune_annotation"

# --- KofamScan database (run_kofam.sh) ---------------------------------------
# config.yml must point profile:/ko_list: at THIS directory (not the old block2 path).
KOFAM_DB_CONFIG="/maps/projects/echo/people/blv970/kofam_db/config.yml"

# --- OrthoFinder (run_orthofind.sh) ------------------------------------------
# Dir holding EXACTLY the 4 protein FASTAs; OrthoFinder writes its Results here.
ORTHO_INPUT_DIR="${WORKSPACE}/ortho_input"
# Fixed results-folder name -> Results_<ORTHO_RUN_NAME> (no date), so the
# OrthoFinder paths in config.yaml stay valid across re-runs.
ORTHO_RUN_NAME="kestrel"
