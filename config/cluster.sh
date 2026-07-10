# =============================================================================
# config/cluster.sh - shell-side config for the HPC SLURM scripts
# -----------------------------------------------------------------------------
# The Python stages read config/config.yaml. The three SLURM wrapper scripts
# (interproscan_run.sh, run_kofam.sh, run_orthofind.sh) cannot cleanly parse
# YAML, so they `source` THIS file instead. Edit it once per project; you do not
# need to edit the SLURM scripts themselves (except the #SBATCH header lines,
# which SLURM requires to be literal: partition, memory, email).
#
# The scripts locate this file via $SLURM_SUBMIT_DIR, so submit them from the
# repo root (e.g. `sbatch 01_interproscan/interproscan_run.sh`), or set the
# CLUSTER_CONFIG environment variable to its absolute path.
#
# Keep the paths below in sync with config.yaml: the per-stage *_WORK_DIR values
# must match the corresponding work_dir in config.yaml so the SLURM raw outputs
# land where the Python stages look for them.
# =============================================================================

# Workspace root: where the INPUT files live (genome, GFF, protein FASTA, etc.).
WORKSPACE="/path/to/your/workspace"

# Target-species protein FASTA (produced by 00_preprocess/extract_proteins.sh).
# Absolute path; keep in sync with reference.protein_fasta in config.yaml.
PROTEIN_FASTA_PATH="${WORKSPACE}/PinkPigeon.faa"

# --- Per-stage output directories (match work_dir in config.yaml) ------------
IPR_WORK_DIR="${WORKSPACE}/outputs/01_interproscan"     # ipr_chunks/, ipr_out/
KOFAM_WORK_DIR="${WORKSPACE}/outputs/02_kofamscan"      # result_kofam_detail.txt etc.

# --- Conda (used by run_kofam.sh) --------------------------------------------
CONDA_SH="/maps/usermodules/shared_software/shared_envmodules/conda/conda-25.1.1/etc/profile.d/conda.sh"
CONDA_ENV="immune_annotation"

# --- KofamScan database (run_kofam.sh) ---------------------------------------
# config.yml must point profile:/ko_list: at your KofamScan DB (see README).
KOFAM_DB_CONFIG="${WORKSPACE}/database/kofam_db/config.yml"

# --- OrthoFinder (run_orthofind.sh) ------------------------------------------
# Must contain exactly the target + reference protein FASTAs, each named to
# match its species name in config.yaml (e.g. PinkPigeon.faa, Mouse.faa, ...).
ORTHO_INPUT_DIR="${WORKSPACE}/ortho_input"
# Fixed results-folder name -> Results_<ORTHO_RUN_NAME> (no date), so the
# OrthoFinder paths in config.yaml stay valid across re-runs.
ORTHO_RUN_NAME="pinkpigeon"
