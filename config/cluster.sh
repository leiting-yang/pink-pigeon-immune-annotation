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
# =============================================================================

# Target-species protein FASTA (produced by 00_preprocess/extract_proteins.sh:
#   gffread -y PinkPigeon.faa -g <genome_fasta> <gff3_raw>).
# Keep this in sync with reference.protein_fasta in config.yaml.
PROTEIN_FASTA="PinkPigeon.faa"

# --- Conda (used by run_kofam.sh) --------------------------------------------
CONDA_SH="/maps/usermodules/shared_software/shared_envmodules/conda/conda-25.1.1/etc/profile.d/conda.sh"
CONDA_ENV="immune_annotation"

# --- KofamScan (run_kofam.sh) ------------------------------------------------
KOFAM_WORK_DIR="/path/to/your/workspace"
KOFAM_DB_CONFIG="${KOFAM_WORK_DIR}/database/kofam_db/config.yml"

# --- OrthoFinder (run_orthofind.sh) ------------------------------------------
# Must contain exactly the target + reference protein FASTAs, each named to
# match its species name in config.yaml (e.g. PinkPigeon.faa, Mouse.faa, ...).
ORTHO_INPUT_DIR="/path/to/your/workspace/orthofinder/ortho_input"
