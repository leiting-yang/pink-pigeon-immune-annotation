#!/usr/bin/env bash
# =============================================================================
# check_deps.sh - verify every external dependency for the Kestrel test run.
# Run on the cluster (login node is fine):  bash test/check_deps.sh
#
# It does NOT change anything. For each tool it starts from a clean module
# environment, loads the module line you specify below (edit the *_MODS vars if
# a name/prereq differs on your cluster), and checks the command resolves.
# =============================================================================

# --- EDIT THESE if `module spider <tool>` shows different names/prereqs -------
JAVA_MODS="openjdk"
IPR_MODS="openjdk interproscan/5.73-104.0"
SEQKIT_MODS="seqkit"
GFFREAD_MODS="gffread"
ORTHO_MODS="orthofinder/3.0.1b1"
KOFAM_MODS=""                     # set if KofamScan is a module, else leave empty
CONDA_SH="/maps/usermodules/shared_software/shared_envmodules/conda/conda-25.1.1/etc/profile.d/conda.sh"
CONDA_ENV="immune_annotation"
KOFAM_DB_CONFIG="/maps/projects/echo/people/blv970/kofam_db/config.yml"
# -----------------------------------------------------------------------------

ok()   { printf '  \033[32mOK\033[0m   %s\n' "$1"; }
bad()  { printf '  \033[31mFAIL\033[0m %s\n' "$1"; }

check_module_tool() {   # name, "module load spec", command-to-find
  local name="$1" mods="$2" cmd="$3"
  ( module purge >/dev/null 2>&1
    if [ -n "$mods" ]; then
      # shellcheck disable=SC2086
      module load $mods >/dev/null 2>&1 || { bad "$name: 'module load $mods' failed"; exit 1; }
    fi
    if command -v "$cmd" >/dev/null 2>&1; then
      ok "$name -> $(command -v "$cmd")"
    else
      bad "$name: '$cmd' not on PATH after 'module load $mods'"
    fi )
}

echo "=== 1. Conda env Python libs ==="
# Purge modules first: a loaded tool module (e.g. orthofinder) prepends its own
# python to PATH and shadows the env's python. Python steps must run clean.
module purge >/dev/null 2>&1
if [ -f "$CONDA_SH" ]; then
  # shellcheck disable=SC1090
  source "$CONDA_SH"
  conda activate "$CONDA_ENV" 2>/dev/null \
    && ok "conda env '$CONDA_ENV' activated" \
    || bad "cannot activate conda env '$CONDA_ENV'"
  py=$(command -v python)
  case "$py" in
    *"/envs/$CONDA_ENV/"*) ok "python -> $py" ;;
    *) bad "python is NOT the env python: $py (activation did not take)" ;;
  esac
  for m in pandas numpy yaml requests matplotlib seaborn matplotlib_venn; do
    python -c "import $m" 2>/dev/null && ok "python: $m" || bad "python: $m missing"
  done
else
  bad "conda.sh not found: $CONDA_SH"
fi

echo "=== 2. Heavy tools via modules ==="
check_module_tool "gffread"      "$GFFREAD_MODS" "gffread"
check_module_tool "seqkit"       "$SEQKIT_MODS"  "seqkit"
check_module_tool "java(openjdk)" "$JAVA_MODS"   "java"
check_module_tool "interproscan" "$IPR_MODS"     "interproscan.sh"
check_module_tool "orthofinder"  "$ORTHO_MODS"   "orthofinder"

echo "=== 3. KofamScan (exec_annotation) ==="
( module purge >/dev/null 2>&1
  [ -n "$KOFAM_MODS" ] && { module load $KOFAM_MODS >/dev/null 2>&1 || bad "module load $KOFAM_MODS failed"; }
  # also try the conda env, in case it was installed there
  [ -f "$CONDA_SH" ] && { source "$CONDA_SH"; conda activate "$CONDA_ENV" 2>/dev/null; }
  if command -v exec_annotation >/dev/null 2>&1; then
    ok "exec_annotation -> $(command -v exec_annotation)"
  else
    bad "exec_annotation NOT found (install kofamscan or load its module)"
  fi )

echo "=== 4. KofamScan database ==="
if [ -f "$KOFAM_DB_CONFIG" ]; then
  ok "config.yml present: $KOFAM_DB_CONFIG"
  prof=$(awk -F': *' '/^profile/{print $2}' "$KOFAM_DB_CONFIG")
  kol=$(awk  -F': *' '/^ko_list/{print $2}' "$KOFAM_DB_CONFIG")
  [ -e "$prof" ] && ok "profiles: $prof" || bad "profiles missing: $prof"
  [ -e "$kol" ]  && ok "ko_list:  $kol"  || bad "ko_list missing: $kol"
else
  bad "KOFAM_DB_CONFIG missing: $KOFAM_DB_CONFIG"
fi

echo "=== 5. KEGG hierarchy JSON (for filter_kofam.py) ==="
JSON="/maps/projects/echo/people/blv970/immune_test/ko00001.json"
[ -s "$JSON" ] && ok "ko00001.json present: $JSON" || bad "ko00001.json missing: $JSON"

echo
echo "Done. Fix every FAIL before running the pipeline."
