#!/bin/bash
# Build the c4 tool env on CSD3 (login node; needs internet).  Idempotent.
#   bash c4_imputation/env/setup_csd3_env.sh
# Installs to $C4_ENV (default ~/rds/hpc-work/envs/c4 -- never /home, it is full).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"
C4_ENV="${C4_ENV:-$HOME/rds/hpc-work/envs/c4}"
if [[ -x "$C4_ENV/bin/bcftools" ]]; then echo "c4 env present: $C4_ENV"; exit 0; fi
MM=""
for c in "$(command -v micromamba || true)" "$REPO/legacy/hpc/work/bin/micromamba" "$REPO/legacy/hpc/work/mamba/bin/micromamba" \
         "$(command -v mamba || true)" "$(command -v conda || true)"; do
  [[ -n "$c" && -x "$c" ]] && { MM="$c"; break; }
done
[[ -n "$MM" ]] || { echo "no micromamba/mamba/conda found; 'module load miniconda/3' or install micromamba, then rerun" >&2; exit 2; }
echo "using $MM -> $C4_ENV"
"$MM" create -y -p "$C4_ENV" -f "$HERE/c4_env.yml"
JAVA="$C4_ENV/bin/java"; [[ -x "$JAVA" ]] || JAVA="$C4_ENV/lib/jvm/bin/java"
"$C4_ENV/bin/bcftools" --version | head -1; "$JAVA" -version 2>&1 | head -1
