#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
ENV_ROOT="$ROOT/.evotrace"
PREFIX="$ENV_ROOT/envs/evotrace"
mkdir -p "$ENV_ROOT/cache/conda-pkgs" "$ENV_ROOT/databases" "$ENV_ROOT/downloads" "$ENV_ROOT/logs" "$ENV_ROOT/envs"

echo "EvoTrace environment setup"
echo "Platform: $(uname -s) $(uname -m)"
echo "Project: $ROOT"

if command -v conda >/dev/null 2>&1; then
    CONDA="$(command -v conda)"
elif [[ -x "$HOME/miniforge3/bin/conda" ]]; then
    CONDA="$HOME/miniforge3/bin/conda"
elif [[ -x "$HOME/miniconda3/bin/conda" ]]; then
    CONDA="$HOME/miniconda3/bin/conda"
else
    echo "Conda/Mamba was not found. Creating a project-local Python environment." >&2
    if [[ -x "$ROOT/.venv/bin/python" ]]; then
        PY="$ROOT/.venv/bin/python"
    elif [[ -e "$ROOT/.venv" ]]; then
        echo "Refusing to replace a non-venv path: $ROOT/.venv" >&2
        exit 2
    else
        python3 -m venv "$ROOT/.venv"
        PY="$ROOT/.venv/bin/python"
    fi
    export PIP_CACHE_DIR="$ENV_ROOT/cache/pip"
    mkdir -p "$PIP_CACHE_DIR"
    "$PY" -m pip install --upgrade pip setuptools wheel
    "$PY" -m pip install -e '.[all]'
    export PATH="$(dirname "$PY"):$PATH"
    MISSING=()
    for tool in mafft hmmscan hmmpress; do command -v "$tool" >/dev/null || MISSING+=("$tool"); done
    command -v iqtree2 >/dev/null || command -v iqtree >/dev/null || command -v iqtree3 >/dev/null || MISSING+=("IQ-TREE")
    if ((${#MISSING[@]})); then
        echo "Python environment is ready, but required native tools are missing: ${MISSING[*]}. Install trusted native builds, then rerun this script." >&2
        exit 3
    fi
    EVOTRACE_PYTHON="$PY" bash "$ROOT/scripts/verify_environment.sh"
    PIP_CACHE_DIR="$ENV_ROOT/cache/pip" "$PY" -m pip cache purge >/dev/null 2>&1 || true
    exit 0
fi

export CONDA_PKGS_DIRS="$ENV_ROOT/cache/conda-pkgs"
export CONDA_ENVS_PATH="$ENV_ROOT/envs"
export CONDA_NOTICES_CACHE_DIR="$ENV_ROOT/cache/notices"
export PIP_CACHE_DIR="$ENV_ROOT/cache/pip"
mkdir -p "$CONDA_ENVS_PATH" "$CONDA_NOTICES_CACHE_DIR" "$PIP_CACHE_DIR"
export CONDA_CHANNEL_PRIORITY=strict
export CONDA_REGISTER_ENVS=false
export CONDA_NOTICES_ENABLED=false
if [[ -f "$PREFIX/conda-meta/history" ]]; then
    "$CONDA" env update --prefix "$PREFIX" --file environment.yml
elif [[ -e "$PREFIX" ]]; then
    echo "Refusing to overwrite non-Conda path: $PREFIX" >&2
    exit 2
else
    "$CONDA" env create --prefix "$PREFIX" --file environment.yml
fi

PY="$PREFIX/bin/python"
"$PY" -m pip install -e '.[all]'
EVOTRACE_PYTHON="$PY" bash "$ROOT/scripts/verify_environment.sh"
# Conda is constrained to this project's cache above; drop downloaded package
# archives only after installation, smoke tests, pytest, and strict doctor pass.
"$CONDA" clean --all --yes
