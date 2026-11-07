#!/usr/bin/env bash
set -Eeuo pipefail

# Install WhatsHap in an isolated venv; do not change the system Python.
python_bin="${PYTHON_BIN:-python3}"
env_dir="${WHATSHAP_ENV_DIR:-${HOME}/.venvs/whatshap-2.8}"

command -v "$python_bin" >/dev/null || {
  echo "Python 3 is required." >&2
  exit 2
}

if [[ ! -x "$env_dir/bin/python" ]]; then
  mkdir -p "$(dirname "$env_dir")"
  "$python_bin" -m venv "$env_dir"
fi

"$env_dir/bin/python" -m pip install --upgrade pip
"$env_dir/bin/python" -m pip install 'whatshap==2.8'
"$env_dir/bin/python" -m pip freeze > "$env_dir/requirements.lock.txt"
"$env_dir/bin/whatshap" --version
echo "WhatsHap environment ready: $env_dir"
