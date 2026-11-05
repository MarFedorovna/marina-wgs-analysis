#!/usr/bin/env bash
# Install report-only dependencies in a private server virtual environment.
set -Eeuo pipefail
umask 077
code_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
data_root=${MARINA_DATA_ROOT:-/srv/marina}
[[ $data_root == /* && $data_root != / ]] || { echo 'Use an absolute project data root.' >&2; exit 2; }
target=${data_root}/tools/report-venv
mkdir -p -- "${data_root}/tools"
exec 9>"${data_root}/tools/.report-env.lock"
flock -n 9 || { echo 'Report environment installation is already running.' >&2; exit 3; }
if [[ ! -x $target/bin/python ]]; then
  if ! python3 -c 'import ensurepip' >/dev/null 2>&1; then
    sudo -n apt-get update
    sudo -n apt-get install -y python3-venv
  fi
  python3 -m venv "$target"
fi
"$target/bin/python" -m pip install --disable-pip-version-check --index-url https://pypi.org/simple \
  -r "${code_root}/requirements-report.txt"
"$target/bin/python" -m pip check
"$target/bin/python" -c 'import importlib.metadata as m; import docx, matplotlib; print({n:m.version(n) for n in ("python-docx","lxml","typing_extensions","matplotlib")})'
"$target/bin/python" -m pip freeze > "${data_root}/tools/report-environment.lock.txt"
sha256sum "${code_root}/requirements-report.txt" > "${data_root}/tools/report-requirements.sha256"
printf 'Report environment READY: %s/bin/python\n' "$target"
