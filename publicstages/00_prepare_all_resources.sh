#!/usr/bin/env bash
# Run on the server. Reuse active downloads; no patient inputs are read.
set -Eeuo pipefail
umask 077
code_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { echo 'Run on the Linux x86_64 server.' >&2; exit 1; }
exec python3 -u "${code_root}/scripts/prepare_all_resources.py" "$@"
