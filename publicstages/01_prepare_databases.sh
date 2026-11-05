#!/usr/bin/env bash
# Run ON THE SERVER. No patient inputs or workflow configuration are required.
set -Eeuo pipefail
umask 077
code_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
data_root=${MARINA_DATA_ROOT:-/srv/marina}
workers=${MARINA_DOWNLOAD_WORKERS:-16}
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { echo 'Run this on the Linux x86_64 server.' >&2; exit 1; }
[[ $data_root == /* && $data_root != / ]] || { echo 'Use an absolute server data root.' >&2; exit 2; }
python3 -c 'import sys; assert sys.version_info >= (3,9), "Python 3.9+ required"'
exec python3 -u "${code_root}/scripts/prepare_annotation_resources.py" \
  --config "${code_root}/config/annotation-resources.example.json" \
  --root "${data_root}/resources/annotation" \
  --direct-cache-dir "${data_root}/resources/.downloads" \
  --download-workers "$workers" --reserve-gib 150 --execute --pull-image
