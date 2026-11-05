#!/usr/bin/env bash
# Continue after a completed filtered cohort VCF and verified offline resources.
set -Eeuo pipefail
code_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
python3 "${code_root}/scripts/workflow.py" annotate "$@"
exec python3 "${code_root}/scripts/workflow.py" report "$@"
