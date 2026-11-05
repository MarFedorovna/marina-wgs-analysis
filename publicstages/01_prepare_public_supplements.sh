#!/usr/bin/env bash
# Additional public-source data only. Default is read-only status.
set -Eeuo pipefail
umask 077
code_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
exec python3 -u "${code_root}/scripts/install_public_supplements.py" "$@"
