#!/usr/bin/env bash
# Reproduce the prepared environment/reference/pilot configuration on an Ubuntu GPU VM.
set -Eeuo pipefail
umask 077
code_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
data_root=/srv/marina
samples=${data_root}/private/samples.tsv
transfer_manifest=
execute=false
start_pilot=false

usage() {
  cat <<'HELP'
Usage: bash scripts/setup_project.sh [--samples /absolute/samples.tsv]
       [--transfer-manifest /absolute/server-transfer-manifest.json]
       [--execute] [--pilot]

Default: print the commands without installing, downloading, writing config or running analysis.
--execute: prepare Docker/GPU tools, the locked Broad GRCh38 reference and a private pilot plan.
--pilot: after successful preparation, start the first sample technical pilot (requires --execute).

Run as the normal sudo-capable VM user, from a reviewed checkout. This does not create a
cloud VM or upload workstation files. Transfer the input files with transfer/ utilities
first. The transfer manifest supplies the expected original read hashes; it stays private.
The sample TSV uses config/samples.example.tsv's columns. It must describe your own inputs.

The pilot aligns one full FASTQ pair and calls diploid chr1 only. It is not a final
whole-genome/cohort/clinical analysis. Billing continues after this command exits.
HELP
}

while (($#)); do
  case "$1" in
    --samples|--transfer-manifest)
      (($# >= 2)) || { echo "Missing value for $1" >&2; exit 2; }
      [[ $2 == /* ]] || { echo "Use an absolute server path for $1." >&2; exit 2; }
      if [[ $1 == --samples ]]; then samples=$2; else transfer_manifest=$2; fi
      shift 2 ;;
    --execute) execute=true; shift ;;
    --pilot) start_pilot=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done
if $start_pilot && ! $execute; then
  echo '--pilot requires --execute; the default is a plan only.' >&2
  exit 2
fi

reference_receipt=${data_root}/references/broad-hg38-v0-altaware-g1575676516681666-prebuilt/reference-receipt.json
configure=(python3 "${code_root}/scripts/configure_pilot.py" --root "$data_root"
  --samples "$samples" --reference-receipt "$reference_receipt")
if [[ -n $transfer_manifest ]]; then configure+=(--transfer-manifest "$transfer_manifest"); fi
configure+=(--execute)

run() {
  printf '+ '
  printf '%q ' "$@"
  printf '\n'
  if $execute; then "$@"; fi
}

if $execute; then
  [[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { echo 'Linux x86_64 GPU VM required.' >&2; exit 1; }
  ((EUID != 0)) || { echo 'Run as the normal VM user; this script invokes sudo only for installation.' >&2; exit 1; }
  [[ -f $samples ]] || { echo 'Supply the private sample TSV before setup.' >&2; exit 1; }
  if [[ -n $transfer_manifest ]]; then
    [[ -f $transfer_manifest ]] || { echo 'The specified server transfer manifest is missing.' >&2; exit 1; }
  else
    [[ -f ${data_root}/private/pilot-inputs-verified.json ]] || {
      echo 'Cold setup requires --transfer-manifest with the expected original input hashes.' >&2; exit 1;
    }
  fi
  [[ -f ${code_root}/scripts/configure_pilot.py ]] || { echo 'Incomplete reviewed code checkout.' >&2; exit 1; }
  sudo -n true || { echo 'Non-interactive sudo is required on the selected VM profile.' >&2; exit 1; }
fi

run sudo -n bash "${code_root}/scripts/bootstrap-server.sh" "$data_root"
run sudo -n apt-get install -y --no-install-recommends bwa
run sudo -n bash "${code_root}/scripts/preflight-server.sh" "$data_root"
run python3 "${code_root}/scripts/prepare_reference.py" \
  --lock "${code_root}/config/broad-hg38-v0.sources.lock.json" \
  --root "${data_root}/references" --execute
run "${configure[@]}"
run bash "${code_root}/scripts/marina.sh" plan
if $start_pilot; then
  run bash "${code_root}/scripts/marina.sh" pilot --wait
else
  echo 'Preparation ends before real analysis. Start the prepared pilot with:'
  printf 'bash %q pilot --wait\n' "${code_root}/scripts/marina.sh"
fi
