#!/usr/bin/env bash
set -Eeuo pipefail
umask 077
task_root=${1:-/srv/marina}
image=${PARABRICKS_IMAGE:-nvcr.io/nvidia/clara/clara-parabricks:4.8.0-1}
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { echo 'Linux x86_64 host required.' >&2; exit 1; }
mkdir -p "${task_root}/logs"
exec > >(tee "${task_root}/logs/preflight-$(date -u +%Y%m%dT%H%M%SZ).log") 2>&1
for executable in docker nvidia-smi python3 samtools bcftools sha256sum flock; do
  command -v "${executable}" >/dev/null || { echo "Missing: ${executable}" >&2; exit 1; }
done
(( $(nproc) >= 24 )) || { echo 'At least 24 available CPU threads required for this configuration.' >&2; exit 1; }
awk '/MemTotal/ {exit !($2 >= 97656250)}' /proc/meminfo || { echo 'At least 100 GB total system RAM required.' >&2; exit 1; }
free_bytes=$(df -PB1 "${task_root}" | awk 'NR==2 {print $4}')
(( free_bytes >= 161061273600 )) || { echo 'Less than 150 GiB free; review scratch/output capacity.' >&2; exit 1; }
gpu_memory=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -n 1 | tr -d ' ')
(( gpu_memory >= 36000 )) || { echo 'Selected profile requires about 40 GB VRAM; GPU 0 is smaller.' >&2; exit 1; }
if nvidia-smi -L | grep -q 'MIG '; then
  echo 'MIG device detected. Parabricks does not support MIG.' >&2; exit 1
fi
nvidia-smi
docker version
docker image inspect "${image}" --format '{{json .RepoDigests}}'
docker run --rm --network none --gpus all "${image}" nvidia-smi
docker run --rm --network none --gpus all "${image}" pbrun fq2bam --version
samtools --version | head -n 2
bcftools --version | head -n 1
df -h "${task_root}"
echo 'Hardware/runtime preflight passed; actual CUDA workload compatibility still needs the first pilot.'
