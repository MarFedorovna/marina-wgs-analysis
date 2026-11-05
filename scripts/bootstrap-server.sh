#!/usr/bin/env bash
# Run on the cloud Ubuntu host, never on the Windows workstation.
set -Eeuo pipefail
umask 077
[[ ${EUID} -eq 0 ]] || { echo 'Run with sudo.' >&2; exit 1; }
[[ $(uname -m) == x86_64 ]] || { echo 'This bootstrap targets Ubuntu x86_64.' >&2; exit 1; }
source /etc/os-release
[[ ${ID} == ubuntu && ${VERSION_ID} =~ ^(22\.04|24\.04|26\.04)$ ]] || {
  echo 'Supported here: Ubuntu 22.04, 24.04, or 26.04 LTS.' >&2; exit 1;
}
task_root=${1:-/srv/marina}
[[ ${task_root} == /* && ${task_root} != / ]] || { echo 'Use an absolute project data directory.' >&2; exit 1; }
install -d -m 0700 "${task_root}" "${task_root}/logs" "${task_root}/private" \
  "${task_root}/data" "${task_root}/references" "${task_root}/outputs" "${task_root}/scratch"
exec > >(tee -a "${task_root}/logs/bootstrap-$(date -u +%Y%m%dT%H%M%SZ).log") 2>&1
trap 'echo "Bootstrap failed at line ${LINENO}; no analysis has been started." >&2' ERR
command -v nvidia-smi >/dev/null || { echo 'GPU driver missing. Select a GPU-ready provider image; this script does not replace drivers.' >&2; exit 1; }
nvidia-smi
apt-get update
apt-get install -y --no-install-recommends ca-certificates curl gnupg python3 git rsync tmux jq pigz samtools bcftools

# Preserve provider Docker installations and running containers.
if ! command -v docker >/dev/null; then
  for package in docker.io podman-docker containerd runc; do
    if dpkg-query -W -f='${Status}' "${package}" 2>/dev/null | grep -q 'install ok installed'; then
      echo "Conflicting package ${package}; review it manually rather than removing provider software." >&2
      exit 1
    fi
  done
  install -m 0755 -d /etc/apt/keyrings
  curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
  chmod a+r /etc/apt/keyrings/docker.asc
  cat > /etc/apt/sources.list.d/docker.sources <<EOF
Types: deb
URIs: https://download.docker.com/linux/ubuntu
Suites: ${UBUNTU_CODENAME:-${VERSION_CODENAME}}
Components: stable
Architectures: $(dpkg --print-architecture)
Signed-By: /etc/apt/keyrings/docker.asc
EOF
  apt-get update
  if [[ -n ${DOCKER_VERSION:-} ]]; then
    apt-get install -y docker-ce="${DOCKER_VERSION}" docker-ce-cli="${DOCKER_VERSION}" containerd.io docker-buildx-plugin docker-compose-plugin
  else
    apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
  fi
fi
systemctl start docker

cuda_image=nvcr.io/nvidia/cuda:13.3.1-base-ubuntu22.04
if ! docker run --rm --gpus all "${cuda_image}" nvidia-smi; then
  if [[ -n $(docker ps -q) ]]; then
    echo 'GPU runtime failed while containers are active. Refusing to restart their Docker daemon.' >&2
    exit 1
  fi
  if ! command -v nvidia-ctk >/dev/null; then
    curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | \
      gpg --batch --yes --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg
    curl -fsSL https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
      sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' \
      > /etc/apt/sources.list.d/nvidia-container-toolkit.list
    apt-get update
    toolkit_version=${NVIDIA_CONTAINER_TOOLKIT_VERSION:-1.20.1-1}
    apt-get install -y nvidia-container-toolkit="${toolkit_version}" \
      nvidia-container-toolkit-base="${toolkit_version}" \
      libnvidia-container-tools="${toolkit_version}" libnvidia-container1="${toolkit_version}"
  fi
  if [[ -f /etc/docker/daemon.json ]]; then
    cp -p /etc/docker/daemon.json "${task_root}/logs/docker-daemon-before-$(date -u +%Y%m%dT%H%M%SZ).json"
  fi
  nvidia-ctk runtime configure --runtime=docker
  systemctl restart docker
  docker run --rm --gpus all "${cuda_image}" nvidia-smi
fi

parabricks_image=nvcr.io/nvidia/clara/clara-parabricks:4.8.0-1
docker pull "${parabricks_image}"
docker image inspect "${parabricks_image}" > "${task_root}/logs/parabricks-image-inspect.json"
docker run --rm --network none --gpus all "${parabricks_image}" pbrun fq2bam --version
dpkg-query -W -f='${Package}\t${Version}\n' > "${task_root}/logs/installed-packages.tsv"
if [[ -n ${SUDO_USER:-} && ${SUDO_USER} != root ]]; then
  chown "${SUDO_UID}:${SUDO_GID}" "${task_root}" "${task_root}/logs" "${task_root}/private" \
    "${task_root}/data" "${task_root}/references" "${task_root}/outputs" "${task_root}/scratch"
  find "${task_root}/logs" -maxdepth 1 -type f -user root -exec chown "${SUDO_UID}:${SUDO_GID}" {} +
fi
echo 'Docker + GPU runtime are ready. No FASTQ processing or data upload was performed.'
echo 'Run scripts/preflight-server.sh next. Copy RepoDigests[0] into the private run.json before processing.'
