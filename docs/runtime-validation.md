# Runtime validation on the selected cloud GPU

Checked on 2026-11-05 after explicit authorization to create/install the VM.

| Component | Observed value |
|---|---|
| Host context | KVM VM, Ubuntu 22.04.5 LTS; not inside Docker |
| GPU | NVIDIA A100-SXM4-40GB, 40960 MiB; MIG disabled |
| CPU | 30 available threads |
| Physical RAM | 232205721600 bytes (about 216 GiB) |
| Root disk | 497 GiB filesystem, about 462 GiB available after container pull |
| Driver | 570.148.08, preserved without replacement |
| Docker Engine | Existing provider installation 28.3.1, preserved |
| NVIDIA Container Toolkit | Existing provider installation 1.17.8, preserved |
| Host samtools / bcftools | Ubuntu packages 1.13 / 1.13 |
| Parabricks | 4.8.0-1 |

The pulled immutable container reference is:

```text
nvcr.io/nvidia/clara/clara-parabricks@sha256:44206e2b522e87370c1107c774f52c5859af9b8f3408aa7232de4bbe3c24953c
```

Hardware preflight passed. `nvidia-smi` inside the CUDA 13.3.1 and Parabricks images
sees the physical GPU. Parabricks `fq2bam --version` and `haplotypecaller --version`
completed; the latter was also checked with a non-root container UID:GID.

The actual tool outputs report these comparison baselines:

- fq2bam: BWA-MEM 0.7.17 and Picard 2.18.25.
- haplotypecaller: GATK 4.3.0.0.

Record the actual command output when comparing CPU/GPU results. A generic online
software-version table may differ from the tool-specific binary report.

## Minimal CUDA compute verification

On 2026-11-05, a separate synthetic compute probe completed with **PASS** inside the
same immutable Parabricks image, using the normal VM user's UID:GID (1000:1000).
It used no genomic inputs and allocated only four bytes of GPU output memory.

The probe initialized CUDA, called the CUDA runtime's `cudaMalloc` and `cudaMemset`,
JIT-compiled a minimal PTX kernel targeting the A100's `sm_80` architecture, launched
that kernel on the GPU, synchronized, and copied the result back with `cudaMemcpy`.
The returned value was **42**, matching the kernel's expected output. CUDA Driver API
and runtime version queries both returned **13030 (CUDA 13.3)**.

The process's `/proc/self/maps` confirmed these actual loaded libraries:

```text
/usr/local/cuda-13.3/compat/libcuda.so.610.57.04
/usr/local/cuda-13.3/compat/libnvidia-ptxjitcompiler.so.610.57.04
/usr/local/cuda-13.3/targets/x86_64-linux/lib/libcudart.so.13.3.29
```

Thus the image's CUDA forward-compatibility libraries were used over the unchanged
host kernel driver **570.148.08**. NVIDIA's [forward compatibility documentation and
support matrix](https://docs.nvidia.com/deploy/cuda-compatibility/forward-compatibility.html),
consulted on 2026-11-05, mark `cuda-compat-13-3` with the 570 driver branch as compatible
for supported Data Center GPUs. This is distinct from plain [minor-version
compatibility](https://docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html),
where CUDA 13.x normally requires at least the 580 driver branch.

The successful probe establishes real minimal GPU execution with this image and host;
it does not establish performance or correctness of a bioinformatics workload. No
FASTQ alignment, variant calling, analysis-workload benchmark, or clinical validation
was performed. A pinned reference, private input manifests, explicit interval/ploidy
policy and a representative pilot are still required for actual analysis. The probe
code/log remain in the ignored local directory; setup logs and complete installed
package versions remain under private `/srv/marina/logs/`, not in git. Neither the
host driver nor the stable analysis code was changed to obtain this compute result.

The root project directory, private config directory and output directory are owned
by the normal VM user with mode 0700. No new SSH users, credentials, public service
ports, Docker-group memberships or firewall changes were created by setup.
