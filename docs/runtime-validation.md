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

No FASTQ alignment, variant calling, CUDA compute benchmark, or clinical validation
was performed by these setup checks. They establish installed tools and container GPU
exposure. Actual analysis still needs a pinned reference, private input manifests,
explicit interval/ploidy policy and a pilot. Setup logs and complete installed-package
versions remain under the private `/srv/marina/logs/`, not in git.

The root project directory, private config directory and output directory are owned
by the normal VM user with mode 0700. No new SSH users, credentials, public service
ports, Docker-group memberships or firewall changes were created by setup.
