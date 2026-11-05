# WGS analysis automation

Private repository for reproducible server setup, data transfer utilities, and a pinned NVIDIA Parabricks workflow. This repository stores code and generic configuration examples. Actual FASTQ, BAM, VCF, database exports, real sample manifests, and connection details live outside Git history.

## Selected environment

The deployment prepared in Brev uses Lambda `gpu_1x_a100_sxm4`: one NVIDIA A100 with 40 GB VRAM, 30 vCPU, 200 GiB system RAM, and a fixed 512 GiB SSD. The console quoted **$2.39/hour including storage** on 2026-11-05; verify the current quote before deployment. Provisioning and processing have not been measured by this repository.

This Brev offering is marked **pre-release** and has **no stop/start**. Charges continue until termination. Brev warns that exhausted organization credits cause the instance and its data to be deleted. Keep original inputs locally, download verified completed outputs during the run, and finish with a budget reserve. A planned budget is not a provider-enforced spending cap.

## Workflow

1. Prepare the raw Linux VM and verify GPU, driver, Docker GPU access, RAM, CPU and free disk space.
2. Build a local private transfer manifest. Review its file list before uploading.
3. Transfer the selected data over SSH/SFTP and verify SHA-256 on the server.
4. Supply a fixed reference FASTA, indexes, known-sites resources and sample/read-group metadata in local configuration. Do not infer the reference or sample metadata from filenames.
5. Process each sample through the same pinned pipeline: alignment, coordinate sorting, duplicate marking, BQSR when configured, and HaplotypeCaller gVCF generation.
6. Complete cohort genotyping and explicit filtering with a compatible, pinned workflow, then review technical QC before annotation or interpretation.
7. Download outputs and logs, verify their checksums locally, and terminate the paid VM after confirming the retained copies.

One GPU is shared through a controlled queue. Four input samples do not imply that four memory-heavy GPU jobs run simultaneously. The scheduler must account for CPU threads, system RAM, GPU memory and scratch disk together.

## Software and reproducibility

The planned Parabricks container is `nvcr.io/nvidia/clara/clara-parabricks:4.8.0-1`. Record the resolved image digest, actual driver/tool versions, reference/resource checksums, command lines and QC metrics with each run. Docker packages the analysis tools; Kubernetes is not required for this single-server workload. Server and transfer instructions are added alongside their respective scripts.

Parabricks 4.8 documents GATK 4.3.0.0 and BWA-MEM 0.7.17 as compatibility baselines. Do not equate this with running an arbitrary newer GATK version on a GPU. Small output differences are documented, so compare normalized variants and genotypes rather than expecting identical compressed files. Technical rerunning cannot replace missing sequencing coverage.

## Server commands

Copy only the reviewed code to the Linux VM and follow [server setup](docs/server-setup.md). Run the bootstrap once as an administrator, then the preflight before processing:

```bash
sudo bash scripts/bootstrap-server.sh /srv/marina
bash scripts/preflight-server.sh /srv/marina
```

Create the private run configuration from the generic examples, using `/srv/marina/private/` for real sample details. Inspect generated commands before execution:

```bash
python3 scripts/run_parabricks.py --config /srv/marina/private/run.json --stage align --sample SAMPLE_01
```

The default is a dry run. Add `--execute` only after supplying and checking the actual reference/resources, data checksums, read groups, explicit call intervals/ploidy, free space and budget. The server guide documents the remaining stages and their ordering.

## Repository boundaries

- Publish only reviewed scripts, generic examples, dependency pins and documentation.
- Keep real manifests and connection settings under ignored `.local/`.
- Do not put raw reads, aligned reads, patient-level variants, phenotype tables, SQL dumps, SSH keys, tokens or passwords in commits or GitHub Actions artifacts.
- `.gitignore` applies to Git operations. Browser file uploads require a separate reviewed file allowlist.
- Do not put the analysis workload on GitHub-hosted runners; the paid GPU VM performs computation.

## Documentation

- [Server setup and analysis commands](docs/server-setup.md)
- [Operations and budget](docs/operations.md)
- [Official Parabricks installation requirements](https://docs.nvidia.com/clara/parabricks/get-started/installation-requirements)
- [Parabricks compatibility baselines](https://docs.nvidia.com/clara/parabricks/about-parabricks/software-overview/output-accuracy-and-compatible-cpu-software-versions)
- [Parabricks FQ2BAM](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/fq2bam)
- [Parabricks HaplotypeCaller](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/haplotypecaller)
