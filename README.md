# WGS analysis automation

Private repository for reproducible server setup, data transfer utilities, and a pinned NVIDIA Parabricks workflow. This repository stores code and generic configuration examples. Actual FASTQ, BAM, VCF, database exports, real sample manifests, and connection details live outside Git history.

## Selected environment

The selected Brev environment uses Lambda `gpu_1x_a100_sxm4`: one NVIDIA A100 with 40 GB VRAM, 30 vCPU, 200 GiB quoted system RAM, and a fixed 512 GiB SSD. The console quoted **$2.39/hour including storage** on 2026-11-05; verify the current quote before any new deployment. Docker, NVIDIA Container Toolkit and the Parabricks image are installed, and hardware/container exposure and tool-version checks passed. See [runtime validation](docs/runtime-validation.md) for observed hardware, installed versions and the immutable image digest. No patient FASTQ alignment or variant calling has been run by these setup checks, and patient-analysis time has not been measured. Synthetic end-to-end validation is recorded separately in the runtime validation report.

This Brev offering is marked **pre-release** and has **no stop/start**. Charges continue until termination. Brev warns that exhausted organization credits cause the instance and its data to be deleted. Keep original inputs locally, download verified completed outputs during the run, and finish with a budget reserve. A planned budget is not a provider-enforced spending cap.

## Workflow

1. Prepare the raw Linux VM and verify GPU, driver, Docker GPU access, RAM, CPU and free disk space.
2. Build a local private transfer manifest. Review its file list before uploading.
3. Transfer the selected data over SSH/SFTP and verify SHA-256 on the server.
4. Supply a fixed reference FASTA, indexes, known-sites resources and sample/read-group metadata in local configuration. Do not infer the reference or sample metadata from filenames.
5. Process each sample through the same pinned pipeline: alignment, coordinate sorting, duplicate marking, BQSR when configured, and HaplotypeCaller gVCF generation.
6. Subsequent stage (not implemented by this runner): cohort joint genotyping, explicit filtering, annotation and interpretation.
7. Download outputs and logs, verify their checksums locally, and terminate the paid VM after confirming the retained copies.

One GPU is shared through a controlled queue. Four input samples do not imply that four memory-heavy GPU jobs run simultaneously. The scheduler must account for CPU threads, system RAM, GPU memory and scratch disk together.

## Software and reproducibility

The pinned Parabricks container is `nvcr.io/nvidia/clara/clara-parabricks:4.8.0-1`; the pulled immutable digest is recorded in [runtime validation](docs/runtime-validation.md). Record the resolved image digest, actual driver/tool versions, reference/resource checksums, command lines and QC metrics with each run. Docker packages the analysis tools; Kubernetes is not required for this single-server workload. Follow the [server instructions](docs/server-setup.md) and [verified transfer instructions](transfer/README.md) alongside their respective scripts.

Parabricks 4.8 documents GATK 4.3.0.0 and BWA-MEM 0.7.17 as compatibility baselines. Do not equate this with running an arbitrary newer GATK version on a GPU. Small output differences are documented, so compare normalized variants and genotypes rather than expecting identical compressed files. Technical rerunning cannot replace missing sequencing coverage.

## Server commands

To reproduce preparation on a fresh Ubuntu GPU VM, copy this reviewed checkout and upload the inputs with the [transfer utilities](transfer/README.md). Keep a private eight-column sample TSV based on `config/samples.example.tsv`, together with the expected-source upload manifest. From the checkout, inspect the plan first:

```bash
bash scripts/setup_project.sh \
  --samples /srv/marina/private/samples.tsv \
  --transfer-manifest /srv/marina/.transfer/PLAN_ID.json
```

Replace `PLAN_ID` with the actual private manifest basename. Repeat with `--execute` to install/verify Docker and GPU tools, install the official Ubuntu BWA package, prepare the locked Broad reference, verify the selected input pair, and create and validate a private technical-pilot plan. This prepares **one first sample** for full-pair alignment followed by diploid chr1 gVCF calling, with BQSR disabled. It does not create the cloud VM or transfer workstation data. Real analysis starts separately:

```bash
bash scripts/marina.sh status
bash scripts/marina.sh pilot --wait
```

For preparation followed immediately by this pilot, add both `--execute --pilot` to the setup command. The [clean notebook template](notebooks/START_HERE.ipynb) provides separate status and manual pilot cells; keep its executed working copy outside the checkout. A four-sample final WGS/cohort analysis requires a separate reviewed plan. The setup wrapper's syntax and default plan have been checked; its complete sequence has not been rerun on another fresh VM. Individual setup/reference stages and synthetic pipeline checks are recorded in the validation documents.

Copy only the reviewed code to the Linux VM and follow [server setup](docs/server-setup.md). Run the bootstrap once as an administrator, then the preflight before processing:

```bash
sudo bash scripts/bootstrap-server.sh /srv/marina
sudo bash scripts/preflight-server.sh /srv/marina
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
- [Pinned Broad GRCh38 reference preparation](docs/reference-setup.md)
- [Observed runtime validation](docs/runtime-validation.md)
- [Verified data transfer](transfer/README.md)
- [Operations and budget](docs/operations.md)
- [Official Parabricks installation requirements](https://docs.nvidia.com/clara/parabricks/get-started/installation-requirements)
- [Parabricks compatibility baselines](https://docs.nvidia.com/clara/parabricks/about-parabricks/software-overview/output-accuracy-and-compatible-cpu-software-versions)
- [Parabricks FQ2BAM](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/fq2bam)
- [Parabricks HaplotypeCaller](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/haplotypecaller)
