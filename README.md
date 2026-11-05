# WGS Analysis

Reproducible workflows and documentation for research WGS analysis. This repository contains code and methods only; it must not contain patient genomic files or identifiable clinical data.

## Project status

Repository scaffold created. The analysis workflow is being defined. Existing VCFs are provisional until read quality, mapped coverage, sample mapping, and a consistent variant-calling workflow have been reviewed.

## Planned workflow

1. Maintain a pseudonymized sample manifest and checksums for approved input files.
2. Run FASTQ quality control and summarize it with FastQC/MultiQC.
3. Pin the GRCh38 reference, tool versions, and workflow configuration.
4. Align reads and measure mapped depth and coverage breadth.
5. Call variants consistently across samples; retain raw and filtered outputs with documented criteria.
6. Annotate variants and prioritize them using phenotype, inheritance, population frequency, and gene-specific guidance.
7. Review candidate LoF, stop-gain, stop-loss, frameshift, splice, and other variants with read-level evidence.
8. Record limitations, provenance, and review decisions in the final report.

## Data handling

Never commit FASTQ, BAM/CRAM, VCF/BCF, reference genomes, sample-level phenotype files, pedigrees, credentials, or cloud keys. The `.gitignore` is a safety net, not a substitute for checking staged files before every commit. Store sensitive inputs and outputs only in an approved restricted location.

## Reproducibility

Each analysis release should record the reference build and checksum, software/container versions, parameters, input manifest checksum, QC results, and output checksums. Do not publish clinical interpretations as automated diagnoses; candidate classifications require appropriate expert review.
