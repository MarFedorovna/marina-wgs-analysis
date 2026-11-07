# Read-backed phasing and annotation coverage

## Phasing needs a verified VCF and matching BAMs

WhatsHap uses aligned reads from BAM/CRAM to phase nearby variants. It writes a new phased VCF; it does not change the source VCF. A phase set does not identify maternal or paternal origin.

Before running, verify all input files against a private transfer receipt. Keep sample IDs, paths, hashes, BAM/VCF files, and patient-derived data out of this public repository.

## Server setup

On Ubuntu/Jupyter Terminal, from the repository root:

```bash
bash scripts/install_whatshap.sh
export PATH="$HOME/.venvs/whatshap-2.8/bin:$PATH"
```

The installer uses a virtual environment and pins WhatsHap 2.8 without replacing system Python. The workflow also needs bcftools, samtools, tabix, a matching reference FASTA plus .fai, an indexed VCF, and one coordinate-sorted BAM per sample. BAM read groups must have an SM value matching the VCF sample name.

## Stop if inputs do not match

Compare the complete VCF SHA-256 to the verified source receipt before analysis. Do not use a VCF whose checksum differs, even if its file size looks correct. Do not run phasing until the correct VCF and all matching BAMs are present.

Create a private TSV manifest with columns `sample_id` and `bam_path`. Never commit this manifest. Inspect each BAM with `samtools quickcheck`, confirm its index, and verify the `@RG SM` tag.

## Annotation fields

A field in a VCF header does not mean every variant has a score. Report missing values as missing; do not replace them with zero or infer that a tool ran from a column name. SpliceAI has a bounded model scope. FATHMM, CADD, REVEL, dbNSFP, AlphaMissense, ACMG, OMIM, LOVD, GeneReviews, and VarSome must each be checked against the exact input artifact and the permitted source/version.

Read-backed phasing is separate from annotation and does not assign ACMG classifications.

Official WhatsHap guide: https://whatshap.readthedocs.io/en/stable/guide.html
