# Gene coverage QC

`scripts/gene_coverage.py --config /absolute/private/config.json --execute`
is a read-only, serial per-sample workflow. The configuration contains
`reference`, `reference_sha256`, `sources_lock`, `mosdepth`, `gtf`, `clingen_csv`,
`output_root`, `status_file`, and `samples` records with `sample_id`, `bam`,
`receipt`. Keep the real configuration outside Git. The reference identity is
carried from the completed alignment manifest; BAM sequence lengths are checked
against its FAI. The BAM is not rehashed or modified.

The source lock pins GENCODE49, an observed official mosdepth0.3.11 binary, and a
dated ClinGen gene-validity CSV. The GTF is downloaded directly onto the server
and checked against the publisher's MD5 and size. The workflow first runs an
actual synthetic BAM test proving mate-overlap correction, BQ0 inclusion and
exclusion of low-MQ/duplicate/secondary/supplementary/QC-failed reads. A four-CPU
affinity mask and16GiB address-space limit bound the process and its children.

The model is the union of all protein-coding CDS intervals within each full
versioned gene identifier and chromosome. GTF1-based coordinates become BED
0-based half-open coordinates. Adjacent intervals merge. It is an all-transcript
CDS model, not a canonical transcript or a clinically validated exon panel.
Only chr1–22 are measured. X/Y and mitochondrial/noncoding/missing-model genes
are not assigned zero coverage; they remain not assessed/NA.

One mosdepth BAM pass uses MQ20 and excluded flags0xF04. It retains default mate
overlap and CIGAR handling, with neither fast-mode nor fragment-mode. **There is
no base-quality filter.** These coverage metrics do not certify a confident
genotype. Per-region means are combined by interval length; threshold counts
are exact integers. A separate global union in the same BED prevents genes
with overlapping CDS from inflating the global denominator.

Each completed sample contains `gene-cds-coverage.tsv`,
`clingen-hcm-cds-coverage.tsv`, `coverage-summary.json`, mosdepth output and a
`workflow-receipt.json`. Gene metrics are CDS bases, mean depth, bases and
percentages at least1/10/20X, and zero-covered bases. The denominator includes
reference gaps within the CDS and is not a callable-genome fraction. HCM rows
use exact MONDO:0005045 and exact HGNC identifiers, retain all evidence
classifications and inheritance modes, and preserve unmatched rows as NA.
There is no automatic diagnosis or sex/ploidy assignment.

Run the small model tests with `python3 scripts/test_gene_coverage.py`.
Receipts are emitted only after the sample completes; failed partial outputs
remain available for diagnosis and are not advertised as completed results.

Method references: [GENCODE GTF format](https://www.gencodegenes.org/pages/data_format.html),
[mosdepth documentation](https://github.com/brentp/mosdepth),
[ClinGen gene validity](https://search.clinicalgenome.org/kb/gene-validity).
