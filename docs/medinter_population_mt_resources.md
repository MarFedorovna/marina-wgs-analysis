# MedInter P01, P02 and P08 public reference supplements

These additions follow the explicit approval of MedInter P01–P11 and the confirmed noncommercial research purpose. They create isolated public-reference sidecars on the VM. They do not change the sealed calling/VEP runtime, the current controllers, stage 6 bindings, or any patient outputs. The reference datasets remain on the server; the public code package contains only scripts, configurations, tests and this documentation.

## P01: population counts, denominators and coverage context

The fixed public panel is the existing ClinGen exact-HCM/Released CSpec panel: 39 genes with GENCODE49 GRCh38 CDS intervals, selected without patient variants. The extractor queries official indexed **gnomAD v4.1.1** exome and genome VCFs using HTTPS byte ranges. It preserves AC, AN, AF, nhomalt, FILTER, ancestry/subset fields and quality fields. It omits only redundant `vep`/`CSQ` consequence payloads from the sidecar; full source headers and cached source byte ranges retain provenance. Missing values remain missing, including SQL NULL for missing AN. A missing variant record does not automatically mean AC=0.

The combined source range budget is 1,900,000,000 bytes, with a 150 GiB analysis reserve and the unchanged absolute deadline `2026-11-06T03:46:48.046025Z`. The source cache and SQLite interval transactions support resume. The final extraction completed at 19:35:55 UTC with **220,733 public variant/source rows across all 2,262 dataset/interval combinations**, using **1,589,726,509 HTTP body bytes**. The first large uncompressed sidecars remain preserved under `medinter-population-hcm39-v1` and `v2`; the compact authoritative successor is **`medinter-population-hcm39-v3`**. Its INFO column contains losslessly zlib-compressed UTF-8 BLOBs, decoded by `medinter_population_mt_common.decode_population_info`; AC/AN/AF/nhomalt/FILTER remain separately queryable. This preserves every population, ancestry and quality field while avoiding a large duplicate text payload. Prior committed intervals and verified source tiles were reused; no current clinical pipeline controller was stopped.

The official exome all-sites AN v4.1 file is **1,149,738,480 bytes**, independently within the modest source cap. Its published MD5 is `ee71ce1ccfa5c1d9dd86e1ee1b1d11e2`. No tabix index is published at the tested path, so it is retained on the server and filtered locally to the same public CDS panel. The resulting sidecar contains **243,349 positions**, including **44 explicitly zero-AN positions**, across a 244,998-base CDS union; positions absent from the subset are not assigned AN=0 or AC=0. Genome all-sites AN is 11.19 GiB and the exome coverage file is 3.77 GiB; those large unindexed sources are not downloaded.

Coverage context instead uses a coordinated public gene-only GraphQL request sequence, at most one request per seven seconds, after the separate pext job finished: **38 nuclear genes, 24,650,850 bytes**, complete at 19:26:06 UTC. The [official source mapping](https://raw.githubusercontent.com/broadinstitute/gnomad-browser/main/graphql-api/src/queries/coverage-queries.ts) maps the `gnomad_r4` coverage query to **v4 exome coverage and v3 genome coverage**. This is a different source version from the v4.1.1 allele-count VCFs. The gene API aggregates across padded merged exons; displayed means/medians and rounded fractions are binned context, not a claim of exact base-level depth. MT-TI is explicitly not applicable to this nuclear coverage dataset. The API does not return an exact coverage release timestamp.

Current official [gnomAD terms](https://raw.githubusercontent.com/broadinstitute/gnomad-browser/main/browser/about/policies/terms.md) identify primary data as CC0 and request attribution; separate annotations can have separate licenses. This corrects the historical ODbL note in the proposal. The source terms and official download manifest are retained verbatim with retrieval hashes. Attribution: gnomAD consortium, exact releases named above; [official downloads](https://gnomad.broadinstitute.org/downloads).

## P02: mtDNA-specific public reference populations

Actual source downloads:

| Source | Bytes | Validation |
| --- | ---: | --- |
| Official HelixMTdb 20200327 TSV | 2,050,188 | Exact byte count, transport identity, SHA256 and tabular integrity |
| Official gnomAD v3.1 chrM VCF | 5,150,029 | Exact byte count, publisher CRC32C, SHA256 and full gzip integrity |
| Its official tabix index | 108 | Exact byte count, publisher CRC32C and tabix magic |
| Official NC_012920.1 rCRS FASTA | 16,864 | Accession/length checked; sequence exactly matches local public GRCh38 chrM |

The ready SQLite sidecar contains **18,164 gnomAD MT records** and **14,104 biallelic Helix records**. All accepted REF alleles match rCRS; no REF mismatches occurred. **220 Helix multiallelic records** are preserved separately with `MULTIALLELIC_COUNTS_NOT_ASSIGNED_TO_SINGLE_ALLELE`, instead of assigning combined counts to each alternate allele. Helix does not supply a per-site AN column, so AN remains NULL. Its total published cohort size is not substituted as a per-site denominator. Homoplasmic/heteroplasmic counts and haplogroup context remain distinct. Neither Helix nor GenBank is asserted to be a proven healthy control cohort.

The [official Helix page](https://www.helix.com/mitochondrial-variant-database) links the versioned source and requests citation/linking; the source page is retained, and this package does not assert an unrestricted redistribution license for Helix content. gnomAD's MT release remains v3.1; it is not silently relabelled as a new v4 MT callset.

MITOMAP's [official Resources endpoint](https://www.mitomap.org/MITOMAP/Resources) returned **HTTP 403** from the VM. Its content and license are not downloaded or assumed. No alternate authentication, browser profile access or access-control bypass was attempted. Consequently P02 has a real MITOMAP evidence gap, despite the two working frequency sources.

## P08: ClinVar submission and literature dossier preparation

Verified official NCBI source files total about 1.128 GB: submission summary, variant summary, variant citations, gene summary, organization summary and README. The main tables' publisher MD5 files are preserved and verified. These are compact tabular references; no full ClinVar XML or paywalled article collection is downloaded. Their late-October/1-November 2026 snapshot is tracked separately from the existing ClinVar 25-October-2026 core snapshot.

The local builder selects records using only the fixed public39-gene panel and GRCh38. It retains public submitted classifications, submitters, dates, conditions and raw source rows, and prepares a deduplicated citation lookup/manual literature queue. Gene overlap or a submitted ClinVar classification is not promoted to HCM variant pathogenicity. Public references may include different variant classes; this is reference preparation and performs no CNV/SV calling or analysis.

The actual selected submission-summary data contain **248,946 versioned SCV rows, with zero missing versions**, covering **109,528 public Variation IDs**. The variant-summary aggregate SCV lists omit versions; they are not substituted for the versioned submission records. The builder prepares **207,404 deduplicated variant/citation links**. Original articles have not been read by this automated source-preparation job, their full texts are not acquired, functional assays are not validated, and overlapping submissions/papers are not counted as independent cases. ClinVar classification is not assigned as independent PP5 evidence. The public dossier template makes these gaps explicit; actual case worksheets belong in the private project area.

NCBI documents the distinction between [ClinVar aggregate and submitted records](https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/) and [review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/). Publication and submitter rights/provenance remain applicable; a citation or public abstract is not permission to redistribute a full article.

## Status and controls

```sh
python3 /srv/marina/private/medinter-population-mt-v1/scripts/medinter_population_mt_status.py
python3 -m unittest discover -s /srv/marina/private/medinter-population-mt-v1/tests -p test_medinter_population_mt.py -v
```

Authoritative new resources are under `/srv/marina/resources/medinter-population-hcm39-v3`, `medinter-population-gene-coverage-v1`, `medinter-population-mt-sidecars-v1`, `medinter-population-mt-sources-v1` and `medinter-population-mt-metadata-v1`. The MT receipt identifies `public-mt-population-v2.sqlite`; the dossier receipt identifies `public-hcm39-clinvar-dossier-v2.sqlite`. Earlier partial databases remain preserved and are not advertised as complete resources. `medinter-population-mt-sidecars-v1/actual-public-source-controls.json` records the passed actual offline source controls and the public sidecars' SHA256 values; these are public reference hashes, never patient hashes. The compact `resource-ready-manifest.json` in that directory is the root-review entrypoint, with exact artifacts, source versions and unresolved gaps. It does not bind any stage.

Controls cover denominator missingness versus zero, multi-allele AC/homozygote bounds, population-field retention when dropping duplicate VEP consequences, BED-to-one-based coordinate boundaries, FASTA line boundaries, the standard CRC32C control, NCBI preamble/gzip validation, exact Content-Range, endpoint allowlisting, long public ClinVar fields and coverage identity/missingness. No control reads patient data or sends patient variants externally. Source-only readiness remains separate from root review and any future annotation binding.
