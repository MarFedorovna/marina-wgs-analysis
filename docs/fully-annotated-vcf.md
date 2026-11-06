# Portable applied-evidence cohort VCF

`scripts/embed_vcf_evidence.py` adds the completed combined annotation master and actually applied public evidence to a new cohort VCF. It does not rerun VEP, change the primary annotation/report, assign ACMG criteria, or change ranking eligibility. The original stage 06/07 and combined master remain immutable controls. Tables and top lists produced by views v3 read the new VCF; comparison with the combined master verifies lossless extraction.

The fixed output directory is `/srv/marina/outputs/fully-annotated-vcf-v1`. It contains `cohort.fully-annotated.vcf.gz`, its `.tbi`, `manifest.json` and `workflow-receipt.json`. Completion requires workflow status `PASS`, stage `fully-annotated-vcf-v1`, and manifest status `PASS_FULLY_ANNOTATED_VCF` with matching artifact hashes. Preparation alone is not completion.

## Embedded schema and evidence

`INFO/WGS_EVIDENCE` is Number=1, Type=String. The codec is canonical JSON (sorted keys, compact separators, UTF-8), zlib, then URL-safe Base64 with padding omitted. Decoding is bounded to 64 MiB per packet and rejects truncated, concatenated and oversized streams. No field is truncated.

The packet shape is:

```json
{"v":1,"allele":["chr1",10,"A","G"],"scope":"cohort","sample":null,
 "rows":[{"i":1,"values":["all original and appended strings in schema order"]}],
 "evidence":{"population":[],"exome_AN":[],"allele_map":[],"SCV":[],"citations":[]},
 "gene_context_keys":["HGNC:1"]}
```

Each `i` is the one-based global combined-master row ordinal. `values` includes every original field and all appended `allref_*` fields, including full JSON strings and missing values. Missing values never become AF=0. Every source allele must have explicit master rows for all four samples, including placeholder rows for variants without consequence annotations. Duplicate VCF alleles, absent master alleles, zero-row alleles, unknown samples, ordinal gaps and repeated rows stop the job.

`##WGS_EVIDENCE_SCHEMA=<Version=1,Codec="json-utf8+zlib+base64url-nopad",Data="TOKEN">` embeds the ordered `original_fields`, `appended_fields`, `full_fields`, row count, original invariant digest, full-field digest and scope. The original invariant is SHA256 of Python `csv.writer(delimiter=TAB,lineterminator=LF)` output for the original ordered header and original field vectors, UTF-8. The full-field digest uses the same convention with every field. Original row ordinals allow exact reconstruction without an external TSV.

The applied evidence contains full exact-allele gnomAD v4.1.1 public39 CDS counts and ancestry/quality INFO; exome v4.1 AN position context; validated source AlleleID/VariationID maps; full applied versioned SCV records; and citation links. Citation presence does not establish article review, independent cases or clinical classification. The v4.1 exome AN source is not genome AN or sequencing depth. Original cache AF remains separately labelled.

Full applied HGNC gene-bin rows are deduplicated in `WGS_GENE_CONTEXT` header chunks. Their decoded object contains `v`, `scope`, an ordered `fields` array and `genes` mapping HGNC IDs to arrays of complete field vectors. These bins are gene context, never proof of depth at the variant locus.

Full compact public39 GTEx, pext and UniProt source objects are included in `WGS_PUBLIC_REFERENCE` chunks, with the pinned context index, source versions, raw UTF-8 content and original byte sizes/SHA256. This includes protein feature coordinates and sequences and pext regions; a source object does not establish an unverified transcript, genome, protein-domain or pext allele mapping. Existing mapping-status fields remain unchanged. This is the approved compact public panel, not a whole database export.

Both chunk types use:

```text
##WGS_PUBLIC_REFERENCE=<Version=1,Part=1,Parts=2,SHA256="SHA_OF_CONCATENATED_TOKEN",Data="CHUNK">
```

Concatenate `Data` in numeric `Part` order, verify the SHA256 of the concatenated ASCII token, then use the packet codec. Parts are at most 60,000 characters. All inflated schema/gene/public header objects together are limited to 64 MiB. `decode_header_bundle(lines,name)` verifies the chunk sequence, checksum and decompression bounds.

## Readable INFO fields

All supplemental readable String values use URL-percent encoding per comma-separated entry. Decode each entry with `urllib.parse.unquote`. The lossless packet is authoritative for complete values and source scope.

| INFO | Contents |
| --- | --- |
| `WGS_GENE`, `WGS_CONSEQUENCE` | Gene symbols and consequences from the full source rows |
| `WGS_CLINVAR`, `WGS_CLINVAR_REVIEW` | Original significance, source IDs/snapshot, review and explicit SCV map/review flags |
| `WGS_GNOMAD_EXOMES`, `WGS_GNOMAD_GENOMES` | Exact public39 source release, match status, AC, AN, AF, nhomalt and population FILTER; absent values stay missing |
| `WGS_GNOMAD_{EXOMES,GENOMES}_{AC,AN,AF,nhomalt}` | Eight Number=1 typed global scalars for unambiguous `EXACT_ALLELE` v4.1.1 matches only: AC/AN/nhomalt Integer, AF Float; missing, unmatched and outside-subset values are `.`; suitable for standard bcftools filtering |
| `WGS_CACHE_AF` | Original cache AF/max/status/version and AF source JSON, explicitly distinct from supplementary AC/AN |
| `WGS_TX_SCORES` | `sample|exact transcript|protein|source-field:value` for SIFT, PolyPhen, FATHMM, CADD and AlphaMissense; transcript versions and score families retained |
| `WGS_LOCAL_SCORES` | Sample and transcript-labelled original local-evidence JSON, including source score arrays; no isoform flattening or inferred reassignment |
| `WGS_TISSUE_STATUS` | Gene/transcript labels and actual GTEx/pext match statuses |
| `WGS_ANN_SCOPE`, `WGS_ANN_STATUS` | Actual availability and limits of the applied reference context |

The original CSQ, INFO order, QUAL, FILTER, FORMAT and all sample strings remain byte-identical. Only new INFO entries are appended. An originally missing INFO `.` is restored by the verification decoder. Original header lines, including CSQ/FORMAT definitions, are independently hashed and verified in sequence. Data record bytes restored after removing the new INFO tags must match the original data-record digest. Indexed seeks test the first available locus on every input contig.

The mitochondrial branch remains a separate appendix and separate MT VCF. Its distinct genotype, heteroplasmy, QC, NUMT and haplocheck schemas are not inserted into nuclear cohort rows. No private per-sample MT metadata is put in the public cohort header. Per-sample VCF packets are filtered to that sample's rows; the global public reference objects may be retained because they contain no cohort genotypes.

## Decoder and controller

Portable extraction needs only Python's standard library and the VCF:

```sh
python3 scripts/embed_vcf_evidence.py decode --vcf cohort.fully-annotated.vcf.gz --master reconstructed.tsv
```

The module exposes `encode_packet`, `decode_packet`, `encode_schema`, `decode_schema`, `schema_header`, `decode_header_bundle`, `select_sample`, `rebuild_readable` and `NAMESPACED_INFO_IDS` for the downstream views exporter. The standalone decoder uses a disk SQLite index to restore global ordinal order.

`scripts/continue_fully_annotated_vcf.py prepare` creates a private candidate only. It requires the actual reviewed combined-master producer rule and existing pinned pysam 0.24.1 environment. The reviewed run command is:

```sh
PYTHON scripts/continue_fully_annotated_vcf.py run --rule ROOT-reviewed-rule.json --reviewed-rule-sha256 EXACT_RULE_SHA --execute
```

Only ROOT changes the candidate review status to `ROOT_REVIEWED_EXACT_FULLY_ANNOTATED_VCF_RULE` after review. The controller waits for the combined-master PASS receipt, locks its private state, and writes status to `/srv/marina/private/fully-annotated-vcf-state-v1/status.json` and phase/count-only progress to `progress.json`. It does not start or retry upstream stages.

Limits are 2 CPUs, 4 GiB memory, 32 GiB compressed SQLite, 16 GiB VCF output, a 24-fold output/source compressed-byte inflation ceiling, 150 GiB remaining disk, 64 MiB record payload and 64 MiB total inflated supplemental header. Original deadline `2026-11-06T03:46:48.046025Z` and at most eight hours cannot be extended; a controller alarm and periodic indexing/verification/checksum guards enforce the deadline. Failures retain hidden staging artifacts and never promote partial success. Successful promotion is atomic and refuses to replace an existing output.

Synthetic controls cover filtered variants, missing and phased calls, multiple transcripts, empty consequence placeholders, explicit AN=0, full ancestry INFO, full SCV/citations, source preservation, decompression attacks, sample isolation and rejected incomplete masters. A separate actual-public-source control embeds the approved GTEx/pext/UniProt objects into a synthetic cohort BGZF VCF and verifies header/source hashes and Tabix reading. No patient job is needed for these controls.
