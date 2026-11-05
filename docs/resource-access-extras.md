# Additional public reference access, 5 November 2026

These source downloads and isolated runtime preparation followed the user's explicit noncommercial research confirmation. The original download receipts remain intact. No patient alleles, genotypes, or clinical records were sent to public APIs. No accounts were created, subscriptions purchased, or provider contracts accepted.

| Resource | Actual result at 18:07 UTC | Evidence on VM |
| --- | --- | --- |
| GeneReviews | Official 640,489,521-byte archive downloaded intact; full gzip integrity and tar member safety validated; 2,940 members. Available as an unmodified noncommercial research reference, with attribution. | `/srv/marina/resources/genereviews-source-only-v1/research-reference-manifest.json` |
| SpliceAI 1.3.1 | Official 16,676,295-byte wheel, publisher SHA256 and ZIP CRC verified; five model weights. Python 3.10 CPU runtime installed in an isolated venv, all dependency wheels individually verified against PyPI metadata, offline hash-locked install and `pip check` passed. | `/srv/marina/resources/spliceai-research-runtime-v1/status.json`, `wheel-manifest.json`, `requirements.lock.txt` |
| SpliceAI synthetic control | Five models loaded; real probability output `(1,101,3)` was finite, normalized and in range; valid synthetic reference accepted and deliberately wrong reference rejected; DS/DP fields emitted. This is a runtime smoke control, not clinical validation. | `/srv/marina/resources/spliceai-research-runtime-v1/synthetic-control-receipt.json` |
| Full RefSeq VEP115 GRCh38 cache | Real queue PID 343859 waits for the primary Ensembl VEP resource manifest to become PASS. Official archive HEAD now succeeds: 24,971,548,035 bytes. Publisher CHECKSUMS entry is BSD sum 00476, 24,386,278 blocks. No second large transfer starts while the primary download is incomplete. | `/srv/marina/resources/refseq-vep115-source-only-v1/status.json` |
| OMIM public identifier mapping | Actual `mim2gene.txt` downloaded, 991,361 bytes. This maps identifiers and does not contain full OMIM descriptions or the licensed gene/phenotype database. | `/srv/marina/resources/reference-api-public-catalogs-v1/omim_mim2gene-receipt.json` |
| OMIM full content/API | ACCESS_REQUIRED. Profile expects an account holder's authorized key in `OMIM_API_KEY`; no key was obtained or used. | `config/reference-api-profiles.research-v1.json` |
| LOVD | Fixed public MYH7 gene endpoint returned an 851-byte Atom XML record. The generic gene catalog endpoint returned no usable body. MYH7 reports `NC_000014.8`, transcript `NM_000257.2`; these source coordinates are explicitly not accepted as GRCh38 variant coordinates. No variant or case data was requested. | `/srv/marina/resources/reference-api-public-catalogs-v1/lovd_myh7_gene-receipt.json` |
| VarSome API | ACCESS_REQUIRED. Profile expects `VARSOME_API_KEY`; test package is limited to chromosome 17. Standard API requires a paid subscription and has metered transfer. No account/token/payment flow was attempted. | `config/reference-api-profiles.research-v1.json` |

## Availability and usage boundaries

The earlier GeneReviews and SpliceAI download obstruction was workflow/access classification, not a genuine account gate: both have official direct public downloads. RefSeq's previously observed connection refusal was transient; the VM returned HTTP 200 on 5 November at 18:03 UTC. Authentication is genuinely required for full OMIM and VarSome API access. LOVD public gene metadata works without an account, while broader data export and reuse depend on curators and submitter licenses.

GeneReviews is owned by the University of Washington, Seattle, © 1993–2026. The [official copyright notice](https://www.ncbi.nlm.nih.gov/sites/books/NBK138602/) grants specified noncommercial research reproduction/distribution rights with source and copyright credit and an original-material link when republished online; further modification is prohibited. Retain attribution to [GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1116/). The unchanged archive is from the [official NLM bulk endpoint](https://ftp.ncbi.nlm.nih.gov/pub/litarch/ca/84/gene_NBK1116.tar.gz). Its inspected official directory publishes no checksum sidecar, so its receipt correctly labels the SHA256 as a local integrity digest, not a publisher checksum.

The [current official SpliceAI repository](https://github.com/Illumina/SpliceAI) identifies code and model licensing separately; the historical PyPI 1.3.1 package metadata also retains its older GPLv3 label. Retain both provenance records and applicable attribution. The prepared runtime uses the official versioned PyPI wheel and does not silently equate its bundled GENCODE v24 annotations to VEP115/GENCODE49. A separate transcript/reference integration and suitable public biological controls are required before using this runtime in the patient analysis pipeline.

RefSeq queue disk guards retain 150 GiB for analysis, a conservative **128 GiB allowance** for unmeasured expanded RefSeq cache, the remaining archive bytes, and 10 GiB headroom. This is a planning ceiling, not a measured expansion size. It enforces the unchanged absolute deadline `2026-11-06T03:46:48.046025Z`, at most four network attempts and a two-hour transfer cap after primary readiness. Successful source-only download will still require archive/expanded-size validation, isolated extraction, and matching VEP115 `--refseq` public controls; it does not install or overwrite the existing Ensembl cache.

The [LOVD FAQ](https://www.lovd.nl/3.0/faq) distinguishes public gene/variant API records from case-level downloads and explains curator and submitter licensing. Its [documentation](https://databases.lovd.nl/shared/docs/) directs automation to APIs. The cached gene record is bibliographic/gene metadata, never an assertion that a patient's allele is pathogenic.

The [VarSome package documentation](https://docs.varsome.com/en/varsome-api-packages-explained) describes chromosome-17-only test queries and paid Standard access. For [OMIM Downloads](https://omim.org/downloads), the account holder must obtain their own approved key and review allowed usage. The provided profiles name secret environment variables but contain no credentials and do not enable authenticated requests.

## Read-only status and reproducibility

Use the already deployed code snapshot on the VM:

```sh
python3 /srv/marina/private/reference-access-extras-v1/scripts/reference_access_snapshot.py
python3 -m unittest discover -s /srv/marina/private/reference-access-extras-v1/tests -p test_reference_access_extras.py -v
```

Six safety tests passed on VM Python 3.10: exact resume Content-Range, symlink ancestor rejection, model archive count/CRC, unsafe tar link rejection, primary manifest compatibility, and disk reserve arithmetic. The Windows run passed five with the host symlink test skipped. Source scripts use distinct resource directories and locks; leave the existing RefSeq queue running instead of starting another transfer. Source receipts and runtime control files distinguish downloaded, runtime-ready, queued and access-required states.

This delta adds only its new source files; it does not modify any pre-existing pipeline stages, active private configurations, frozen bootstrap files, or patient outputs.
