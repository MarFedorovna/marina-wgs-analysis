# Public LOVD reference preparation

This is a separate reference snapshot, **not an annotation already applied to the final VCF**.
The existing final VCF and its active download remain immutable. A future version requires
reviewed local allele matching and its own receipts; this collector has no VCF input.

## Access and rights

- [LOVD official FAQ](https://www.lovd.nl/3.0/faq) distinguishes free basic API retrieval,
  curator-controlled downloads, and licensed full case-level JSON exports. Data rights can
  differ by installation and submitter. Software/manual licences do not license all records.
- [Official manual, section 11.1](https://www.lovd.nl/3.0/docs/LOVD_manual_3.0.pdf)
  documents unauthenticated basic genes/variants retrieval and a maximum of five requests
  per second per server. This profile uses **one request per second**, never more than four.
  Unavailable entries are not retried for at least four hours.
- [Official public example](https://databases.lovd.nl/shared/api/rest.php/variants/IVD/0000000563?format=application/json)
  was retrieved once on 2026-11-06: 624 bytes. The actual response is one JSON object,
  with explicit `variant_genomic.hg19` and `variant_genomic.hg38` descriptions and
  versioned transcript HGVS. This differs from the older manual's flat position example.
- [HGMD public registration](https://www.hgmd.cf.ac.uk/docs/register.html) currently
  redirects to the [Cardiff service disruption page](https://outage.cardiff.ac.uk/service-disruption/).
  Indexed registration information describes academic registration and delayed public data;
  the current live service/terms cannot be verified while unavailable. No anonymous
  HGMD bulk/API access was established. No HGMD licence was supplied.
- [VarSome API packages](https://docs.varsome.com/en/varsome-api-packages-explained)
  describe chromosome-17-only Test access and subscription-based Standard access.
  No VarSome licence/token was supplied. Free website access is not a bulk API subscription.

HGMD and VarSome therefore remain `ACCESS_REQUIRED_NOT_APPLIED`. No registration,
subscription, payment, licence acceptance, cookie extraction, or patient API query is performed.

## Fixed source scope

The collector accepts only the existing, public 39-gene panel, including TNNT2 and MT-TI.
Its exact panel SHA256 is
`09c518f2659bd2443ee024d65dd83070922864304ed92f705271e250c227fedb`.
There are 78 basic requests: one gene metadata and one variants response per public gene.
Maximum response 16 MiB, total 256 MiB, duration 10 minutes, reserve 150 GiB, and the
original compute deadline 2026-11-06T03:46:48.046025Z. No pagination/completeness is inferred.
Redirects are refused; default TLS verification remains enabled. No API key or login is used.
Raw responses, source URL/hash, source build, observation ID, stable LOVD DBID, transcript
versions and raw assertions are retained. `Times_reported` is not independent case evidence.

All rows remain **unapplied**. Explicit `variant_genomic.hg38` single-nucleotide changes
are only `GRCH38_SNV_CANDIDATE_PENDING_LOCAL_FASTA_VALIDATION`. Build-37 descriptions,
position-only descriptions and unsupported complex HGVS remain `UNMAPPED_*`.
No transcript-to-genome guess, liftover, automatic classification or ACMG/PP5 is performed.

## Reviewable commands

Prepare a new source profile, without network access:

```sh
python3 collect_lovd_public_basic_reference.py prepare \
  --panel /srv/marina/resources/public-hcm-cds-v2/panel-manifest.json \
  --profile /srv/marina/private/lovd-public-basic-candidate-v1/profile-candidate.json \
  --output /srv/marina/resources/lovd-public-basic-public39-v1
```

Root must review the exact code, panel pin, scope, rights and limits. A reviewed copy sets
`review_status=ROOT_REVIEWED_PUBLIC_LOVD_BASIC_REFERENCE_PROFILE` and its SHA256 is
passed to `collect --reviewed-profile-sha256 HASH --execute`. Without `--execute`,
`collect --profile PROFILE` prints only the plan. The current deadline must not be extended
merely to fetch this optional source.

The collector's final source manifest is
`PUBLIC_LOVD_BASIC_REFERENCE_SNAPSHOT_NOT_APPLIED`; failed HTTP responses are explicitly
unavailable, never an empty evidence set. A separate future local application must require
build GRCh38, exact normalized CHROM/POS/REF/ALT and FASTA REF validation, source/version
hashes and provenance, and must preserve all original VCF and genotype fields.
