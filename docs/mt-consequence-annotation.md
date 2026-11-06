# Mitochondrial consequence annotation

`scripts/annotate_mitochondrial_consequences.py` annotates four completed, unshifted,
filtered Mutect2 mitochondrial callsets. It does not call variants again. This is
separate from nuclear X/Y calling and from the MT functional table exporter.

The inputs are the actual successful MT calling receipts and the completed
`mitochondrial-reference-context-v1` manifest. The latter supplies exact rCRS
gnomAD-MT and HelixMTdb reference context and the original NUMT, haplocheck and
contamination QC. Original files remain unchanged.

## Scientific contract

- GRCh38 `chrM`, 16,569 bases; the selected analysis sequence must match the pinned
  rCRS sequence. The public control compares it to RefSeq NC_012920.1.
- Offline VEP 115 with verified MT cache data, global cache metadata and exact-allele
  ClinVar custom annotation. It does not label the entire human cache as newly applied.
- Vertebrate mitochondrial genetic code 2. Actual controls distinguish it from the
  nuclear code: rCRS 3399A>G changes ATA to ATG and retains Met; 3564A>G changes TGA
  to TGG and retains Trp. Both produce synonymous MT-ND1 consequences in the pinned
  runtime. See [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG2).
- Every original record, REF/ALT vector, QUAL, FILTER, INFO token, FORMAT field and
  raw sample string is preserved in the final VCF. AF, HL, AD and DP remain source
  values. No nuclear GQ threshold, diploid carrier rule or calling-ploidy change is applied.
- Multiallelic source records remain intact. `ALLELE_NUM` identifies a consequence's
  ALT when available; population context retains its original exact-match/unmapped
  status. Missing frequency is not zero, and Helix arrays are not flattened.
- AlphaMissense and SpliceAI are explicitly not applied as MT models. Consequences
  and ClinVar assertions are research evidence, without automatic ACMG classification.

## VCF and table outputs

The new output root is `mitochondrial-consequence-annotation-v1`. Each sample has:

| Artifact | Meaning |
|---|---|
| `mitochondria.fully-consequence-annotated.vcf.gz` and `.tbi` | Original calls plus CSQ and complete applied evidence |
| `MT-functional-all-transcripts.tsv.gz` | Every consequence row, original heteroplasmy vectors and frequency/QC context |
| `master.from-vcf.tsv.gz` | Full field vectors reconstructed from the persisted final VCF |
| `source-QC/` | Copies of that sample's original QC records |

`WGS_MT_CONSEQUENCE_EVIDENCE` stores complete string field vectors and row ordinals.
`WGS_MT_CONSEQUENCE_SCHEMA` declares their order. Encoding is UTF-8 JSON, zlib and
unpadded base64url, using the pinned bounded `encode_packet`/`decode_packet` functions
from `embed_vcf_evidence.py`. Each plaintext packet is limited to 64 MiB; overflow
fails instead of truncating evidence. This MT namespace does not claim nuclear
core33 or splice-model annotation. Roundtrip checks read the actual BGZF file and
compare every reconstructed row and original record field.

The root manifest status is `PASS_MT_CONSEQUENCE_ANNOTATION`. Its actual completion
receipt has `status=PASS`, `stage=mitochondrial-consequence-annotation-v1` and artifact
pins. A directory or a progress status alone is not completion.

## Commands

Variables below refer to existing reviewed private inputs and a new output. Expired profiles are historical; a new run needs a newly reviewed deadline within the code limit.

```sh
python scripts/annotate_mitochondrial_consequences.py prepare \
  --helper scripts/annotate_supplementary_callsets.py \
  --base-config "$BASE_CONFIG" \
  --mt-context "$MT_CONTEXT_MANIFEST" \
  --codec scripts/embed_vcf_evidence.py \
  --output "$NEW_MT_CONSEQUENCE_OUTPUT" \
  --profile "$NEW_MT_PROFILE_CANDIDATE" \
  --deadline-utc "$REVIEWED_DEADLINE_UTC"
```

Preparation starts no annotation. A separately reviewed profile sets
`root_reviewed=true`, `status=ROOT_REVIEWED_MT_CONSEQUENCE_ANNOTATION` and pins the
actual `PASS_ACTUAL_MT_CODE2_OFFLINE_CONTROL` receipt for the same producer and VEP
image. Execution uses its exact file SHA256:

```sh
python scripts/annotate_mitochondrial_consequences.py run \
  --helper scripts/annotate_supplementary_callsets.py \
  --profile "$REVIEWED_PROFILE" --reviewed-profile-sha256 "$REVIEWED_PROFILE_SHA256" --execute
```

Limits are two CPUs, 8 GiB address space, 150 GiB disk reserve and the reviewed
absolute deadline, bounded by this run's 2026-11-06 06:25 UTC cutoff. Containers use
an immutable image, no network and read-only reference/cache mounts.

To organize individual and shared functional types, use the separate
[MT functional table exporter](mt-functional-tables.md). It reads these completed
VCF packets and verifies every vector against the original consequence tables.
