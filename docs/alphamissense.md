# AlphaMissense: installed public scores and strict offline matching

AlphaMissense is a mandatory table/report source by the user's explicit requirement. It remains separate from the uncertain fourth predictor in the original four-tool discussion and is not counted as another independent ACMG evidence item.

Install directly on the server:

```bash
python3 scripts/prepare_alphamissense.py --execute
```

The default lock is `config/alphamissense-resources.lock.json`. It pins Google DeepMind's public GCS object generations for both GRCh38 files, whose byte counts and publisher MD5 match [Zenodo record 8360242 v2](https://zenodo.org/records/8360242): canonical 642,961,469 bytes and isoforms 1,177,361,934 bytes. The updated licensing README is retained. Current [official project licensing](https://github.com/google-deepmind/alphamissense#alphamissense-predictions-license) permits predictions under CC BY 4.0; the older 2023 Zenodo README retains historical wording. The installer records the current licensing source and exact downloaded file hashes. Model weights are not released and are not used here; these are published precomputed predictions.

The installer builds local tabix indexes, checks schema and build, performs a public first-record lookup for both files, records tabix version/binary hash, and atomically publishes a separate supplement. It requires the 150 GiB reserve and preserves failed partials. This does not modify the running VEP base/v2 recipe or its download.

READY location:

```text
/srv/marina/resources/annotation-supplements/alphamissense-2023-v2-grch38-v3/resource-manifest.json
```

`annotation.alphamissense_manifest` can select another reviewed supplement. The default annotation stage requires this resource; it is queried only through local tabix against locally generated positions. No patient coordinates are sent to a remote service. Canonical and isoform source versions are both recorded.

## Identity and interpretation

The source uses GRCh38.p13 and GENCODE v32 transcripts, while VEP 115 uses newer transcript models. Matching requires an exact normalized genomic allele, the **full Ensembl transcript ID including version**, and the same reference/position/alternate protein substitution. Canonical rows additionally require the source UniProt accession to occur in the current annotation. The isoform file does not contain UniProt accessions; that missing source field is preserved, and exact transcript/protein-change matching is still required.

A transcript version, protein substitution or accession mismatch does not assign the score to the current transcript. The source's allele-level rows remain in `AlphaMissense_allele_predictions_json` for manual review. Conflicting exact predictions also remain unassigned. Non-missense transcripts are N/A; absent predictions are missing, not benign. Scores are never broadcast to every transcript at a genomic position.

Long, main, supplementary and top tables contain separate score, provider model label, matching status, source transcript/version, protein substitution/accession, dataset/release and raw source-context fields. A compatible score can contribute a separate, explicit research-ranking component (default threshold 0.564). It does not change genotype quality, prove mechanism, count in the original 3-of-4 vote, or assign PVS1/ACMG/GoF/LoF. The model labels `likely_pathogenic`, `likely_benign`, and `ambiguous` describe the predictor; they are **not clinical ACMG classifications**.

Validation on the selected server: both full files passed publisher MD5/size and local SHA256 recording; both indexes passed local probes. An additional offline test reproduced 20 public canonical/isoform score rows and rejected a deliberately changed transcript version. These tests used public variants and synthetic control genotypes, not patient data. This proves the AlphaMissense local adapter; the separate full VEP/ClinVar control test still depends on the VEP cache becoming ready.

Citation: [Cheng et al., Science 2023](https://doi.org/10.1126/science.adg7492). The [official code repository](https://github.com/google-deepmind/alphamissense) distinguishes implementation code, precomputed predictions and unavailable trained weights.
