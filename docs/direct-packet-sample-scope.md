# Direct packet sample scope

The direct producer adds a cohort `membership` object containing all four sample states. This is appropriate in the joint VCF and joint intersection tables. It must not be copied into a single-sample VCF.

`spliceai_vcf_codec_v2.py` preserves the reviewed first SpliceAI codec and changes only `select_sample`: it retains the selected sample's complete original/context/score rows and the closed public packet keys `v`, `allele`, `rows`, `evidence`, `gene_context_keys`, `scope`, `sample`. It drops cohort membership and unknown future private metadata. The own original GT/FT remains in the selected rows and the one-sample VCF FORMAT columns. Full public reference headers and exact-allele public evidence may remain.

The cohort producer still retains its complete membership. No new cohort uniqueness claim is inferred from a single-sample export. INFO IDs, field order, score definitions and portable encoding remain unchanged. The earlier codec source files and original primary workflow stay immutable. Tests explicitly inject another sample's states and unknown private annotations, and verify that neither survives the single-sample packet or readable tags.
