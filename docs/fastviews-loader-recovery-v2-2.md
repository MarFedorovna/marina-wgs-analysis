# Views controller loader recovery v2.2

The earlier dynamically loaded exporter was missing from `sys.modules`.
`ProcessPoolExecutor` could not pickle its `process_batch` function, so the
report failed before producing any batch outputs. Its partial directory contains
only the incomplete marker and is preserved.

The new controller adds `import sys` and registers the exporter by its module
specification name before `exec_module`. The frozen exporter, source authority,
ranking, genotype handling, packet pruning, memory limits and absolute deadline
are unchanged. The new reviewed rule uses its own state directory.

The Linux regression runs the actual controller loader and frozen exporter with
24 processes on 30 synthetic public alleles/240 rows. Only the production path
validator is replaced by a strict temporary-fixture validator. It checks the
lossless master, indexed individual VCF packets and completion receipt. The three
parent deadline and memory-limit tests run against the new controller as well.

This recovery restarts only the independent views producer after review. It
does not rerun VEP, change the completed cohort VCF or delete old partials.
