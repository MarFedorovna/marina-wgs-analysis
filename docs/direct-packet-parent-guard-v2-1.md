# Direct packet views: parent guard revision 2.1

Use `continue_packet_views_fast_v2_1.py` for the reviewed 24-worker report
continuation. Earlier controller packages are retained as history. The exporter
and its scientific ranking, packet projection and sample filtering are unchanged.

```sh
python3 scripts/continue_packet_views_fast_v2_1.py \
  --rule "$FAST_VIEWS_RULE" --rule-sha256 "$FAST_VIEWS_RULE_SHA256"
```

The default is status only. After review, append `--execute --watch` to wait and
produce the independent reports. The exact rule binds the code, original direct
producer, accepted input and public-control identities. No primary process is
started or stopped.

A POSIX real-time timer bounds the whole parent operation, including waiting,
worker dispatch, table concatenation, sample VCF assembly, artifact hashes and
atomic promotion. Its deadline is the earliest of the unchanged VM allowance,
the reviewed absolute deadline and the bounded watch duration. Parent address
space is limited to 2 GiB. Prior signal handler/timer and memory limits are restored
on success or failure; worker limits remain 1.5 GiB each. At most 24 workers run,
with at most 128 MiB queued source text plus one 64 MiB current batch.

The VM tests interrupt a real two-second parent sleep after 50 ms, reject an already
expired deadline, and verify signal/memory restoration after success and failure.
The separate 2-versus-24-worker proof compares all decoded output tables and every
sample packet. All transcripts remain in the final master; main top20 ranking
continues to use the unchanged frozen policy.
