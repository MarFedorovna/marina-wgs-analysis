# Bounded automatic database resume

Run this command **on the server** from the reviewed checkout:

```bash
python3 scripts/supervise_databases.py --created-utc VM_CREATED_UTC
```

Replace the timestamp with the actual timezone-qualified VM creation time; it is never reset to the moment this command starts. Alternatively put `vm_created_utc` in the private server-budget settings. The supervisor starts the standalone database installer with 32 direct-server connections, monitors it, and allows at most two total attempts. A six-hour download timeout or a recognized transient network failure can trigger one automatic resume using the same verified range receipts. Checksum/identity/version/unsafe-path failures are not retried blindly.

The estimated VM allowance remains $40 at the selected $2.39/hour rate, measured from VM creation, within the user's $50 overall budget. At the allowance deadline the supervisor stops only its verified installer process group. A monotonic timer prevents a backward system-clock change from extending the allowance. These controls **do not terminate the VM or stop provider billing**, and they do not reserve future analysis time. Retain results and terminate the paid VM separately when appropriate.

To adopt an already running direct-server installer without restarting its unfinished ranges, add `--adopt-pid PID`. UID, exact installer/preparation path, process start time and separate process-group identity are checked before adopting or stopping a process. Calling/genotyping jobs are outside that group and are not controlled by this supervisor. Process state and total attempt count persist under `/srv/marina/private/databases-supervisor-status.json`; a process lock prevents duplicate supervisors. A repeated command cannot silently reset the allowance or attempt count.

For a background run:

```bash
nohup python3 -u scripts/supervise_databases.py --created-utc VM_CREATED_UTC \
  > /srv/marina/logs/annotation-supervisor.log 2>&1 < /dev/null &
```

Installer output remains in `/srv/marina/logs/annotation-databases.log`. `READY_VALIDATED` means the final resource manifest and all its recorded artifact hashes match the configured recipe. Partial downloads, retries, installed containers and individual small-database receipts do not establish full annotation readiness. No workstation download route is used.
