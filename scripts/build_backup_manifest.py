#!/usr/bin/env python3
"""Create a private SHA-256 manifest for a project snapshot.

The generated manifest can contain sensitive relative paths. Store it outside
Git in a directory ignored by .gitignore, for example local/manifests/.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--exclude", action="append", default=[".git", "__pycache__"])
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    out = args.out.resolve()
    excluded = set(args.exclude)
    rows = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        relative = path.relative_to(root)
        if path.resolve() == out or any(part in excluded for part in relative.parts):
            continue
        rows.append({
            "relative_path": relative.as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": digest(path),
        })
    payload = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "root_label": root.name,
        "file_count": len(rows),
        "files": rows,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(rows)} records to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

