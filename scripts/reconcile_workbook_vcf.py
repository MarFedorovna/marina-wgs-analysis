import gzip
import json
import struct
from collections import defaultdict

from openpyxl import load_workbook

def key(chrom, pos, ref, alt):
    c = str(chrom).removeprefix("chr")
    return ("MT" if c in {"M", "MT"} else c, int(pos), str(ref).upper(), str(alt).upper())


def read_linear_index(path):
    with gzip.open(str(path) + ".tbi", "rb") as index:
        assert index.read(4) == b"TBI\x01"
        count, _, _, _, _, _, _, names_len = struct.unpack("<8i", index.read(32))
        names = index.read(names_len).rstrip(b"\x00").decode().split("\x00")
        assert len(names) == count
        offsets = {}
        for name in names:
            bins, = struct.unpack("<i", index.read(4))
            for _ in range(bins):
                _, chunks = struct.unpack("<Ii", index.read(8))
                index.read(16 * chunks)
            intervals, = struct.unpack("<i", index.read(4))
            offsets[name] = struct.unpack(f"<{intervals}Q", index.read(8 * intervals))
    return offsets


class BgzfLines:
    def __init__(self, handle):
        self.handle = handle
        self.block = b""
        self.offset = 0
        self.cursor = 0

    def read_block(self):
        self.offset = self.handle.tell()
        fixed = self.handle.read(12)
        if not fixed:
            self.block = b""
            return
        if len(fixed) != 12 or fixed[:4] != b"\x1f\x8b\x08\x04":
            raise ValueError("Invalid BGZF block")
        extra = self.handle.read(struct.unpack_from("<H", fixed, 10)[0])
        start = 0
        while start + 4 <= len(extra):
            tag, size = struct.unpack_from("<2sH", extra, start)
            if tag == b"BC" and size == 2:
                block_size = struct.unpack_from("<H", extra, start + 4)[0] + 1
                rest = self.handle.read(block_size - len(fixed) - len(extra))
                self.block = gzip.decompress(fixed + extra + rest)
                self.cursor = 0
                return
            start += size + 4
        raise ValueError("Missing BGZF block size")

    def seek(self, virtual):
        if self.offset != virtual >> 16 or not self.block:
            self.handle.seek(virtual >> 16)
            self.read_block()
        self.cursor = virtual & 0xFFFF

    def lines(self):
        fragments = []
        while self.block:
            end = self.block.find(b"\n", self.cursor)
            if end == -1:
                fragments.append(self.block[self.cursor:])
                self.read_block()
                continue
            part = self.block[self.cursor:end]
            self.cursor = end + 1
            yield b"".join([*fragments, part]) if fragments else part
            fragments.clear()


def indexed_records(path, wanted):
    names = read_linear_index(path)
    grouped = defaultdict(set)
    for k in wanted:
        grouped[k[0]].add(k[1])
    with path.open("rb") as handle:
        reader = BgzfLines(handle)
        for name, intervals in names.items():
            positions = grouped.get(key(name, 1, "", "")[0])
            if not positions or not intervals:
                continue
            starts = defaultdict(set)
            for pos in sorted(positions):
                window = min((pos - 1) >> 14, len(intervals) - 1)
                while window and not intervals[window]:
                    window -= 1
                starts[intervals[window] or next(x for x in intervals if x)].add(pos)
            last_record = None
            for offset, candidates in sorted(starts.items()):
                reader.seek(offset)
                end = max(candidates)
                for line in reader.lines():
                    chrom, position, _, ref, alts, _ = line.split(b"\t", 5)
                    if chrom.decode() != name:
                        break
                    pos = int(position)
                    if pos > end:
                        break
                    if pos not in candidates or line == last_record:
                        continue
                    ref = ref.decode()
                    if any(key(name, pos, ref, alt.decode()) in wanted for alt in alts.split(b",")):
                        last_record = line
                        yield line.decode("utf-8")


def vcf_header(path):
    data = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ,"):
                data["csq"] = line.split("Format: ", 1)[1].rstrip().removesuffix('">').split("|")
            if line.startswith("##WGS_EVIDENCE_SCHEMA="):
                data["evidence_schema"] = line.split('Data="', 1)[1].split('"', 1)[0]
            if line.startswith("#CHROM"):
                data["samples"] = line.rstrip().split("\t")[9:]
                break
    return data


def collect_targets(workbook_path):
    workbook = load_workbook(workbook_path, read_only=True, data_only=True)
    targets = {}
    for name in ("Все_прошедшие", "X_Y_все", "MT_все"):
        sheet = workbook[name]
        rows = iter(sheet.values)
        columns = next(rows)
        indexes = [columns.index(h) for h in ("Хромосома", "Позиция", "REF", "ALT")]
        targets[name] = {key(*(row[i] for i in indexes)) for row in rows if all(row[i] not in (None, '', '.') for i in indexes)}
    workbook.close()
    return list(targets.values())


def extract_vcf(path, wanted):
    records = {}
    for line in indexed_records(path, wanted):
        fields = line.split("\t", 8)
        for allele_number, alt in enumerate(fields[4].split(","), 1):
            k = key(fields[0], fields[1], fields[3], alt)
            if k not in wanted:
                continue
            value = {"id": fields[2], "qual": fields[5], "filter": fields[6], "info": fields[7], "allele_number": allele_number}
            if k in records and records[k] != value:
                raise ValueError(f"Conflicting indexed VCF records at {k}")
            records[k] = value
    return records


def write_extract(vcf_path, output_path, targets):
    auto, xy, mt = targets
    merged = extract_vcf(vcf_path, auto)
    print(f"Fully annotated VCF: {len(merged)}/{len(auto)} autosomal alleles", flush=True)
    if set(merged) != auto:
        raise ValueError("Fully annotated VCF did not account for all autosomal report alleles")
    with output_path.open("w", encoding="utf-8") as output:
        output.write(json.dumps({"header": vcf_header(vcf_path), "targets": [len(auto), len(xy), len(mt)]}) + "\n")
        for k, value in merged.items():
            output.write(json.dumps([k, value], ensure_ascii=False) + "\n")
