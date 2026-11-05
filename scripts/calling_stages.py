#!/usr/bin/env python3
"""Private germline call/joint/filter stages; no clinical interpretation.

References for the command/filter policy:
https://docs.nvidia.com/clara/parabricks/tool-reference/tools/haplotypecaller
https://gatk.broadinstitute.org/hc/en-us/articles/9570489472411-GenotypeGVCFs
https://gatk.broadinstitute.org/hc/en-us/articles/4404607512475-GenomicsDBImport
https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
https://samtools.github.io/bcftools/bcftools.html
"""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
from collections import Counter

GATK_VERSION = '4.3.0.0'

IDENTIFIER = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,99}\Z')
DEFAULT_FILTERS = {
    'min_dp': 10, 'min_gq': 20, 'min_qual': 30,
    'snp': {'QD': ['lt', 2], 'FS': ['gt', 60], 'SOR': ['gt', 3], 'MQ': ['lt', 40],
            'MQRankSum': ['lt', -12.5], 'ReadPosRankSum': ['lt', -8]},
    'indel': {'QD': ['lt', 2], 'FS': ['gt', 200], 'ReadPosRankSum': ['lt', -20]},
    'missing_annotations': 'flag',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def number(value):
    try:
        parsed = float(value)
        return parsed if math.isfinite(parsed) else None
    except (TypeError, ValueError):
        return None


def filter_policy(overrides=None):
    policy = json.loads(json.dumps(DEFAULT_FILTERS))
    for key, value in (overrides or {}).items():
        if key == 'site_qual':
            policy['min_qual'] = value
            continue
        aliases = {'min_qd': ('QD', 'lt'), 'max_fs': ('FS', 'gt'), 'max_sor': ('SOR', 'gt'),
                   'min_mq': ('MQ', 'lt'), 'min_mq_rank_sum': ('MQRankSum', 'lt'),
                   'min_read_pos_rank_sum': ('ReadPosRankSum', 'lt')}
        if key in aliases:
            annotation, operator = aliases[key]
            policy['snp'][annotation] = [operator, value]
            if key == 'min_qd':
                policy['indel']['QD'] = [operator, value]
            continue
        indel_aliases = {'indel_max_fs': ('FS', 'gt'), 'indel_min_read_pos_rank_sum': ('ReadPosRankSum', 'lt')}
        if key in indel_aliases:
            annotation, operator = indel_aliases[key]
            policy['indel'][annotation] = [operator, value]
            continue
        require(key in policy, 'Unknown filter-policy field.')
        if key in ('snp', 'indel'):
            policy[key].update(value)
        else:
            policy[key] = value
    require(policy['missing_annotations'] in ('flag', 'ignore'), 'Invalid missing-annotation policy.')
    require(all(number(policy[k]) is not None and float(policy[k]) >= 0
                for k in ('min_dp', 'min_gq', 'min_qual')), 'Filter thresholds must be nonnegative numbers.')
    for kind in ('snp', 'indel'):
        for annotation, rule in policy[kind].items():
            require(IDENTIFIER.fullmatch(annotation) and isinstance(rule, list) and len(rule) == 2
                    and rule[0] in ('lt', 'gt') and number(rule[1]) is not None, 'Invalid annotation threshold.')
    return policy


def genotype_flags(keys, values, policy):
    fields = dict(zip(keys, values))
    flags = [f for f in fields.get('FT', '').split(';') if f not in ('', '.', 'PASS')]
    gt = fields.get('GT', '.')
    alleles = re.split(r'[/|]', gt)
    if not alleles or any(a == '.' or not a.isdigit() for a in alleles):
        flags.append('NoCall')
    dp, gq = number(fields.get('DP')), number(fields.get('GQ'))
    if dp is None or dp < float(policy['min_dp']):
        flags.append('LowDP')
    if gq is None or gq < float(policy['min_gq']):
        flags.append('LowGQ')
    return list(dict.fromkeys(flags)), any(a.isdigit() and int(a) > 0 for a in alleles)


def soft_filter_vcf(source, tagged, selected, overrides=None):
    """Preserve all normalized rows/GT values; add site FILTER and individual FT.

    The selected subset retains site-PASS records with at least one non-reference
    genotype passing FT. Other genotypes remain present and explicitly tagged;
    consumers must respect FT. Neither GT masking nor row deletion occurs in tagged.
    Missing DP/GQ fail genotype confidence. Missing site metrics receive a separate
    flag by default, rather than silently passing unavailable annotations.
    """
    policy = filter_policy(overrides)
    counts = Counter()
    existing_filters, existing_ft = set(), False
    # A header is small; data records are streamed and never accumulated in memory.
    opener = gzip.open if str(source).endswith('.gz') else open
    with opener(source, 'rt') as stream:
        headers = []
        for line in stream:
            headers.append(line)
            if line.startswith('##FILTER=<ID='):
                existing_filters.add(line.split('ID=', 1)[1].split(',', 1)[0])
            if line.startswith('##FORMAT=<ID=FT,'):
                existing_ft = True
            if line.startswith('#CHROM\t'):
                break
        require(headers and headers[-1].startswith('#CHROM\t'), 'VCF lacks a valid sample header.')
        header_fields = headers[-1].rstrip('\n').split('\t')
        require(len(header_fields) >= 10, 'A genotyped VCF with samples is required.')
        sample_count = len(header_fields) - 9
        definitions = {
            'LowQUAL': 'QUAL is missing/nonfinite or below the configured threshold',
            'MissingAnnotation': 'A configured site annotation is unavailable/nonfinite',
            'UnsupportedAllele': 'Symbolic, spanning-deletion, or other unclassified allele; retained for review',
        }
        for kind in ('snp', 'indel'):
            for annotation, rule in policy[kind].items():
                definitions['Fail_' + annotation] = 'Fails configured ' + annotation + ' threshold; see filtering-policy.json'
        additions = [f'##FILTER=<ID={name},Description="{description}">\n'
                     for name, description in definitions.items() if name not in existing_filters]
        if not existing_ft:
            additions.append('##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filters: LowDP includes missing DP; LowGQ includes missing GQ; NoCall includes partial/no GT. GT is preserved.">\n')
        additions.append('##marina_soft_filter_policy=' + json.dumps(policy, separators=(',', ':')) + '\n')
        with Path(tagged).open('x') as all_rows, Path(selected).open('x') as chosen:
            for output in (all_rows, chosen):
                output.writelines(headers[:-1] + additions + headers[-1:])
            for line in stream:
                require(not line.startswith('#'), 'Unexpected VCF header within records.')
                fields = line.rstrip('\n').split('\t')
                require(len(fields) == 9 + sample_count, 'Malformed VCF sample columns.')
                flags = [f for f in fields[6].split(';') if f not in ('', '.', 'PASS')]
                quality = number(fields[5])
                if quality is None or quality < float(policy['min_qual']):
                    flags.append('LowQUAL')
                info = dict(piece.split('=', 1) if '=' in piece else (piece, True)
                            for piece in fields[7].split(';') if piece != '.')
                alts = fields[4].split(',')
                if any(not re.fullmatch('[ACGTNacgtn]+', a) for a in [fields[3], *alts]):
                    flags.append('UnsupportedAllele')
                    kind = None
                else:
                    if len(fields[3]) == 1 and all(len(a) == 1 for a in alts):
                        kind = 'snp'
                    elif any(len(a) != len(fields[3]) for a in alts):
                        kind = 'indel'
                    else:
                        kind = None
                        flags.append('UnsupportedAllele')
                if kind:
                    for annotation, (operator, threshold) in policy[kind].items():
                        value = number(info.get(annotation))
                        if value is None:
                            counts['missing_' + annotation] += 1
                            if policy['missing_annotations'] == 'flag':
                                flags.append('MissingAnnotation')
                        elif (operator == 'lt' and value < float(threshold)) or (operator == 'gt' and value > float(threshold)):
                            flags.append('Fail_' + annotation)
                flags = list(dict.fromkeys(flags))
                fields[6] = ';'.join(flags) if flags else 'PASS'
                keys = fields[8].split(':')
                require(len(set(keys)) == len(keys), 'Duplicate FORMAT keys.')
                original_keys = keys[:]
                if 'FT' not in keys:
                    keys.append('FT')
                passing_alt = False
                for index in range(9, len(fields)):
                    values = fields[index].split(':')
                    require(len(values) <= len(original_keys), 'Malformed genotype FORMAT values.')
                    values += ['.'] * (len(original_keys) - len(values))
                    gt_flags, has_alt = genotype_flags(original_keys, values, policy)
                    if len(keys) > len(original_keys):
                        values.append('PASS')
                    values[keys.index('FT')] = ';'.join(gt_flags) if gt_flags else 'PASS'
                    fields[index] = ':'.join(values)
                    for flag in gt_flags:
                        counts['genotype_' + flag] += 1
                    passing_alt |= not gt_flags and has_alt
                fields[8] = ':'.join(keys)
                text = '\t'.join(fields) + '\n'
                all_rows.write(text)
                counts['all_records'] += 1
                for flag in flags:
                    counts['site_' + flag] += 1
                if not flags and passing_alt:
                    chosen.write(text)
                    counts['selected_records'] += 1
    return {'policy': policy, 'counts': dict(counts), 'samples': sample_count,
            'selection': 'site PASS and at least one nonreference FT-PASS genotype; other GT/FT preserved',
            'raw_records_preserved_in_separate_input': True, 'normalized_records_all_preserved': True,
            'vqsr_applied': False}


def interval_rows(path, lengths):
    """Return checked 1-based closed intervals; never infer ploidy from contig names."""
    result = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith(('#', '@', 'track ', 'browser ')):
            continue
        fields = line.split()
        if len(fields) >= 3:
            contig, start, end = fields[:3]
            start, end = int(start), int(end)
            if Path(path).suffix.lower() == '.bed':
                start += 1
        else:
            match = re.fullmatch(r'([^:]+)(?::([0-9,]+)-([0-9,]+))?', fields[0])
            require(match, 'Unsupported interval syntax.')
            contig = match[1]
            start = int(match[2].replace(',', '')) if match[2] else 1
            end = int(match[3].replace(',', '')) if match[3] else lengths.get(contig, 0)
        require(contig in lengths and 1 <= start <= end <= lengths[contig], 'Interval is outside the selected reference.')
        result.append((contig, start, end))
    require(result, 'Empty call intervals.')
    return result


def context(cfg):
    reference = Path(cfg['reference']['fasta']).resolve(strict=True)
    require(sha256(reference) == cfg['reference']['sha256'], 'Reference SHA256 mismatch.')
    lengths = {line.split('\t')[0]: int(line.split('\t')[1])
               for line in Path(str(reference) + '.fai').read_text().splitlines()}
    require(reference.with_suffix('.dict').is_file(), 'Reference dictionary is required.')
    samples = read_tsv(cfg['samples_manifest'])
    ids = [s['sample_id'] for s in samples]
    require(ids and len(set(ids)) == len(ids) and all(IDENTIFIER.fullmatch(s) for s in ids), 'Invalid or duplicate sample IDs.')
    callsets = cfg.get('callsets', [])
    require(callsets and len({c['name'] for c in callsets}) == len(callsets), 'Provide unique explicit callsets.')
    seen = []
    metadata = []
    for callset in callsets:
        require(IDENTIFIER.fullmatch(callset['name']), 'Invalid callset name.')
        selected = callset.get('samples', ids)
        require(selected and len(set(selected)) == len(selected) and set(selected).issubset(ids), 'Invalid callset sample selection.')
        intervals = Path(callset['intervals']).resolve(strict=True)
        require(sha256(intervals) == callset['intervals_sha256'], 'Call intervals SHA256 mismatch.')
        spans = interval_rows(intervals, lengths)
        require(not any(c == other_c and start <= other_end and other_start <= end
                        for c, start, end in spans for other_c, other_start, other_end in seen),
                'Callsets overlap; use per_sample_ploidy within one interval callset instead.')
        seen.extend(spans)
        contigs = sorted({r[0] for r in spans}, key=lambda x: list(lengths).index(x))
        require(not any(c in ('MT', 'M', 'chrMT', 'chrM') for c in contigs),
                'Mitochondrial calling requires a separate reviewed policy; it is not implemented here.')
        sex = any(c in ('X', 'Y', 'chrX', 'chrY') for c in contigs)
        scope = callset.get('scope', 'autosomes')
        if scope == 'autosomes':
            require(all(re.fullmatch(r'(?:chr)?(?:[1-9]|1[0-9]|2[0-2])', c) for c in contigs),
                    'Autosome scope contains non-autosomal contigs.')
        if sex:
            require(scope != 'autosomes' and all(s in callset.get('per_sample_ploidy', {}) for s in selected),
                    'X/Y require explicit per-sample ploidy and a separate callset scope; no sex is inferred.')
        for sample in selected:
            ploidy = callset.get('per_sample_ploidy', {}).get(sample, callset.get('ploidy', 2))
            require(str(ploidy) in ('1', '2'), 'Only explicit haploid/diploid call partitions are supported.')
        metadata.append({'callset': callset['name'], 'scope': scope, 'covered_contigs': contigs,
                         'excluded_contigs': list(dict.fromkeys(callset.get('excluded_contigs', []) + [c for c in lengths if c not in contigs])),
                         'intervals_sha256': callset['intervals_sha256'],
                         'included_samples': selected, 'excluded_samples': [sid for sid in ids if sid not in selected],
                         'full_whole_WGS_claimed': False, 'sex_inferred': False, 'MT_CNV_analysis': False})
    return reference, samples, callsets, metadata


def mounts(*paths):
    """Same absolute host/container names keep private manifests reproducible."""
    # Specific writable stage mount must not be nested below a read-only mount.
    unique = {}
    for path, readonly in paths:
        path = Path(path).resolve()
        key = str(path)
        unique[key] = {'source': key, 'target': key, 'readonly': bool(readonly)}
    return list(unique.values())


def checked_vcf(runtime, cfg, path, sample_ids, logname):
    tool = cfg['tools'].get('bcftools', 'bcftools')
    runtime.run([tool, 'index', '--force', '--tbi', str(path)], logname + '-index')
    actual = runtime.run([tool, 'query', '-l', str(path)], logname + '-samples', capture=True).splitlines()
    if not runtime.dry_run:
        require(actual == list(sample_ids), 'VCF sample order/identity differs from the planned cohort.')
    runtime.run([tool, 'index', '--stats', str(path)], logname + '-index-check', capture=True)


def expand_sample_columns(source, destination, sample_ids):
    """Absent scoped samples remain missing; never invent reference GT/ploidy."""
    opener = gzip.open if str(source).endswith('.gz') else open
    with opener(source, 'rt') as stream, Path(destination).open('x') as output:
        original = None
        for line in stream:
            if line.startswith('#CHROM\t'):
                fields = line.rstrip('\n').split('\t')
                original = fields[9:]
                require(len(set(original)) == len(original) and set(original).issubset(sample_ids),
                        'Cannot expand a VCF with unknown/duplicate sample columns.')
                excluded = [sid for sid in sample_ids if sid not in original]
                output.write('##marina_scope_excluded_samples=' + json.dumps(excluded) + '\n')
                output.write('\t'.join(fields[:9] + list(sample_ids)) + '\n')
            elif line.startswith('#'):
                output.write(line)
            else:
                require(original is not None, 'Missing VCF sample header.')
                fields = line.rstrip('\n').split('\t')
                require(len(fields) == 9 + len(original), 'Malformed VCF sample columns.')
                values = dict(zip(original, fields[9:]))
                unknown = ':'.join('.' for _ in fields[8].split(':'))
                output.write('\t'.join(fields[:9] + [values.get(sid, unknown) for sid in sample_ids]) + '\n')


def alignment_rows(cfg, ids):
    rows = read_tsv(cfg['alignment_manifest'])
    require(rows and all(set(('sample_id', 'bam', 'receipt')).issubset(row) for row in rows), 'Invalid alignment manifest.')
    require(len({row['sample_id'] for row in rows}) == len(rows), 'Duplicate alignment sample IDs.')
    selected = {row['sample_id']: row for row in rows}
    require(set(ids).issubset(selected), 'Missing alignment for a requested sample.')
    for sid in ids:
        row = selected[sid]
        bam = Path(row['bam']).resolve(strict=True)
        receipt_path = Path(row['receipt']).resolve(strict=True)
        receipt = json.loads(receipt_path.read_text())
        require(receipt.get('sample') == sid and bam.parent == receipt_path.parent,
                'BAM does not belong to the original sample receipt directory.')
        require(receipt.get('config', {}).get('reference', {}).get('sha256') == cfg['reference']['sha256'],
                'BAM source receipt uses a different reference; refusing reuse.')
        outputs = receipt.get('outputs', {})
        require(bam.name in outputs and sha256(bam) == outputs[bam.name], 'BAM source receipt SHA256 mismatch.')
        bai = Path(str(bam) + '.bai')
        require(bai.name in outputs and bai.is_file() and sha256(bai) == outputs[bai.name], 'BAM index source receipt SHA256 mismatch.')
        row.update(bam=str(bam), receipt=str(receipt_path), bai=str(bai))
    return selected


def call_stage(cfg, workdir, runtime):
    reference, samples, callsets, scope = context(cfg)
    ids = [s['sample_id'] for s in samples]
    requested = [sid for sid in ids if any(sid in c.get('samples', ids) for c in callsets)]
    alignments = alignment_rows(cfg, requested)
    resources = cfg.get('resources', {})
    cpus = int(resources.get('call_cpus', 24))
    memory = int(resources.get('call_memory_gib', 160))
    require(1 <= cpus <= 30 and 32 <= memory <= 180, 'Invalid HaplotypeCaller CPU/memory limits.')
    rows = []
    for callset in callsets:
        for sid in callset.get('samples', ids):
            alignment = alignments[sid]
            ploidy = int(callset.get('per_sample_ploidy', {}).get(sid, callset.get('ploidy', 2)))
            parameters = {'sample': sid, 'callset': callset, 'ploidy': ploidy,
                          'image': cfg['tools']['parabricks_image'], 'reference_sha256': cfg['reference']['sha256'],
                          'source_receipt_sha256': sha256(alignment['receipt'])}
            def produce(directory, sid=sid, alignment=alignment, callset=callset, ploidy=ploidy):
                directory = Path(directory)
                if not runtime.dry_run:
                    (directory / 'tmp').mkdir()
                target = directory / 'variants.g.vcf.gz'
                args = ['pbrun', 'haplotypecaller', '--ref', str(reference), '--in-bam', alignment['bam'],
                        '--gvcf', '--out-variants', str(target), '--num-gpus', '1', '--ploidy', str(ploidy),
                        '--interval-file', str(Path(callset['intervals']).resolve()), '--num-htvc-threads', '8',
                        '--tmp-dir', str(directory / 'tmp'), '--logfile', str(directory / 'parabricks.log')]
                command = runtime.docker(cfg['tools']['parabricks_image'], args,
                    mounts=mounts((reference.parent, True), (Path(alignment['bam']).parent, True),
                                  (Path(callset['intervals']).parent, True), (directory, False)),
                    gpu=True, cpus=cpus, memory_gib=memory)
                runtime.run(['samtools', 'quickcheck', '-v', alignment['bam']], 'source-bam-check')
                runtime.run(command, 'haplotypecaller')
                checked_vcf(runtime, cfg, target, [sid], 'gvcf')
                return {'outputs': ['variants.g.vcf.gz', 'variants.g.vcf.gz.tbi'], 'metadata': parameters}
            directory = Path(runtime.checkpoint('gvcf/' + callset['name'] + '/' + sid,
                [reference, Path(str(reference) + '.fai'), reference.with_suffix('.dict'),
                 alignment['bam'], alignment['bai'], alignment['receipt'], callset['intervals']], parameters, produce))
            rows.append({'sample_id': sid, 'callset': callset['name'], 'ploidy': str(ploidy),
                         'gvcf': str(directory / 'variants.g.vcf.gz'), 'intervals': str(Path(callset['intervals']).resolve()),
                         'intervals_sha256': callset['intervals_sha256'], 'reference_sha256': cfg['reference']['sha256']})
    manifest = workdir / 'gvcfs.tsv'
    with manifest.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    (workdir / 'scope.json').write_text(json.dumps(scope, indent=2) + '\n')
    return {'outputs': ['gvcfs.tsv', 'scope.json'], 'metadata': {'scope': scope, 'gvcfs': len(rows), 'source_alignments_reused_without_receipt_changes': True}}


def joint_stage(cfg, workdir, runtime):
    reference, samples, callsets, scope = context(cfg)
    ids = [s['sample_id'] for s in samples]
    call_root = Path(runtime.stage_path('call'))
    rows = read_tsv(call_root / 'gvcfs.tsv')
    tool = cfg['tools'].get('bcftools', 'bcftools')
    resources = cfg.get('resources', {})
    cpus = int(resources.get('joint_cpus', 8))
    memory = int(resources.get('joint_memory_gib', 48))
    heap = int(resources.get('joint_java_heap_gib', 32))
    require(1 <= cpus <= 30 and 4 <= heap < memory <= 180, 'Invalid joint CPU/memory limits.')
    cohort_parts = []
    for callset in callsets:
        selected = callset.get('samples', ids)
        candidates = [r for r in rows if r['callset'] == callset['name']]
        require(len(candidates) == len(selected) and {r['sample_id'] for r in candidates} == set(selected), 'Incomplete gVCF partition manifest.')
        partitions = []
        for ploidy in sorted({r['ploidy'] for r in candidates}):
            group = [next(r for r in candidates if r['sample_id'] == sid) for sid in ids
                     if any(r['sample_id'] == sid and r['ploidy'] == ploidy for r in candidates)]
            for row in group:
                require(row['reference_sha256'] == cfg['reference']['sha256']
                        and row['intervals_sha256'] == callset['intervals_sha256'], 'Incompatible gVCF reference/interval partition.')
                expected_ploidy = str(callset.get('per_sample_ploidy', {}).get(row['sample_id'], callset.get('ploidy', 2)))
                require(row['ploidy'] == expected_ploidy, 'gVCF ploidy differs from the current explicit policy.')
                actual = runtime.run([tool, 'query', '-l', row['gvcf']], 'source-gvcf-sample', capture=True).splitlines()
                if not runtime.dry_run:
                    require(actual == [row['sample_id']], 'Source gVCF sample identity mismatch; no sample-map renaming is allowed.')
            parameters = {'callset': callset, 'ploidy': ploidy, 'samples': [r['sample_id'] for r in group],
                          'gatk_image': cfg['tools']['gatk_image'], 'gatk_version': GATK_VERSION,
                          'variant_call_confidence': 0, 'missing_to_reference': False}
            def produce(directory, group=group, ploidy=ploidy, callset=callset):
                directory = Path(directory)
                if not runtime.dry_run:
                    (directory / 'tmp').mkdir()
                sample_map = directory / 'sample-map.tsv'
                sample_map.write_text(''.join(r['sample_id'] + '\t' + r['gvcf'] + '\n' for r in group))
                db = directory / 'genomicsdb'
                target = directory / 'partition.raw.vcf.gz'
                common_mounts = mounts((reference.parent, True), (Path(callset['intervals']).parent, True),
                                       *((Path(r['gvcf']).parent, True) for r in group), (directory, False))
                version = runtime.run(runtime.docker(cfg['tools']['gatk_image'], ['gatk', '--version'],
                    mounts=common_mounts, gpu=False, cpus=2, memory_gib=4), 'gatk-version', capture=True)
                if not runtime.dry_run:
                    require(re.search(r'(?<![0-9])4\.3\.0\.0(?![0-9])', version), 'Actual GATK version is not 4.3.0.0.')
                commands = [
                    ['gatk', '--java-options', '-Xmx32g -XX:ParallelGCThreads=4', 'GenomicsDBImport',
                     '--genomicsdb-workspace-path', str(db), '--sample-name-map', str(sample_map),
                     '--reader-threads', '4', '--batch-size', '50', '-L', str(callset['intervals']),
                     '--tmp-dir', str(directory / 'tmp')],
                    ['gatk', '--java-options', '-Xmx32g -XX:ParallelGCThreads=4', 'GenotypeGVCFs',
                     '-R', str(reference), '-V', 'gendb://' + str(db), '-O', str(target),
                     '--sample-ploidy', str(ploidy), '--standard-min-confidence-threshold-for-calling', '0',
                     '-L', str(callset['intervals']), '--tmp-dir', str(directory / 'tmp')],
                ]
                for index, command in enumerate(commands):
                    command[2] = f'-Xmx{heap}g -XX:ParallelGCThreads={min(4, cpus)}'
                    runtime.run(runtime.docker(cfg['tools']['gatk_image'], command, mounts=common_mounts,
                                gpu=False, cpus=cpus, memory_gib=memory), 'gatk-' + str(index))
                checked_vcf(runtime, cfg, target, [r['sample_id'] for r in group], 'joint-partition')
                return {'outputs': ['partition.raw.vcf.gz', 'partition.raw.vcf.gz.tbi', 'sample-map.tsv'],
                        'metadata': parameters}
            directory = Path(runtime.checkpoint('joint/' + callset['name'] + '/p' + ploidy,
                [reference, callset['intervals'], *[p for r in group for p in (r['gvcf'], r['gvcf'] + '.tbi')]], parameters, produce))
            partitions.append(directory / 'partition.raw.vcf.gz')
        target = workdir / ('callset-' + callset['name'] + '.raw.vcf.gz')
        if len(partitions) > 1:
            runtime.run([tool, 'merge', '-m', 'all', '-Oz', '-o', str(target), *map(str, partitions)], 'merge-ploidy-partitions')
        else:
            runtime.run([tool, 'view', '-Oz', '-o', str(target), str(partitions[0])], 'copy-joint-partition')
        # Reorder samples explicitly after a partition merge; never fill absent GT with reference.
        ordered = workdir / ('callset-' + callset['name'] + '.ordered.vcf.gz')
        runtime.run([tool, 'view', '-s', ','.join(s for s in ids if s in selected), '-Oz', '-o', str(ordered), str(target)], 'order-cohort-samples')
        checked_vcf(runtime, cfg, ordered, [s for s in ids if s in selected], 'joint-callset')
        if set(selected) != set(ids):
            expanded_plain = workdir / ('callset-' + callset['name'] + '.expanded.vcf')
            expanded = Path(str(expanded_plain) + '.gz')
            if not runtime.dry_run:
                expand_sample_columns(ordered, expanded_plain, ids)
            runtime.run([tool, 'view', '-Oz', '-o', str(expanded), str(expanded_plain)], 'expand-scope-samples')
            checked_vcf(runtime, cfg, expanded, ids, 'scope-samples')
            ordered = expanded
        cohort_parts.append((ordered, ids))
    cohort = workdir / 'cohort.raw.vcf.gz'
    if len(cohort_parts) == 1:
        runtime.run([tool, 'view', '-Oz', '-o', str(cohort), str(cohort_parts[0][0])], 'cohort-raw')
    else:
        combined = workdir / 'cohort.concat.vcf.gz'
        runtime.run([tool, 'concat', '-a', '-Oz', '-o', str(combined), *[str(p[0]) for p in cohort_parts]], 'concat-disjoint-callsets')
        runtime.run([tool, 'sort', '-Oz', '-o', str(cohort), str(combined)], 'sort-cohort')
    checked_vcf(runtime, cfg, cohort, ids, 'cohort-raw')
    (workdir / 'scope.json').write_text(json.dumps(scope, indent=2) + '\n')
    return {'outputs': ['cohort.raw.vcf.gz', 'cohort.raw.vcf.gz.tbi', 'scope.json'],
            'metadata': {'scope': scope, 'gatk_version': GATK_VERSION, 'standard_call_confidence': 0,
                         'mixed_ploidy_joint_genotyping': 'separate compatible partitions, then merge without missing-to-ref'}}


def filter_stage(cfg, workdir, runtime):
    reference, samples, callsets, scope = context(cfg)
    source = Path(runtime.stage_path('joint')) / 'cohort.raw.vcf.gz'
    tool = cfg['tools'].get('bcftools', 'bcftools')
    normalized = workdir / 'cohort.normalized.vcf.gz'
    runtime.run([tool, 'norm', '-f', str(reference), '--check-ref', 'e', '-m', '-any', '--old-rec-tag', 'ORIG',
                 '-Oz', '-o', str(normalized), str(source)], 'normalize')
    checked_vcf(runtime, cfg, normalized, [s['sample_id'] for s in samples], 'normalized')
    tagged, selected = workdir / 'cohort.filtered.vcf', workdir / 'cohort.quality.vcf'
    if not runtime.dry_run:
        summary = soft_filter_vcf(normalized, tagged, selected, cfg.get('filters'))
    else:
        summary = {'policy': filter_policy(cfg.get('filters')), 'status': 'dry-run; no VCF records processed'}
    for plain, name in [(tagged, 'cohort.filtered.vcf.gz'), (selected, 'cohort.quality.vcf.gz')]:
        target = workdir / name
        runtime.run([tool, 'view', '-Oz', '-o', str(target), str(plain)], 'compress-' + plain.stem)
        checked_vcf(runtime, cfg, target, [s['sample_id'] for s in samples], plain.stem)
    (workdir / 'filtering-policy.json').write_text(json.dumps(summary, indent=2) + '\n')
    (workdir / 'scope.json').write_text(json.dumps(scope, indent=2) + '\n')
    return {'outputs': ['cohort.normalized.vcf.gz', 'cohort.normalized.vcf.gz.tbi',
                        'cohort.filtered.vcf.gz', 'cohort.filtered.vcf.gz.tbi',
                        'cohort.quality.vcf.gz', 'cohort.quality.vcf.gz.tbi', 'filtering-policy.json', 'scope.json'],
            'metadata': {'scope': scope, 'filtering': summary, 'raw_joint_vcf_preserved': str(source)}}


def run_stage(stage, cfg, workdir, runtime):
    """Entry point; root runtime owns all execution, deadlines, locks and receipts."""
    require(stage in ('call', 'joint', 'filter'), 'Unsupported calling stage.')
    require(re.fullmatch(r'nvcr\.io/nvidia/clara/clara-parabricks@sha256:[0-9a-f]{64}', cfg['tools']['parabricks_image']),
            'Parabricks requires an immutable image digest.')
    require(re.fullmatch(r'(?:docker\.io/)?broadinstitute/gatk@sha256:[0-9a-f]{64}', cfg['tools']['gatk_image']),
            'GATK requires an immutable official broadinstitute/gatk image digest for version 4.3.0.0.')
    workdir = Path(workdir)
    require(workdir.is_dir(), 'Runtime must supply a fresh private stage directory.')
    return {'call': call_stage, 'joint': joint_stage, 'filter': filter_stage}[stage](cfg, workdir, runtime)
