#!/usr/bin/env python3
"""Network-denied reference-only HGVS proof and strictly local validation.

--public-controls never reads project variants. Optional --input-hgvs is private
local input only; references absent from the manifest are explicit unavailable.
"""
import argparse
import hashlib
import json
from pathlib import Path
import socket


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--public-controls', action='store_true')
    group.add_argument('--input-hgvs', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = json.loads((root / 'public-source-manifest.json').read_text())
    references = set(manifest['fixture_reference_ids'])
    for entry in manifest['files']:
        if entry.get('reference_id'):
            payload = (root / 'public-reference-fixtures' / entry['name']).read_bytes()
            if hashlib.sha256(payload).hexdigest() != entry['sha256']:
                raise RuntimeError('Reference digest mismatch')

    from mutalyzer.normalizer import normalize
    import mutalyzer.description
    import mutalyzer_retriever.retriever
    from mutalyzer_retriever.retriever import NoReferenceError

    network_attempts = []

    def deny(*arguments, **keywords):
        network_attempts.append('denied')
        raise RuntimeError('Offline HGVS: network access denied')

    socket.socket = deny
    socket.create_connection = deny
    socket.getaddrinfo = deny

    def local_raw(reference_id, reference_source=None, reference_type=None,
                  size_off=True, configuration_path=None, timeout=1, **unused):
        if reference_id not in references:
            raise NoReferenceError({}, [])
        extension = 'fasta' if reference_type == 'fasta' else 'gff3'
        target = root / 'public-reference-fixtures' / (reference_id + '.' + extension)
        if not target.is_file():
            raise NoReferenceError({}, [])
        return target.read_text(), extension, 'ncbi'

    # Same local fixture boundary used by the upstream test suite, pinned in manifest.
    mutalyzer_retriever.retriever.retrieve_raw = local_raw
    mutalyzer.description.get_cds_to_mrna = lambda *a, **k: None
    rows = []
    if args.public_controls:
        cases = [
            ('mt_public_snv', 'NC_012920.1:m.3243A>G', 'NC_012920.1:m.3243A>G', None),
            ('transcript_public_deletion', 'NM_003002.4:c.274del', None, None),
            ('wrong_reference_base', 'NC_012920.1:m.3243C>G', None, 'ESEQUENCEMISMATCH'),
            ('missing_reference_is_explicit', 'NM_999999999.1:c.1A>T', None, 'ERETR'),
        ]
        for pair in manifest.get('mane14_public_pairs', []):
            description = pair['RefSeq_nuc'] + ':c.1A>G'
            cases.append(('mane14_' + pair['symbol'] + '_synthetic_start_substitution', description, description, None))
    else:
        cases = [('local_input', value.strip(), None, None)
                 for value in args.input_hgvs.read_text().splitlines()
                 if value.strip() and not value.lstrip().startswith('#')]
    for name, description, expected, expected_error in cases:
        result = normalize(description)
        codes = [error['code'] for error in result.get('errors', [])]
        normalized = result.get('normalized_description')
        passed = (not codes and bool(normalized)) if expected_error is None else expected_error in codes
        if expected is not None:
            passed = passed and normalized == expected
        rows.append(dict(control=name, input=description, normalized=normalized,
                         error_codes=codes, passed=passed, result=result))
    # Explicit guard self-test; this attempt is separate from normalizer attempts.
    normalizer_network_attempts = len(network_attempts)
    try:
        socket.create_connection(('example.invalid', 443))
    except RuntimeError:
        guard_pass = True
    else:
        guard_pass = False
    passed = all(row['passed'] for row in rows) and guard_pass and normalizer_network_attempts == 0
    output = dict(schema_version=1, status='PASS' if passed else 'FAIL',
                  mode='public_controls' if args.public_controls else 'private_local_input',
                  network_denied=True, network_guard_selftest_pass=guard_pass,
                  normalizer_network_attempts=normalizer_network_attempts,
                  reference_ids=sorted(references), rows=rows,
                  project_annotation_compatibility='NOT_YET_VALIDATED: fixture controls do not certify current HCM MANE1.4/GENCODE49/VEP115 mapping',
                  clinical_classification=False)
    target = args.output or root / ('public-control-receipt.json' if args.public_controls else 'private-offline-hgvs-results.json')
    target.write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({key: value for key, value in output.items() if key != 'rows'}))
    raise SystemExit(0 if passed else 1)


if __name__ == '__main__':
    main()
