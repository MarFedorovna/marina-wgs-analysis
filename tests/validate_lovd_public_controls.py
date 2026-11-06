"""Root-authorised small public-only API controls; never patient variants."""
import argparse
import datetime as dt
import importlib.util
import json
from pathlib import Path
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--collector', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--fetch-approved-mybpc3-controls', action='store_true')
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('lovd_control', args.collector)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    args.output.mkdir(parents=True, exist_ok=False)
    original = Path('/srv/marina/resources/lovd-public-basic-control-v1/IVD-0000000563.verbatim.json')
    if mod.sha(original) != '9cd638e1a859f3033bbb80242eafa450a2943d3160c457218e420a61b923fe99':
        raise ValueError('Exact previously fetched public manual control required')
    parsed = mod.rows(original.read_bytes(), 'variants', 'IVD')
    report = dict(status='PASS_ACTUAL_PUBLIC_RESPONSE_SCHEMA', patient_data_used=False,
                  collector=mod.pin(args.collector), manual_control=mod.pin(original),
                  rows=len(parsed), mapping_status=mod.mapping_status(parsed[0], None),
                  additional_controls=[], source_application_status='NOT_APPLIED')
    if args.fetch_approved_mybpc3_controls:
        for kind in ['genes', 'variants']:
            body, receipt = mod.fetch(mod.url(kind, 'MYBPC3'), 8 * 1024**2,
                                      min(mod.END, mod.now() + dt.timedelta(seconds=45)))
            if body is not None:
                target = args.output / ('MYBPC3.' + kind + '.verbatim.json')
                target.write_bytes(body)
                value = json.loads(body)
                receipt['raw'] = mod.pin(target)
                receipt['json_top_level_type'] = type(value).__name__
                try:
                    records = mod.rows(body, kind, 'MYBPC3')
                    receipt.update(schema_status='PASS', record_count=len(records),
                                   field_names=sorted({k for r in records for k in r}))
                    if kind == 'genes':
                        receipt['source_build_declaration'] = [r.get('refseq_build') for r in records]
                    else:
                        receipt['explicit_hg38_description_count'] = sum(
                            isinstance(r.get('variant_genomic'), dict) and
                            'hg38' in r['variant_genomic'] for r in records)
                except ValueError:
                    receipt['schema_status'] = 'REJECTED_NOT_APPLIED'
                receipt['coverage'] = 'BASIC_API_RESPONSE_ONLY_COMPLETENESS_NOT_ESTABLISHED'
            report['additional_controls'].append(receipt)
            time.sleep(1)
    mod.dump(args.output / 'public-control-receipt.json', report)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
