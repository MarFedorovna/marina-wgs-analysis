#!/usr/bin/env python3
"""Render unique-allele count histograms using the pinned report environment."""
import argparse
import hashlib
import json
from pathlib import Path
import re


def render(counts_path, output):
    import matplotlib
    if matplotlib.__version__ != '3.10.8':
        raise ValueError('Use the pinned report Python with matplotlib==3.10.8.')
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    plt.rcParams.update({'svg.hashsalt': 'marina-count-hist-v1', 'font.family': 'DejaVu Sans', 'font.size': 10})
    counts = json.loads(counts_path.read_text())
    output.mkdir(parents=True, exist_ok=True)
    scopes = [('cohort', counts['cohort'])] + [('sample_' + name, value) for name, value in sorted(counts['samples'].items())]
    kinds = [('variant_type_histogram', 'Variant types (exclusive)'),
             ('exclusive_primary_consequence', 'Primary consequence (exclusive)'),
             ('nonexclusive_flags', 'Consequence flags (overlap; do not sum)')]
    levels = ['raw', 'qc', 'candidate']
    colors = ['#64748b', '#2563eb', '#c2410c']
    plots = []
    for scope, counters in scopes:
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,120}', scope):
            raise ValueError('Unsafe plot scope name.')
        for level in levels:
            item = counters[level]
            total = item['unique_normalized_alleles']
            if not isinstance(total, int) or total < 0:
                raise ValueError('Invalid unique-allele count.')
            for exclusive in ('variant_type_histogram', 'exclusive_primary_consequence'):
                if sum(item[exclusive].values()) != total:
                    raise ValueError('Exclusive histogram does not sum to its unique-allele denominator.')
            if any(not isinstance(n, int) or n < 0 or n > total for n in item['nonexclusive_flags'].values()):
                raise ValueError('Nonexclusive flag count exceeds the allele denominator.')
        for kind, title in kinds:
            labels = sorted(set().union(*(set(counters[level][kind]) for level in levels))) or ['no observed alleles']
            fig, ax = plt.subplots(figsize=(12, max(3.2, len(labels) * .65 + 1.5)), layout='constrained')
            for j, (level, color) in enumerate(zip(levels, colors)):
                values = [counters[level][kind].get(label, 0) for label in labels]
                bars = ax.barh([i + (j - 1) * .25 for i in range(len(labels))], values, height=.24, color=color,
                    label=level + ' (N=' + str(counters[level]['unique_normalized_alleles']) + ')')
                ax.bar_label(bars, labels=[str(n) if n else '' for n in values], padding=3, fontsize=8)
            ax.set_yticks(range(len(labels)), labels)
            ax.invert_yaxis()
            ax.set_xlabel('Unique normalized alleles; transcript rows counted once')
            ax.set_title(scope + ' — ' + title)
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.grid(axis='x', alpha=.15)
            ax.set_axisbelow(True)
            ax.legend(loc='lower right', fontsize=8)
            ax.spines[['top', 'right']].set_visible(False)
            filename = scope + '__' + kind
            for suffix in ('svg', 'png'):
                path = output / (filename + '.' + suffix)
                if path.exists():
                    raise ValueError('Refusing to overwrite a figure.')
                metadata = {'Date': None, 'Creator': 'Matplotlib 3.10.8; unique normalized allele counts'} if suffix == 'svg' else {'Software': 'Matplotlib 3.10.8'}
                fig.savefig(path, dpi=160, metadata=metadata)
            plt.close(fig)
            plots.append(dict(scope=scope, kind=kind, svg=filename + '.svg', png=filename + '.png',
                denominator={level: counters[level]['unique_normalized_alleles'] for level in levels}))
    manifest = dict(schema_version=1, matplotlib_version=matplotlib.__version__,
        count_file_sha256=hashlib.sha256(counts_path.read_bytes()).hexdigest(), plots=plots,
        note='Raw/QC/research-candidate levels are distinct; consequence flags overlap and must not be summed.')
    target = output / 'plot_manifest.json'
    with target.open('x') as stream:
        json.dump(manifest, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'plots': len(plots), 'formats': ['svg', 'png'], 'matplotlib': matplotlib.__version__}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--counts', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    render(args.counts, args.output_dir)
