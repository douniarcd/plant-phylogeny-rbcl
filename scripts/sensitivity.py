#!/usr/bin/env python3
"""Optional exploratory site-filtered IQ-TREE reruns, separate from the primary analysis.

The alignment-column rules are explicit. Support is reported for splits actually
present in each inferred tree, not for absent splits. Original outputs are untouched.
"""
from __future__ import annotations

import argparse
import hashlib
import re
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess

from Bio import AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--output', type=Path, required=True, help='New output directory')
    args = parser.parse_args()
    if args.threads < 1 or args.seed < 1:
        parser.error('Threads and seed must be positive')
    binary = shutil.which('iqtree2') or shutil.which('iqtree')
    if not binary:
        raise RuntimeError('IQ-TREE 2 is required')
    version = subprocess.run([binary, '--version'], check=True, capture_output=True, text=True)
    if 'version 2.' not in version.stdout + version.stderr:
        raise RuntimeError('IQ-TREE major version 2 is required')
    aln = AlignIO.read(ROOT / 'results/rbcl.aligned.fasta', 'fasta')
    n = aln.get_alignment_length()
    if n != 1443 or len(aln) != 5:
        raise ValueError('Expected the original five-taxon, 1443-column alignment')
    variants = {
        'gap_free': [i for i in range(n) if '-' not in aln[:, i]],
        'tail_removed': list(range(1400)),
        'positions_1_2': [i for i in range(n) if i % 3 != 2],
        'position_3': [i for i in range(n) if i % 3 == 2],
    }
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    summary = {}
    provenance = {'alignment_sha256': hashlib.sha256((ROOT / 'results/rbcl.aligned.fasta').read_bytes()).hexdigest(), 'iqtree_version': (version.stdout + version.stderr).strip(), 'threads': args.threads, 'seed': args.seed, 'bootstrap_replicates': 1000, 'started_utc': datetime.now(timezone.utc).isoformat()}
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    for name, columns in variants.items():
        folder = out / name
        folder.mkdir()
        fasta = folder / 'alignment.fasta'
        records = [SeqRecord(Seq(''.join(str(rec.seq[i]) for i in columns)), id=rec.id, description='') for rec in aln]
        SeqIO.write(records, fasta, 'fasta')
        (folder / 'retained_columns.json').write_text(json.dumps([i+1 for i in columns]) + '\n')
        command = [binary, '-s', str(fasta), '-st', 'DNA', '-m', 'MFP', '-B', '1000', '-T', str(args.threads), '-seed', str(args.seed), '-pre', str(folder / 'tree')]
        (folder / 'command.json').write_text(json.dumps(command, indent=2) + '\n')
        with (folder / 'stdout.log').open('w') as stdout, (folder / 'stderr.log').open('w') as stderr:
            subprocess.run(command, check=True, stdout=stdout, stderr=stderr)
        tree = Phylo.read(folder / 'tree.treefile', 'newick')
        labels = {c.name for c in tree.get_terminals()}
        splits = {}
        for clade in tree.get_nonterminals():
            side = {c.name for c in clade.get_terminals()}
            other = labels - side
            if min(len(side), len(other)) >= 2:
                key = '; '.join(sorted(side if len(side) < len(other) else other))
                splits[key] = clade.confidence
        if labels != {rec.id for rec in aln} or len(splits) != 2 or any(v is None for v in splits.values()):
            raise ValueError(f'{name}: unexpected tree taxa or missing support')
        model = re.search(r'Best-fit model according to BIC:\s*([^\r\n]+)', (folder / 'tree.iqtree').read_text())
        if not model:
            raise ValueError(f'{name}: missing selected model')
        summary[name] = {'columns': len(columns), 'model': model.group(1).strip(), 'observed_splits_ufboot_pct': splits}
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')


if __name__ == '__main__':
    main()
