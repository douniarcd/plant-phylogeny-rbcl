"""Summarize codon-position variation and model selection from actual outputs."""
import json
import re
from collections import Counter
from pathlib import Path
from Bio import AlignIO, Phylo

def compute_diagnostics(root: Path) -> dict:
    results = root / 'results'
    alignment = AlignIO.read(results / 'rbcl.aligned.fasta', 'fasta')
    counts = Counter({1: 0, 2: 0, 3: 0})
    for i in range(alignment.get_alignment_length()):
        col = alignment[:, i].upper()
        if set(col) <= set('ACGT') and len(set(col)) > 1:
            counts[i % 3 + 1] += 1
    models = []
    in_table = False
    for line in (results / 'rbcl_ml.iqtree').read_text().splitlines():
        if line.startswith('Model') and 'w-BIC' in line and 'LogL' in line:
            in_table = True
            continue
        if not in_table:
            continue
        if not line.strip():
            if models:
                break
            continue
        fields = [field for field in line.split() if field not in ('+', '-')]
        try:
            if len(fields) not in (8, 9):
                raise ValueError('Unexpected model table row')
            models.append({'model': fields[0], 'bic': float(fields[-2]), 'weight': float(fields[-1])})
        except ValueError:
            if models:
                break
    if len(models) < 2:
        raise ValueError('ModelFinder candidate table was not recognized.')
    models.sort(key=lambda item: item['bic'])
    tree = Phylo.read(results / 'rbcl_ml.treefile', 'newick')
    summary = {'gap_free_variable_sites': sum(counts.values()),
               'gap_free_variable_by_alignment_codon_position': {str(k): v for k, v in counts.items()},
               'top_bic_models': models[:2],
               'delta_bic': round(models[1]['bic'] - models[0]['bic'], 6),
               'terminal_branch_lengths': {leaf.name: leaf.branch_length for leaf in tree.get_terminals()}}
    (results / 'site_diagnostics.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary

if __name__ == '__main__':
    compute_diagnostics(Path(__file__).resolve().parents[1])
