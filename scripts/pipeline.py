#!/usr/bin/env python3
"""Reproducible five-taxon rbcL gene-tree analysis; no simulated inputs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import Bio
from Bio import AlignIO, Entrez, Phylo, SeqIO
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw'
PROCESSED = ROOT / 'data/processed'
RESULTS = ROOT / 'results'
FIGURES = ROOT / 'figures'
LOGS = ROOT / 'logs'
METADATA = ROOT / 'metadata'
FASTA = PROCESSED / 'rbcl.fasta'
ALIGNMENT = RESULTS / 'rbcl.aligned.fasta'
PREFIX = RESULTS / 'rbcl_ml'
TREE = Path(f'{PREFIX}.treefile')
LOG = logging.getLogger('rbcl')


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def taxa() -> list[dict[str, str]]:
    records = read_tsv(METADATA / 'accessions.tsv')
    if len(records) != 5 or len({r['accession'] for r in records}) != 5:
        raise ValueError('Exactly five unique accession versions are required.')
    return records


def executable(*names: str) -> str:
    for name in names:
        if path := shutil.which(name):
            return path
    raise RuntimeError(f'Required executable not found: {" or ".join(names)}')


def preflight() -> None:
    LOG.info('Python %s; Biopython %s', sys.version.split()[0], Bio.__version__)
    missing = []
    versions = {'python': sys.version.split()[0], 'biopython': Bio.__version__}
    for names, flag in [(('seqkit',), 'version'), (('mafft',), '--version'), (('iqtree2', 'iqtree'), '--version')]:
        try:
            path = executable(*names)
            proc = subprocess.run([path, flag], capture_output=True, text=True, check=True)
            versions[names[0]] = (proc.stdout + proc.stderr).strip()
            if names[0] == 'iqtree2' and not re.search(r'version\s+2\.', versions[names[0]], re.I):
                raise RuntimeError('IQ-TREE major version 2 is required.')
            LOG.info('%s: %s', names[0], versions[names[0]])
        except (RuntimeError, subprocess.CalledProcessError) as exc:
            missing.append(str(exc))
            LOG.error('%s', exc)
    import matplotlib
    versions.update(matplotlib=matplotlib.__version__)
    (LOGS / 'software_versions.json').write_text(json.dumps(versions, indent=2) + '\n')
    if missing:
        raise RuntimeError('Preflight failed: ' + '; '.join(missing))


def run_command(args: list[str], name: str, stdout_path: Path | None = None) -> None:
    LOG.info('Executing %s', subprocess.list2cmdline(args))
    with (LOGS / 'commands.jsonl').open('a', encoding='utf-8') as handle:
        handle.write(json.dumps({'utc': utc(), 'argv': args, 'cwd': str(ROOT)}) + '\n')
    with (LOGS / f'{name}.stderr.log').open('w', encoding='utf-8') as err:
        with (stdout_path or LOGS / f'{name}.stdout.log').open('w', encoding='utf-8') as out:
            subprocess.run(args, cwd=ROOT, stdout=out, stderr=err, check=True)


def validate_record(record: SeqRecord, item: dict[str, str]) -> None:
    if record.id != item['accession']:
        raise ValueError(f"Expected {item['accession']}, received {record.id}")
    if not record.annotations.get('organism', '').startswith(item['species']):
        raise ValueError(f'Unexpected organism in {record.id}: {record.annotations.get("organism")}')
    if len(record.seq) < 100_000 or not record.features:
        raise ValueError(f'{record.id}: expected complete annotated chloroplast genome.')


def fetch_records(email: str | None, offline: bool) -> None:
    if not offline:
        if not email or '@' not in email:
            raise ValueError('Provide a real contact email via --email or NCBI_EMAIL, per NCBI policy.')
        Entrez.email = email
        Entrez.tool = 'cbb26_rbcl_course_project'
        Entrez.max_tries = 2
        Entrez.sleep_between_tries = 2
        if os.getenv('NCBI_API_KEY'):
            Entrez.api_key = os.environ['NCBI_API_KEY']
    manifest = []
    for item in taxa():
        accession = item['accession']
        gb = RAW / f'{accession}.gb'
        source = 'existing local GenBank file'
        if not gb.exists():
            if offline:
                raise FileNotFoundError(f'Offline input missing: {gb}')
            LOG.info('Fetching %s from NCBI EFetch', accession)
            with Entrez.efetch(db='nuccore', id=accession, rettype='gbwithparts', retmode='text') as response:
                text = response.read()
            record = SeqIO.read(io.StringIO(text), 'genbank')
            validate_record(record, item)
            if not text.rstrip().endswith('//'):
                raise ValueError(f'{accession}: truncated GenBank response')
            with tempfile.NamedTemporaryFile(mode='w', dir=RAW, delete=False, encoding='utf-8') as handle:
                handle.write(text)
                temp = Path(handle.name)
            temp.replace(gb)
            source = 'NCBI EFetch (gbwithparts)'
        record = SeqIO.read(gb, 'genbank')
        validate_record(record, item)
        genome = RAW / f'{accession}.fasta'
        SeqIO.write(record, genome, 'fasta')
        manifest.append({**item, 'record_organism': record.annotations['organism'], 'genome_bp': len(record.seq),
                         'source': source, 'record_date': record.annotations.get('date', ''), 'verified_utc': utc(),
                         'genbank_sha256': sha256(gb), 'fasta_sha256': sha256(genome),
                         'ncbi_url': f'https://www.ncbi.nlm.nih.gov/nuccore/{accession}'})
    write_tsv(METADATA / 'record_manifest.tsv', manifest, list(manifest[0]))


def extract_rbcl() -> None:
    sequences, proteins, rows = [], [], []
    for item in taxa():
        record = SeqIO.read(RAW / f"{item['accession']}.gb", 'genbank')
        validate_record(record, item)
        features = [f for f in record.features if f.type == 'CDS' and
                    any(g.casefold() == 'rbcl' for g in f.qualifiers.get('gene', []))]
        if len(features) != 1:
            raise ValueError(f'{record.id}: expected exactly one annotated rbcL CDS; found {len(features)}')
        feature = features[0]
        q = feature.qualifiers
        seq = feature.extract(record.seq)
        codon_start = int(q.get('codon_start', ['1'])[0])
        table = int(q.get('transl_table', ['11'])[0])
        coding = seq[codon_start - 1:]
        if len(coding) % 3 or set(str(coding).upper()) - set('ACGT'):
            raise ValueError(f'{record.id}: CDS has incomplete codons or ambiguous bases; review manually.')
        translated = coding.translate(table=table)
        internal_stops = str(translated).rstrip('*').count('*')
        if internal_stops or not str(translated).endswith('*'):
            raise ValueError(f'{record.id}: unexpected stop codon pattern; review annotation.')
        # cds=True also validates the start codon and translates permitted alternative starts as M.
        protein = coding.translate(table=table, cds=True)
        annotation = ''.join(q.get('translation', [])).replace(' ', '').rstrip('*')
        if not annotation or str(protein) != annotation:
            raise ValueError(f'{record.id}: computed translation does not match annotation; review manually.')
        label = item['label']
        sequences.append(SeqRecord(seq, id=label, description=f'{record.id} rbcL CDS'))
        proteins.append(SeqRecord(protein, id=label, description=f'{record.id} rbcL protein'))
        rows.append({**item, 'location_biopython_0_based': str(feature.location),
                     'start_1_based': int(feature.location.start) + 1, 'end_1_based_inclusive': int(feature.location.end),
                     'strand': feature.location.strand, 'cds_bp': len(seq), 'protein_aa': len(protein),
                     'codon_start': codon_start, 'translation_table': table, 'internal_stops': internal_stops,
                     'terminal_stop': True, 'ambiguous_bases': 0, 'annotation_translation_matches': True,
                     'protein_id': q.get('protein_id', [''])[0]})
    SeqIO.write(sequences, FASTA, 'fasta')
    SeqIO.write(proteins, PROCESSED / 'rbcl.proteins.fasta', 'fasta')
    write_tsv(RESULTS / 'sequence_qc.tsv', rows, list(rows[0]))
    LOG.info('Extracted and validated %d real rbcL CDS records', len(sequences))


def align_sequences(threads: int) -> None:
    run_command([executable('seqkit'), 'stats', '-a', '-T', str(FASTA)], 'seqkit', RESULTS / 'seqkit_stats.tsv')
    # --auto selects the MAFFT algorithm; no trimming or column removal.
    temp = RESULTS / 'rbcl.aligned.tmp.fasta'
    run_command([executable('mafft'), '--auto', '--thread', str(threads), str(FASTA)], 'mafft', temp)
    aln = AlignIO.read(temp, 'fasta')
    originals = {r.id: str(r.seq).upper() for r in SeqIO.parse(FASTA, 'fasta')}
    if len(aln) != 5 or {r.id for r in aln} != set(originals):
        raise ValueError('Alignment has unexpected or duplicate sequence identifiers.')
    for record in aln:
        if str(record.seq).upper().replace('-', '') != originals[record.id]:
            raise ValueError(f'{record.id}: ungapped alignment differs from the original CDS.')
    temp.replace(ALIGNMENT)
    n = aln.get_alignment_length()
    variable, informative, gap_columns = 0, 0, 0
    for i in range(n):
        col = aln[:, i].upper()
        counts = Counter(base for base in col if base in 'ACGT')
        variable += len(counts) > 1
        informative += sum(count >= 2 for count in counts.values()) >= 2
        gap_columns += '-' in col
    stats = {'taxa': len(aln), 'alignment_bp': n, 'variable_sites': variable,
             'parsimony_informative_sites': informative, 'gap_columns': gap_columns,
             'trimming': 'none', 'alignment_sha256': sha256(ALIGNMENT)}
    (RESULTS / 'alignment_qc.json').write_text(json.dumps(stats, indent=2) + '\n')
    rows = [{'label': r.id, 'alignment_bp': n, 'gap_count': str(r.seq).count('-'),
             'gap_pct': round(100 * str(r.seq).count('-') / n, 3)} for r in aln]
    write_tsv(RESULTS / 'alignment_per_taxon.tsv', rows, list(rows[0]))
    pairs = []
    for i in range(len(aln)):
        for j in range(i + 1, len(aln)):
            valid = [(a, b) for a, b in zip(str(aln[i].seq).upper(), str(aln[j].seq).upper()) if a in 'ACGT' and b in 'ACGT']
            if not valid:
                raise ValueError('No comparable bases in pairwise alignment.')
            pairs.append({'taxon_a': aln[i].id, 'taxon_b': aln[j].id, 'compared_sites': len(valid),
                          'identity_pct': round(100 * sum(a == b for a, b in valid) / len(valid), 4)})
    write_tsv(RESULTS / 'pairwise_identity.tsv', pairs, list(pairs[0]))
    LOG.info('Alignment QC: %s', stats)


def infer_tree(threads: int, seed: int) -> None:
    if not ALIGNMENT.exists():
        raise FileNotFoundError('Run alignment first.')
    # Refuse to overwrite previous runs or inadvertently resume a different analysis.
    if list(RESULTS.glob('rbcl_ml.*')) or (RESULTS / 'inference_provenance.json').exists():
        raise RuntimeError('Existing IQ-TREE outputs detected. Archive them before a new inference.')
    run_command([executable('iqtree2', 'iqtree'), '-s', str(ALIGNMENT), '-st', 'DNA',
                 '-m', 'MFP', '-B', '1000', '-T', str(threads), '-seed', str(seed),
                 '-pre', str(PREFIX)], 'iqtree')
    if not TREE.exists() or not Path(f'{PREFIX}.iqtree').exists():
        raise RuntimeError('IQ-TREE did not produce the required tree and model report.')
    tree = Phylo.read(TREE, 'newick')
    if len(tree.get_terminals()) != 5 or {x.name for x in tree.get_terminals()} != {r['label'] for r in taxa()}:
        raise ValueError('IQ-TREE output has unexpected taxa.')
    (RESULTS / 'inference_provenance.json').write_text(json.dumps({
        'alignment_sha256': sha256(ALIGNMENT), 'tree_sha256': sha256(TREE),
        'iqtree_report_sha256': sha256(Path(f'{PREFIX}.iqtree')),
        'seed': seed, 'threads': threads, 'bootstrap_replicates': 1000,
        'completed_utc': utc()}, indent=2) + '\n')
    LOG.info('IQ-TREE inference completed.')


def tree_evidence() -> dict:
    provenance = json.loads((RESULTS / 'inference_provenance.json').read_text())
    alignment_qc = json.loads((RESULTS / 'alignment_qc.json').read_text())
    if (provenance['alignment_sha256'] != sha256(ALIGNMENT)
            or alignment_qc['alignment_sha256'] != sha256(ALIGNMENT)
            or provenance['tree_sha256'] != sha256(TREE)
            or provenance['iqtree_report_sha256'] != sha256(Path(f'{PREFIX}.iqtree'))):
        raise ValueError('Evidence changed after inference; summary refused.')
    tree = Phylo.read(TREE, 'newick')
    all_labels = {r['label'] for r in taxa()}
    if {x.name for x in tree.get_terminals()} != all_labels:
        raise ValueError('Unexpected terminal taxa.')
    expected = {r['label'] for r in taxa() if r['group'] == 'Poaceae'}
    expected_key = frozenset(all_labels - expected)
    splits = {}
    for clade in tree.find_clades(order='preorder'):
        if clade is tree.root or clade.is_terminal():
            continue
        side = {x.name for x in clade.get_terminals()}
        other = all_labels - side
        if min(len(side), len(other)) < 2:
            continue
        key = frozenset(side if len(side) < len(other) else other)
        support = clade.confidence
        if support is None and clade.name:
            try:
                support = float(clade.name)
            except ValueError:
                pass
        if support is None:
            raise ValueError('Internal branch lacks ultrafast bootstrap support.')
        if key in splits and splits[key] != float(support):
            raise ValueError('Conflicting support for the same unrooted split.')
        splits[key] = float(support)
    if len(splits) != 2:
        raise ValueError(f'Expected two nontrivial splits in a resolved five-taxon tree; found {len(splits)}')
    rows = [{'side_a': '; '.join(sorted(k)), 'side_b': '; '.join(sorted(all_labels - k)),
             'ufbootstrap_pct': v, 'predicted_split': k == expected_key} for k, v in splits.items()]
    write_tsv(RESULTS / 'branch_support.tsv', rows, list(rows[0]))
    report = Path(f'{PREFIX}.iqtree').read_text(encoding='utf-8')
    match = re.search(r'Best-fit model according to BIC:\s*([^\r\n]+)', report)
    if not match:
        raise ValueError('Could not find the BIC-selected model in IQ-TREE report.')
    evidence = {'model': match.group(1).strip(), 'predicted_split_present': expected_key in splits,
                'predicted_split_support': splits.get(expected_key), 'bootstrap_replicates_requested': 1000,
                'tree_sha256': sha256(TREE), 'alignment_sha256': sha256(ALIGNMENT)}
    (RESULTS / 'tree_evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    return evidence


def plot_unrooted() -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    tree = Phylo.read(TREE, 'newick')
    leaves = tree.get_terminals()
    # Layout only: preserve every inferred branch length and the unrooted topology.
    preferred = {'Arabidopsis_thaliana': 55, 'Brassica_rapa': -55,
                 'Oryza_sativa': 105, 'Zea_mays': 160, 'Triticum_aestivum': 245}
    angles = {leaf: math.radians(preferred[leaf.name]) if leaf.name in preferred else
              2 * math.pi * i / len(leaves) for i, leaf in enumerate(leaves)}
    positions = {tree.root: (0., 0.)}
    fig, ax = plt.subplots(figsize=(10, 4.6))
    for parent in tree.find_clades(order='preorder'):
        x, y = positions[parent]
        for child in parent.clades:
            angle = sum(angles[leaf] for leaf in child.get_terminals()) / len(child.get_terminals())
            length = child.branch_length
            if length is None or length < 0:
                raise ValueError('Missing or invalid inferred branch length.')
            xx, yy = x + length * math.cos(angle), y + length * math.sin(angle)
            positions[child] = (xx, yy)
            ax.plot([x, xx], [y, yy], color='#24495a', lw=1.6)
            if not child.is_terminal() and child.confidence is not None:
                offset = (6, -19) if length < .04 else (0, 9)
                ax.annotate(f'{child.confidence:g}', ((x + xx) / 2, (y + yy) / 2), xytext=offset,
                            textcoords='offset points', ha='center', fontsize=12, color='#194b39',
                            bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 0.8})
    scale = max(max(abs(x), abs(y)) for x, y in positions.values())
    if scale <= 0:
        raise ValueError('All branch lengths are zero; cannot construct scale bar.')
    for leaf in leaves:
        x, y = positions[leaf]
        left = math.cos(angles[leaf]) < 0
        dy = 9 if leaf.name == 'Oryza_sativa' else 0
        ax.annotate(leaf.name.replace('_', ' '), (x, y), xytext=(-5 if left else 5, dy),
                    textcoords='offset points', ha='right' if left else 'left', va='center',
                    fontsize=12, fontstyle='italic')
    step = 10 ** math.floor(math.log10(scale / 3))
    sx = min(x for x, _ in positions.values()) + .035
    sy = min(y for _, y in positions.values()) - .018
    ax.plot([sx, sx + step], [sy, sy], color='black', lw=2)
    ax.annotate(f'{step:g} substitutions/site', (sx + step/2, sy), xytext=(0,-8),
                textcoords='offset points', ha='center', va='top', fontsize=10)
    ax.set_aspect('equal'); ax.margins(x=.27, y=.30); ax.axis('off')
    ax.set_title('rbcL maximum-likelihood gene tree', loc='left', fontsize=14, fontweight='bold')
    for extension in ('png', 'pdf'):
        fig.savefig(FIGURES / f'rbcl_tree.{extension}', dpi=300, bbox_inches='tight', pad_inches=.15)
    plt.close(fig)


def summarize() -> None:
    tree_evidence()
    plot_unrooted()
    from diagnostics import compute_diagnostics
    compute_diagnostics(ROOT)
    LOG.info('Tree evidence, site diagnostics and figure generated.')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['preflight', 'download', 'extract', 'align', 'infer', 'summarize', 'all'])
    parser.add_argument('--email', default=os.getenv('NCBI_EMAIL'))
    parser.add_argument('--offline', action='store_true', help='Require five local GenBank files; do not fetch.')
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    if args.threads < 1 or args.seed < 1:
        parser.error('Threads and seed must be positive.')
    for path in (RAW, PROCESSED, RESULTS, FIGURES, LOGS, METADATA):
        path.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',
                        handlers=[logging.StreamHandler(), logging.FileHandler(LOGS / 'pipeline.log', encoding='utf-8')])
    os.chdir(ROOT)
    stages = {'preflight': preflight, 'download': lambda: fetch_records(args.email, args.offline),
              'extract': extract_rbcl, 'align': lambda: align_sequences(args.threads),
              'infer': lambda: infer_tree(args.threads, args.seed), 'summarize': summarize}
    try:
        for stage in stages if args.stage == 'all' else [args.stage]:
            LOG.info('Starting stage: %s', stage)
            stages[stage]()
            LOG.info('Completed stage: %s', stage)
    except Exception:
        LOG.exception('Pipeline stopped; no subsequent stages were executed.')
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
