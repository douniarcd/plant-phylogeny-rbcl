# CBB26 — Project 3: Plant phylogeny using chloroplast rbcL

## Objective
Test whether rice, maize and wheat form a group separate from Arabidopsis and Brassica using a five-taxon maximum-likelihood rbcL gene tree.

## Project contents

| Directory or file | Contents |
|---|---|
| `docs/scientific_report.pdf` | Scientific report |
| `docs/methods.md`, `docs/visual_review.md` | Methods, references and visual quality control |
| `data/raw/`, `data/processed/` | Versioned GenBank records, genome FASTA, extracted CDS and protein FASTA |
| `metadata/` | Accessions, record provenance and SHA-256 manifests |
| `results/` | Alignment, quality-control tables, model selection, tree and branch support |
| `figures/` | Unrooted tree and AliView evidence |
| `logs/` | Commands, software versions and execution logs |
| `scripts/` | Analysis and reproducibility workflow |
| `environment.yml` | Recorded analysis dependency specification |

## Reproduction
Use a fresh copy or the notebook, which creates a separate working directory. The archived original results must not be overwritten. The included GenBank records permit offline analysis after software installation.

```bash
conda env create -f environment.yml
conda activate rbcl-phylogeny
python scripts/pipeline.py preflight
python scripts/pipeline.py download --offline
python scripts/pipeline.py extract
python scripts/pipeline.py align --threads 2
python scripts/pipeline.py infer --threads 2 --seed 42
python scripts/pipeline.py summarize
```

For a fresh run in a copy of the project, first remove the existing `results/rbcl_ml.*` files and `results/inference_provenance.json` from that copy. The notebook performs this isolation automatically. Optional fresh acquisition uses NCBI EFetch with a contact email supplied privately through `NCBI_EMAIL`.

Additional exploratory site-filtering analyses can be executed without changing the primary alignment or tree:

```bash
python scripts/sensitivity.py --threads 1 --seed 42 --output results/sensitivity
```

The output directory must be new. The four subsets use gap-free sites, columns 1–1400, alignment codon positions 1+2, and position 3 (including gaps). Each run saves the alignment, column indices, command, model, tree, logs and bootstrap support. These supplementary runs are not included among the completed results in this archive.

## Integrity and limitations

```bash
sha256sum -c metadata/artifact_checksums.sha256
```

The primary workflow was executed using the versions in `logs/software_versions.json`; an exact recreation of the environment specification has not been tested. AliView screenshots cover the start and end regions; full-length alignment checks were computational. FigTree was not used; the unrooted figure was rendered in Python.
