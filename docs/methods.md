# Scientific methods

**Course:** CBB26, Session 12, Project 3. Submission instructions: pages 9–10 and 17–20.

**Question:** Do rice, maize and wheat group separately from Arabidopsis and Brassica in the chloroplast rbcL gene tree?

**Hypothesis:** The maximum-likelihood tree contains the unrooted split Poaceae | Brassicaceae. This is a prediction, not a result.

## Data

| Species | Versioned NCBI accession |
|---|---|
| Arabidopsis thaliana | [NC_000932.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_000932.1) |
| Oryza sativa | [NC_001320.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001320.1) |
| Zea mays | [NC_001666.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_001666.2) |
| Triticum aestivum | [NC_002762.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_002762.1) |
| Brassica rapa | [NC_040849.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_040849.1) |

## Workflow

1. Download annotated GenBank records using NCBI EFetch. Verify accession versions and organisms; export matching full-genome FASTA and SHA-256 checksums.
2. Extract one annotated rbcL CDS per genome using Biopython feature extraction, preserving biological orientation. Check sequence length, ambiguous bases, translation, terminal stop and agreement with the GenBank protein annotation. Record coordinates and strand.
3. Run SeqKit statistics and MAFFT `--auto` on the five nucleotide CDS. Verify taxa and ungapped sequence equality. Record alignment length, variable and parsimony-informative sites, gaps and pairwise identity. Inspect in AliView. Trimming: none; annotated stop codons were retained. AliView screenshots document the start, internal gap and terminal regions.
4. Run IQ-TREE 2 with `-m MFP -B 1000 -st DNA -seed 42 -T 2`. Save the selected BIC model, logs and Newick tree.
5. Visualize the unrooted tree with species names, branch support and substitutions/site scale. The tree was rendered and inspected using Python; FigTree was not used. Test the family-level bipartition and report both internal branches. Use 95% UFBoot as a common strong-support guideline, not as a probability of correctness.
6. Summarize the observed splits and site diagnostics; evaluate the hypothesis in the scientific report. No outgroup is specified. The analysis does not estimate divergence dates or establish a species tree.

## References and software

- NCBI E-utilities: https://www.ncbi.nlm.nih.gov/books/NBK25501/
- Cock et al. (2009). Biopython. https://doi.org/10.1093/bioinformatics/btp163
- Shen et al. (2016). SeqKit. https://doi.org/10.1371/journal.pone.0163962
- Katoh & Standley (2013). MAFFT v7. https://doi.org/10.1093/molbev/mst010
- Larsson (2014). AliView. https://doi.org/10.1093/bioinformatics/btu531
- Kalyaanamoorthy et al. (2017). ModelFinder. https://doi.org/10.1038/nmeth.4285
- Hoang et al. (2018). UFBoot2. https://doi.org/10.1093/molbev/msx281
- Minh et al. (2020). IQ-TREE 2. https://doi.org/10.1093/molbev/msaa015
- FigTree: https://github.com/rambaut/figtree/releases

- Grass Phylogeny Working Group II (2012). New grass phylogeny resolves deep evolutionary relationships and discovers C4 origins. https://doi.org/10.1111/j.1469-8137.2011.03972.x
- Nuclear grass phylogenomics: https://doi.org/10.1111/nph.20263
