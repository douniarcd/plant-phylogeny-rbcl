# Visual and computational quality control

## Alignment — 26 November 2026

The MAFFT alignment contains the five expected sequences and 1,443 columns. Ungapped sequences match the extracted CDS, and all sequences begin with ATG. No ambiguous bases or all-gap columns were detected. There are 258 variable sites, 177 parsimony-informative sites and 15 gap-containing columns.

AliView inspection is documented by three screenshots:

- `figures/aliview_start.png`: sequence labels, conserved ATG start codons and nucleotide substitutions.
- `figures/aliview_internal_gap.png`: columns 1404–1406 contain a gap in all taxa except maize.
- `figures/aliview_end.png`: terminal gaps in maize (1432–1443) and rice and wheat (1438–1443).

No unexpected anomalies were visible in the displayed regions. The screenshots do not show the middle of the alignment; full-length integrity was checked computationally. All gap runs are multiples of three. The alignment was retained without trimming, including annotated stop codons. The AliView version was not visible in the screenshots.

## Tree

The Newick tree was parsed and the unrooted Python rendering visually inspected. Five species labels, nonnegative branch lengths, 100% and 80% UFBoot support values, and a scale bar in substitutions per site were verified. FigTree was not used.
