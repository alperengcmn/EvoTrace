# Changelog

## Unreleased

- Added method-specific HyPhy parsers for FEL, MEME, FUBAR, aBSREL, and BUSTED, preserving posterior evidence, corrected branch p-values, and distinct site/branch/gene result levels.
- Added conservative codeml output parsing for labeled likelihood/omega values and BEB posterior sites; uncomputed model comparisons remain explicitly null.
- Added alignment QC, codon QC output, linked motif/domain/selection/structure evidence in candidate sites, execution metadata, a canonical YAML example, and a real MAFFT-backed workflow integration test.
- Expanded dependency diagnostics and the Streamlit workspace; clarified stage cache and score interpretation docs.

## 0.1.0

Initial working release: FASTA validation, MAFFT alignment, column summaries, descriptive ETS, optional descriptive pairwise dN/dS, user-defined regex motifs, provenance, SVG profile and HTML report.
