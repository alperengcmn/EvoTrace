# Methodology

## Input and alignment

FASTA IDs must be unique; records must be non-empty and consist of recognized nucleotide or protein symbols. `auto` selects nucleotide only when every symbol is compatible with the nucleotide alphabet, otherwise protein. MAFFT is invoked as a fixed executable with an argument vector (no shell expansion), and its output is checked for matching record count and equal aligned lengths.

## Alignment summaries

Column occupancy is the fraction of sequences with a canonical non-gap symbol. Gap fraction is the fraction of `-` symbols. Pairwise identity is exact matches divided by positions where both sequences have canonical symbols. A BLOSUM62 similarity summary is also provided for protein alignments; positive substitution scores count as similar. These summaries do not replace model-based evolutionary distance.

## Entropy and conservation

At each alignment column, frequencies are calculated among observed residues only. Shannon entropy is H = -Σ pᵢ log₂(pᵢ), in bits. Gaps and `?` are omitted from this distribution and separately quantified. Normalized entropy is H/log₂(K), with K=4 for nucleotide and K=20 for protein data; it can be below 1 when ambiguity/stop symbols are present. Effective alphabet size is 2ᴴ. Dominant-residue frequency is the conservation proportion. These are descriptive alignment statistics and depend on sampling and alignment quality.

## Pairwise dN/dS

The optional pairwise estimator requires a validated codon alignment. Protein alignments can be back-translated when matching CDS records are supplied and their translation exactly agrees. It averages synonymous/nonsynonymous changes over minimum mutation paths, then applies a Jukes–Cantor correction. Gaps, ambiguous codons, and stop-containing codons are excluded. Saturation and zero dS can yield null estimates. This is descriptive, not a likelihood-ratio test or evidence of positive selection. HyPhy/PAML adapters run only when selected and their external programs are installed.

## Evolutionary Trace Score (ETS)

For at least three aligned records, the pipeline additionally writes a Neighbor Joining tree based on uncorrected pairwise p-distance. The tree has no bootstrap or other branch support and is exploratory only.

The current alignment score output consists of two transparent descriptive values per position:

- Constraint = dominant-residue frequency × occupancy.
- Divergence = normalized entropy × occupancy.

The values lie in [0,1] under the supported alphabets. They are correlated views of the alignment and are not independent evidence streams or measures of adaptation.

`ets_descriptive_score` is an export alias for divergence, not a third measure. Selection tests are retained at their native site/branch/gene levels and grouped as one dependent family; domain, motif, and structure data add annotation context only. No composite adaptation metric is generated from these heterogeneous inputs.

## Motifs, domains, structures

Motifs are user-supplied regular expressions and mapped through the alignment coordinate map. HMMER searches run only with a protein input, configured profile database, and installed executable. PDB/mmCIF or AlphaFold DB models can be mapped to a selected reference sequence with explicit chain/residue coordinates and confidence metadata. Missing tools, databases, structures, or hits are represented with stage status and are not interpreted as biological negatives.
