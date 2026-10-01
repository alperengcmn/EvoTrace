# Troubleshooting

- **MAFFT not found:** install MAFFT and make sure `mafft` is on `PATH`; validation and `--dry-run` do not require invoking it.
- **YAML unavailable:** use JSON configuration or install PyYAML.
- **Duplicate IDs / invalid characters:** give every FASTA record a unique first header token and remove unsupported symbols.
- **Codon analysis not run:** enable `selection` only for nucleotide alignments. For interpretation, supply a codon-aware alignment; EvoTrace skips indel/ambiguous/stop codons and does not provide a positive-selection test.
- **Resume appears stale:** `--resume` requires matching input/settings/software fingerprints and output SHA-256 hashes. If an output artifact was edited, the affected stage will be recomputed.
