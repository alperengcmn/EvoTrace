# CLI

- `evotrace validate --input sequences.fasta [--alphabet auto|dna|protein]`
- `evotrace analyze --input sequences.fasta --output results/ [--config config.json] [--threads N] [--dry-run] [--resume] [--verbose]`
- `evotrace pipeline` is an alias for `analyze`.
- `evotrace doctor` lists required and optional tools with versions, paths, and install guidance.
- `evotrace phylogeny --alignment aligned.fasta` runs IQ-TREE when installed.
- `evotrace selection --alignment codons.fasta --engine hyphy --method FEL` or `--engine paml --tree tree.nwk` runs a selected method.
- `evotrace domains --sequences proteins.fasta --database profiles.hmm` runs HMMER when installed.

Use `evotrace --help` and `evotrace analyze --help` for current flags.
