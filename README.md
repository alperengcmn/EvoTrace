# EvoTrace

EvoTrace is a modular platform for comparative sequence analysis. It validates inputs, aligns sequences with MAFFT, computes conservation and evolutionary summaries, builds a phylogeny, and exports machine-readable results, provenance, and a portable HTML report. Optional modules connect to IQ-TREE, HyPhy, PAML/codeml, HMMER, NCBI E-utilities, and PDB/AlphaFold structure mapping.

## Install

Python 3.9+ and MAFFT are required for the integrated workflow.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
# Install MAFFT with your platform package manager (for example: conda install -c bioconda mafft)
evotrace doctor
```

Optional web interface: `python -m pip install -e '.[web]'`, then `streamlit run scripts/streamlit_app.py`.

## Quick start

```bash
evotrace validate --input examples/tp53_mammals.fasta --alphabet protein
evotrace analyze --output results/tp53 --config configs/example.yaml
```

The main example is real mammalian TP53 protein data from NCBI; accession provenance and links are in [examples/README.md](examples/README.md). `homologs.fasta` is a synthetic software fixture only. For coding-sequence selection analysis, provide a codon-aware alignment or protein sequences plus matching CDS records for validated back-translation.

`--dry-run` validates inputs and reports the planned alignment step without running it. `--resume` reuses stage outputs only when inputs, settings, software version, and output hashes match. `--accession NP_000537.3 NP_035770.2` fetches and caches NCBI records with provenance. Config paths are resolved relative to the config file.

## Configuration and outputs

Nested configuration supports input paths/accessions, alignment, phylogeny, selection, domains, structures, reports, and workers. YAML is shown in `configs/example.yaml`; JSON is also supported. Run `evotrace doctor` to inspect optional executables and dependencies.

Outputs include the aligned FASTA, per-column CSV metrics, alignment statistics, pairwise identity matrix, descriptive constraint/divergence summaries, candidate site CSV/TSV/JSON, tree Newick/SVG, conservation SVG, HTML report, provenance, stage manifest, and run summary. Positions are 1-based. Entropy is over observed canonical non-gap residues; gaps and ambiguous symbols are reported separately. Normalized entropy uses log2(4) for nucleotide data and log2(20) for protein data. These alignment summaries do not establish adaptation or significance. Pairwise Nei–Gojobori dN/dS is descriptive, not a positive-selection test.

## Methods and limitations

See [methodology](docs/methodology.md), [architecture](docs/architecture.md), [reproducibility](docs/reproducibility.md), [configuration](docs/configuration.md), [pipeline](docs/pipeline.md), and [troubleshooting](docs/troubleshooting.md). EvoTrace does not perform orthology validation. Neighbor Joining uses uncorrected p-distance and has no branch support. Optional external analyses are reported as unavailable if their executable or database is missing; never interpret missing results as evidence of no effect.

An independent ProteinGym v1.3 DMS validation of the alignment-derived constraint component, including its scope and limitations, is documented in [validation](docs/validation.md).

## Development

```bash
python -m pip install -e '.[dev]'
python -m pytest
ruff check .
mypy src
python -m build
```

Docker: `docker build -t evotrace .` then `docker run --rm -v "$PWD:/work" evotrace analyze --input /work/examples/tp53_mammals.fasta --output /work/results`.

## Citation and license

MIT License. Cite the software version and the methods/tools actually run in each analysis. This project does not yet have a peer-reviewed publication or archived release DOI; see `CITATION.cff`.
