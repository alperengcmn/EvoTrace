# Local environment

## Setup

On macOS or Linux, run `bash scripts/setup_environment.sh`. The preferred setup is an isolated Conda prefix under `.evotrace/envs/evotrace`; downloaded packages are cached under `.evotrace/cache/conda-pkgs`. The script creates `.evotrace/databases`, `.evotrace/downloads`, and `.evotrace/logs` for local project resources. It will not replace an unrelated path at the selected environment location.

The Conda environment includes Python, the EvoTrace package with web and development extras, and the core executables MAFFT, IQ-TREE 2, and HMMER. HyPhy, PAML/codeml, and BLAST+ are included when the configured Conda channels provide compatible builds. With no Conda installation, the script falls back to `.venv` for Python packages; install the core command-line tools separately from trusted platform packages.

Run `bash scripts/verify_environment.sh` after activation, or set `EVOTRACE_PYTHON` to the target Python executable. It checks imports, CLI/doctor, real MAFFT/IQ-TREE/HMMER smoke workflows, and the project's test suite. `evotrace doctor --strict` fails for any missing required component; plain `evotrace doctor` reports missing optional components without failing.

## External databases

No large reference database is downloaded by environment setup. Domain annotation is optional and requires a user-selected HMM profile library; place it under `.evotrace/databases/`, run `hmmpress -f <profile.hmm>` for HMMER's pressed index files, and pass the profile path with `evotrace domains --database`. EvoTrace currently does not require a pre-downloaded BLAST database. NCBI and AlphaFold retrieval use the configured online sources and the analysis output's cache; those are not prerequisites for local FASTA analysis.

## Reproducibility and cleanup

Keep `environment.yml` with project changes to record the requested environment. The isolated environment, package cache, and local database/download/log directories are ignored by Git. To remove this project-local environment, delete `.evotrace/` after confirming that it contains no user-supplied databases or other files you want to keep. This cleanup is intentionally manual.
