# Reproducibility

Each run records UTC time, EvoTrace/Python/platform and available tool versions, resolved input/configuration, input and normalized configuration hashes, stage statuses, outputs and hashes, command metadata, working directory, runtime, and bounded external stdout/stderr in `provenance.json` and `pipeline_stages.json`. Output files are deterministic for fixed input/software/settings except timestamp/runtime metadata and tool-specific behavior across versions. Preserve the FASTA, config, and provenance files with reported results.
