# Pipeline execution

`evotrace analyze` is the canonical API used by both the CLI and Streamlit. The core dependency chain is validation → alignment → alignment QC → phylogeny and downstream per-position summaries. Motifs are mapped through the alignment; configured HMMER domains and structure mappings are projected onto the same reference/alignment positions before the candidate table and report are generated. Codon-dependent work is gated on a validated codon alignment, then native site/branch/gene-level selection outputs flow into the evidence ledger without converting branch or gene results into residue hits.

Every managed stage records status, inputs, outputs, parameters, software version, fingerprint, runtime, and SHA-256 hashes for outputs. When an external tool returns logs, the stage manifest retains bounded stdout/stderr and command metadata. Resume verifies stage fingerprints and output hashes; a changed/corrupt output is recomputed. `--dry-run` validates and reports the input plan before external alignment.

The HMMER job may overlap with downstream analysis when enabled. Most lightweight calculations run serially. The configured worker count applies to independent HyPhy method processes. Optional stage errors are reported and do not erase successful independent outputs. Missing tools/databases are `SKIPPED` or `FAILED`, never successful scientific analyses.

See `pipeline_stages.json`, `run_summary.json`, `provenance.json`, and `report.html` in an output directory for the actual execution record.

The Streamlit interface uses this same pipeline, keeps uploaded runs under `.evotrace/ui-runs/`, shows stage-manifest progress, filters candidate rows, and downloads either the report or a ZIP of run outputs.
