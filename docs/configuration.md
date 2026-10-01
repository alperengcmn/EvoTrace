# Configuration

The canonical nested YAML example is `configs/example.yaml`; equivalent JSON is in `configs/example.json`. Input and external database/file paths in configuration resolve relative to the configuration file. The top-level configuration accepts `input`, `alphabet`, `max_ambiguity`, `alignment`, `phylogeny`, `selection`, `domains`, `structure`, `pipeline`, `report`, and `motifs` sections. The loader also accepts legacy flat settings for compatibility.

Selection inference needs `selection.enabled: true` and either a validated `selection.codon_alignment` or protein input plus matching `selection.coding_sequences` for exact back-translation. `selection.engine` chooses `hyphy` or `paml`; HyPhy methods are listed in `selection.analyses`. HMMER needs `domains.enabled: true` and a local profile HMM database. Structure mapping uses a local `structure.pdb_path` or `structure.alphafold_accession`; network retrieval uses the public AlphaFold DB API.

PAML models are `branch`, `site`, and `branch_site`; branch and branch-site models require a foreground `#1` branch label in the selected tree. PAML currently accepts `selection.genetic_code: 1` (standard code) only and rejects unsupported codes before launching codeml.

Unknown or invalid choices should fail with a direct configuration error. Optional programs/databases that are unavailable are reported as skipped/failed, not successful analyses.
