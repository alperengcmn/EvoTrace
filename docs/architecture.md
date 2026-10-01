# Architecture

EvoTrace uses a `src/` package layout. `cli` and the Streamlit presentation layer call the same pipeline API; `io` and `validation` handle input; `alignment` wraps MAFFT; `conservation` computes descriptive metrics; `selection` contains pairwise dN/dS; `tools` adapts external phylogeny, selection, and domain software; `motifs`, `structure`, and `mapping` project annotations into shared coordinates; `evidence` builds per-site summaries; and `pipeline` records stage execution and provenance.

External tools are called through fixed argument lists with captured diagnostics, never shell-interpolated user commands. Optional dependencies are not imported unless selected. Scientific modules remain usable independently of the CLI.
