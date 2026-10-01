"""Streamlit presentation layer over the same analyze() API as the CLI."""

import json
import concurrent.futures
import io
import time
import uuid
import zipfile
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

from evotrace.pipeline.analyze import analyze

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RUNS_ROOT = PROJECT_ROOT / ".evotrace" / "ui-runs"
RUNS_ROOT.mkdir(parents=True, exist_ok=True)

PAGES = [
    "Home",
    "Input",
    "Configuration",
    "QC",
    "Alignment",
    "Phylogeny",
    "Conservation",
    "Selection",
    "Domains",
    "Motifs",
    "Structure",
    "Candidate Sites",
    "Integrated View",
    "Report",
]

st.set_page_config(page_title="EvoTrace", layout="wide")
st.title("EvoTrace · comparative sequence analysis")
st.caption("A local analysis interface. Scientific calculations use the same pipeline as the CLI.")
page = st.sidebar.radio("Analysis workspace", PAGES)

if page in ("Home", "Input", "Configuration"):
    st.subheader(page)
    if page == "Home":
        st.write(
            "Upload homologous sequences, review the pipeline configuration, then run the analysis. Results are descriptive unless a method-specific inferential test is configured."
        )
        st.code(
            "evotrace analyze --input sequences.fasta --config config.yaml --output results/run"
        )
    fasta = st.file_uploader(
        "FASTA sequences", type=["fa", "fasta", "fna", "faa", "txt"], key="fasta"
    )
    config = st.file_uploader(
        "Optional JSON/YAML configuration", type=["json", "yaml", "yml"], key="config"
    )
    threads = st.number_input("MAFFT threads", min_value=1, max_value=64, value=1)
    if fasta and st.button("Validate and analyze", type="primary"):
        work = RUNS_ROOT / time.strftime("%Y%m%d-%H%M%S") / uuid.uuid4().hex[:8]
        work.mkdir(parents=True, exist_ok=False)
        input_path, output_path = work / Path(fasta.name).name, work / "results"
        input_path.write_bytes(fasta.getvalue())
        config_path = None
        if config:
            config_path = work / Path(config.name).name
            config_path.write_bytes(config.getvalue())
        try:
            with st.status("Running EvoTrace workflow", expanded=True) as status:
                progress = st.progress(0, text="Validating input")
                progress_value = 0
                with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
                    future = pool.submit(
                        analyze, str(input_path), str(output_path),
                        str(config_path) if config_path else None, threads=int(threads)
                    )
                    while not future.done():
                        state_path = output_path / "pipeline_stages.json"
                        stages = {}
                        if state_path.exists():
                            try:
                                stages = json.loads(state_path.read_text()).get("stages", {})
                            except (OSError, ValueError):
                                pass
                        finished = sum(s.get("status") in ("SUCCESS", "SKIPPED") for s in stages.values())
                        fraction = min(0.92, finished / max(1, len(stages)))
                        progress_value = max(progress_value, int(fraction * 100))
                        current = next((n for n, s in stages.items() if s.get("status") == "RUNNING"), "Preparing workflow")
                        progress.progress(progress_value, text="{} · {} / {} stages complete".format(current, finished, len(stages)))
                        time.sleep(0.35)
                    summary = future.result()
                progress.progress(100, text="All configured stages finished")
                status.update(label="Analysis finished", state="complete")
            st.session_state["evotrace_results"] = str(output_path)
            st.session_state["evotrace_summary"] = summary
        except Exception as exc:
            st.error("Analysis failed: {}".format(exc))

root_value = st.session_state.get("evotrace_results")
if not root_value:
    if page not in ("Home", "Input", "Configuration"):
        st.info("Run an analysis from the Input page to populate this view.")
    st.stop()

root = Path(root_value)
summary = st.session_state.get("evotrace_summary", {})
stats = json.loads((root / "alignment_stats.json").read_text(encoding="utf-8"))
alignment = (root / "alignment.fasta").read_text(encoding="utf-8")
metrics = pd.read_csv(root / "alignment_metrics.csv")
candidates = pd.read_csv(root / "candidate_sites.csv")
st.caption(
    "{} sequences · {} alignment columns · {} invariant columns".format(
        stats.get("sequence_count"), stats.get("alignment_length"), stats.get("conserved_columns")
    )
)

if page == "Home":
    cols = st.columns(4)
    cols[0].metric("Sequences", stats["sequence_count"])
    cols[1].metric("Alignment columns", stats["alignment_length"])
    cols[2].metric("Invariant columns", stats["conserved_columns"])
    cols[3].metric("Mean identity", "{:.1%}".format(stats["mean_pairwise_identity"]))
    st.json(summary)
elif page == "Input":
    st.subheader("Validated input and aligned records")
    st.text(alignment[:30000])
elif page == "Configuration":
    st.subheader("Configuration and provenance")
    st.json(json.loads((root / "provenance.json").read_text(encoding="utf-8")))
elif page == "QC":
    for name in ("validation.json", "alignment_qc.json", "codon_alignment_qc.json"):
        path = root / name
        if path.exists():
            st.subheader(name.replace(".json", ""))
            st.json(json.loads(path.read_text(encoding="utf-8")))
    st.json(summary.get("stage_status", {}))
elif page == "Alignment":
    st.text(alignment[:50000])
    st.download_button(
        "Download aligned FASTA", (root / "alignment.fasta").read_bytes(), "alignment.fasta"
    )
elif page == "Phylogeny":
    svg = root / "phylogeny.svg"
    if svg.exists():
        st.image(str(svg), use_container_width=True)
    nwk = root / "phylogeny.nwk"
    if nwk.exists():
        st.code(nwk.read_text(encoding="utf-8"), language="text")
    st.json(json.loads((root / "phylogeny_summary.json").read_text(encoding="utf-8")))
elif page == "Conservation":
    st.image(str(root / "conservation.svg"), use_container_width=True)
    st.dataframe(metrics, use_container_width=True)
elif page == "Selection":
    path = root / "selection.json"
    st.json(json.loads(path.read_text(encoding="utf-8")) if path.exists() else [])
    st.info(
        "Pairwise dN/dS is descriptive; HyPhy and PAML tests have method-specific statistical interpretations."
    )
elif page in ("Domains", "Motifs", "Structure"):
    names = {
        "Domains": "domains.json",
        "Motifs": "motif_hits.json",
        "Structure": "structure_mapping.json",
    }
    path = root / names[page]
    result = (
        json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"status": "unavailable"}
    )
    st.json(result)
elif page == "Candidate Sites":
    max_position = int(candidates["alignment_position"].max()) if not candidates.empty else 1
    start, end = st.slider(
        "Alignment position range", 1, max(1, max_position), (1, max(1, max_position))
    )
    filtered = candidates[candidates["alignment_position"].between(start, end)]
    domains_available = sorted(x for x in candidates.get("domain", pd.Series(dtype=str)).dropna().unique() if str(x).strip())
    selected_domains = st.multiselect("Domains", domains_available)
    if selected_domains:
        filtered = filtered[filtered["domain"].isin(selected_domains)]
    selection_methods = sorted(x for x in candidates.get("selection_method", pd.Series(dtype=str)).dropna().unique() if str(x).strip())
    selected_methods = st.multiselect("Selection methods", selection_methods)
    if selected_methods:
        filtered = filtered[filtered["selection_method"].fillna("").apply(lambda v: any(m in v.split(";") for m in selected_methods))]
    min_score = float(candidates.get("ets_descriptive_score", pd.Series([0.0])).fillna(0).min()) if not candidates.empty else 0.0
    cutoff = st.slider("Minimum descriptive ETS proxy", min_value=0.0, max_value=1.0, value=min_score, step=0.01)
    if "ets_descriptive_score" in filtered:
        filtered = filtered[filtered["ets_descriptive_score"].fillna(0) >= cutoff]
    st.dataframe(filtered, use_container_width=True)
    st.download_button(
        "Download filtered candidate sites",
        filtered.to_csv(index=False),
        "candidate_sites.csv",
        "text/csv",
    )
elif page == "Integrated View":
    st.subheader("Position-level evidence map")
    st.caption(
        "Sequence variability and annotations share the same alignment coordinate. Missing annotations remain missing."
    )
    plot = px.scatter(
        metrics,
        x="position",
        y="conservation_score",
        color="entropy",
        hover_data=["dominant_residue", "gap_fraction", "occupancy"],
        labels={
            "position": "Alignment position",
            "conservation_score": "Conservation",
            "entropy": "Entropy (bits)",
        },
        title="Conservation, entropy, and mapped candidate annotations",
    )
    plot.add_scatter(
        x=candidates.loc[candidates["domain"].notna(), "alignment_position"],
        y=candidates.loc[candidates["domain"].notna(), "conservation"],
        mode="markers",
        name="Domain overlap",
        marker={"symbol": "diamond-open", "size": 12},
    )
    plot.add_scatter(
        x=candidates.loc[candidates["motif"].notna(), "alignment_position"],
        y=candidates.loc[candidates["motif"].notna(), "conservation"],
        mode="markers",
        name="Motif overlap",
        marker={"symbol": "x", "size": 12},
    )
    if "selection_pvalue" in candidates:
        selected_sites = candidates[candidates["selection_pvalue"].notna()]
        plot.add_scatter(
            x=selected_sites["alignment_position"],
            y=selected_sites["conservation"],
            mode="markers",
            name="Selection evidence (raw p reported)",
            marker={"symbol": "star", "size": 13, "color": "#bd3a3a"},
            customdata=selected_sites[["selection_method", "selection_pvalue"]],
            hovertemplate="Position %{x}<br>Conservation %{y:.3f}<br>%{customdata[0]} p=%{customdata[1]}<extra></extra>",
        )
    if "structure_residue" in candidates:
        mapped_sites = candidates[candidates["structure_residue"].notna()]
        plot.add_scatter(
            x=mapped_sites["alignment_position"],
            y=mapped_sites["conservation"],
            mode="markers",
            name="Structure mapped",
            marker={"symbol": "triangle-up-open", "size": 12, "color": "#3d6688"},
            customdata=mapped_sites[["structure_residue", "structure_confidence", "structure_b_factor"]],
            hovertemplate="Position %{x}<br>Conservation %{y:.3f}<br>%{customdata[0]}<br>AlphaFold confidence %{customdata[1]}<br>Mean B-factor %{customdata[2]}<extra></extra>",
        )
    st.plotly_chart(plot, use_container_width=True)
    position = st.number_input(
        "Inspect alignment position",
        min_value=1,
        max_value=max(1, int(metrics["position"].max())),
        value=1,
    )
    selected = candidates[candidates["alignment_position"] == position]
    if not selected.empty:
        st.json(selected.iloc[0].to_dict())
elif page == "Report":
    st.components.v1.html(
        (root / "report.html").read_text(encoding="utf-8"), height=1000, scrolling=True
    )
    st.download_button(
        "Download HTML report", (root / "report.html").read_bytes(), "report.html", "text/html"
    )

with st.sidebar.expander("Download outputs"):
    for file in sorted(root.iterdir()):
        if file.is_file() and file.name != "report.html":
            st.download_button(file.name, file.read_bytes(), file.name, key="download-" + file.name)
    bundle = io.BytesIO()
    with zipfile.ZipFile(bundle, "w", zipfile.ZIP_DEFLATED) as archive:
        for file in root.rglob("*"):
            if file.is_file():
                archive.write(file, file.relative_to(root))
    st.download_button("Download all run outputs (.zip)", bundle.getvalue(), "evotrace-results.zip", "application/zip")
