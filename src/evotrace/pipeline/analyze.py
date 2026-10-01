"""Content-addressed, scientifically explicit evolutionary analysis workflow."""

import csv
import hashlib
import json
import logging
import platform
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from Bio import Phylo

from evotrace import __version__
from evotrace.alignment.mafft import run_mafft
from evotrace.config import load_config
from evotrace.conservation.core import alignment_stats, columns, identity_matrix
from evotrace.codon.backtranslate import validate_codon_alignment
from evotrace.evidence.sites import candidate_sites, write_candidate_sites
from evotrace.io.fasta import read_fasta, write_fasta
from evotrace.io.ncbi import fetch_accessions
from evotrace.mapping.residues import AlignmentMap
from evotrace.motifs.scan import scan_motifs
from evotrace.phylogeny.nj import neighbor_joining
from evotrace.pipeline.engine import PipelineEngine
from evotrace.reporting.html import render_report
from evotrace.scoring.ets import evolutionary_scores
from evotrace.selection.ng import pairwise_dnds
from evotrace.structure.providers import AlphaFoldProvider, PDBProvider, map_sequence_to_structure
from evotrace.tools.common import ToolUnavailable
from evotrace.tools.hmmer import run_hmmscan
from evotrace.tools.hyphy import run_hyphy
from evotrace.tools.iqtree import run_iqtree
from evotrace.tools.paml import run_codeml
from evotrace.validation.fasta import validate_records
from evotrace.visualization.svg import profile_svg

log = logging.getLogger("evotrace")


def _sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True, default=str), encoding="utf-8")


def _configured_path(value, config_path):
    if not value:
        return None
    path = Path(value).expanduser()
    if not path.is_absolute() and config_path:
        path = Path(config_path).resolve().parent / path
    return path.resolve()


def _tool_version(name):
    exe = shutil.which(name)
    if not exe:
        return "not-installed"
    if Path(exe).name == "codeml":
        return "installed; codeml has no non-interactive version flag"
    version_flag = "-version" if Path(exe).name in ("iqtree", "iqtree2", "iqtree3") else "--version"
    try:
        proc = subprocess.run(
            [exe, version_flag], capture_output=True, text=True, timeout=8, check=False
        )
        return ((proc.stdout or proc.stderr).strip().splitlines() or ["unknown"])[0]
    except Exception as exc:
        return "version-query-failed: {}".format(exc)


def _draw_tree_svg(tree, output):
    """Render a compact rectangular phylogram with branch lengths and supports."""
    terminals = tree.get_terminals()
    if not terminals:
        raise ValueError("Tree contains no terminal nodes.")
    max_depth = max((tree.distance(t) for t in terminals), default=0.0) or 1.0
    width, row_h, left, right, top = 1100, 30, 70, 240, 35
    height = max(100, top * 2 + row_h * len(terminals))
    ys = {t: top + row_h * (i + 0.5) for i, t in enumerate(terminals)}
    depths = tree.depths()
    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" width="100%">'.format(
            width, height
        ),
        '<rect width="100%" height="100%" fill="white"/>',
    ]

    def x(node):
        return left + (width - left - right) * depths.get(node, 0.0) / max_depth

    def y(node):
        if node in ys:
            return ys[node]
        vals = [y(child) for child in node.clades]
        return sum(vals) / len(vals) if vals else top

    for node in tree.get_nonterminals(order="postorder"):
        child_y = [y(child) for child in node.clades]
        parts.append(
            '<path d="M {:.2f} {:.2f} V {:.2f}" stroke="#55727c" fill="none"/>'.format(
                x(node), min(child_y), max(child_y)
            )
        )
        for child in node.clades:
            parts.append(
                '<path d="M {:.2f} {:.2f} H {:.2f}" stroke="#55727c" fill="none"/>'.format(
                    x(node), y(child), x(child)
                )
            )
    for node in terminals:
        # terminal label coordinates derive from its cumulative branch depth
        parts.append(
            '<text x="{:.2f}" y="{:.2f}" font-family="sans-serif" font-size="13">{}</text>'.format(
                x(node) + 5, y(node) + 4, __import__("html").escape(node.name or "unnamed")
            )
        )
    for node in tree.get_nonterminals():
        if node.confidence is not None:
            parts.append(
                '<text x="{:.2f}" y="{:.2f}" font-family="sans-serif" font-size="10" fill="#9a4c13">{}</text>'.format(
                    x(node) + 3, y(node) - 3, __import__("html").escape(str(node.confidence))
                )
            )
    parts.append("</svg>")
    Path(output).write_text("\n".join(parts), encoding="utf-8")


def _is_codon_alignment(alignment):
    try:
        validate_codon_alignment(alignment)
        return True, ""
    except ValueError as exc:
        return False, str(exc)


def analyze(
    input_path=None,
    output_dir="results",
    config_path=None,
    dry_run=False,
    threads=1,
    resume=False,
    workers=1,
    accessions=None,
):
    started = time.time()
    cfg = load_config(config_path)
    if workers == 1 and int(cfg.get("pipeline", {}).get("workers", 1)) > 1:
        workers = int(cfg["pipeline"]["workers"])
    out = Path(output_dir).expanduser().resolve()
    input_cfg = cfg.get("input_config", {})
    if input_path and accessions:
        raise ValueError("Choose either --input or --accession, not both.")
    if input_path is None and not accessions:
        input_path = input_cfg.get("path")
        if input_path and config_path and not Path(input_path).expanduser().is_absolute():
            input_path = Path(config_path).resolve().parent / input_path
        if not input_path and input_cfg.get("type") == "accession":
            accessions = input_cfg.get("ids", input_cfg.get("accessions", []))
    retrieval_info = None
    out.mkdir(parents=True, exist_ok=True)
    if accessions:
        downloaded_path = out / "downloaded_accessions.fasta"
        retrieval_info = fetch_accessions(
            accessions,
            downloaded_path,
            out / ".cache" / "ncbi",
            database=input_cfg.get("database", "protein"),
            cache_ttl_days=float(input_cfg.get("cache_ttl_days", 7)),
            timeout=int(input_cfg.get("timeout", 30)),
        )
        input_path = downloaded_path
    if not input_path:
        raise ValueError("Provide --input, --accession, input.path, or input.ids in configuration.")
    inp = Path(input_path).expanduser().resolve()
    if not inp.is_file():
        raise ValueError("Input FASTA does not exist or is not a file: {}".format(inp))
    if threads < 1 or workers < 1:
        raise ValueError("threads and workers must be positive integers.")
    engine = PipelineEngine(out, resume=resume)
    input_hash = _sha(inp)
    records = read_fasta(inp)
    alphabet, validation_warnings = validate_records(
        records, cfg.get("alphabet", "auto"), float(cfg.get("max_ambiguity", 0.5))
    )
    validation_path = out / "validation.json"
    if dry_run:
        print(
            "Validated {} {} sequences; planned MAFFT strategy {}; output {}".format(
                len(records), alphabet, cfg.get("mafft_strategy"), out
            )
        )
        return {"dry_run": True, "sequence_count": len(records), "alphabet": alphabet}
    validation = {
        "sequence_count": len(records),
        "alphabet": alphabet,
        "warnings": validation_warnings,
        "input_sha256": input_hash,
        "sequence_lengths": {r.id: len(r.sequence) for r in records},
    }
    engine.run_stage(
        "validation",
        lambda: _write_json(validation_path, validation),
        {"input_sha256": input_hash},
        [validation_path],
        parameters={"alphabet": cfg.get("alphabet"), "max_ambiguity": cfg.get("max_ambiguity")},
    )
    log.info("Validated %d %s sequences", len(records), alphabet)

    aligned_path = out / "alignment.fasta"
    mafft_version = _tool_version("mafft")
    command_holder = {}

    def align_action():
        command_holder["command"] = run_mafft(
            records,
            aligned_path,
            cfg.get("mafft_strategy", "auto"),
            threads,
            timeout=int(cfg.get("alignment_timeout", 86400)),
        )
        return command_holder["command"]

    align_result = engine.run_stage(
        "alignment",
        align_action,
        {"input_sha256": input_hash},
        [aligned_path],
        dependencies=["validation"],
        parameters={"strategy": cfg.get("mafft_strategy", "auto"), "threads": threads},
        software_version=mafft_version,
    )
    engine.state["stages"]["alignment"]["command"] = command_holder.get("command", align_result)
    engine._save()
    aligned = read_fasta(aligned_path)
    qc_path = out / "alignment_qc.json"

    def qc_alignment():
        identifiers = [record.id for record in aligned]
        lengths = {len(record.sequence) for record in aligned}
        if len(identifiers) != len(set(identifiers)):
            raise ValueError("Alignment QC failed: sequence identifiers are not unique.")
        if len(lengths) != 1 or not lengths or next(iter(lengths)) == 0:
            raise ValueError(
                "Alignment QC failed: aligned sequences must have one non-zero length."
            )
        if {record.id for record in aligned} != {record.id for record in records}:
            raise ValueError("Alignment QC failed: sequence IDs differ from the input dataset.")
        _write_json(
            qc_path,
            {
                "status": "PASS",
                "sequence_count": len(aligned),
                "alignment_length": len(aligned[0].sequence),
                "unique_ids": True,
                "input_ids_preserved": True,
            },
        )

    engine.run_stage(
        "alignment_qc",
        qc_alignment,
        {"alignment_sha256": _sha(aligned_path)},
        [qc_path],
        dependencies=["alignment"],
        parameters={"checks": ["unique_ids", "equal_nonzero_length", "input_ids_preserved"]},
    )
    stats = alignment_stats(aligned, alphabet, cfg.get("gap_policy", "ignore"))
    stats["alphabet"] = alphabet
    stats["sequence_lengths"] = {r.id: len(r.sequence.replace("-", "")) for r in records}
    col = columns(aligned, alphabet, cfg.get("gap_policy", "ignore"))
    matrix = identity_matrix(aligned, alphabet)
    scores = evolutionary_scores(col)
    warnings = list(validation_warnings)
    phylo_info, tree_path = {}, None
    phy_cfg = cfg.get("phylogeny", {})
    if len(aligned) < 3:
        engine.skip_stage(
            "phylogeny",
            "At least three aligned records are required.",
            dependencies=["alignment_qc"],
        )
        warnings.append("Phylogeny was skipped: at least three sequences are required.")
    else:
        iqtree_requested = phy_cfg.get("tool", "nj") == "iqtree"
        if iqtree_requested:
            iqtree_dir = out / "iqtree"
            iqtree_tree_path = iqtree_dir / "iqtree.treefile"
            iqtree_report_path = iqtree_dir / "iqtree.iqtree"
            try:
                iqtree_holder = {}

                def iqtree_action():
                    iqtree_holder.update(
                        run_iqtree(
                            aligned_path,
                            iqtree_dir,
                            phy_cfg.get("model", "MFP"),
                            int(phy_cfg.get("bootstrap", 1000)),
                            threads,
                        )
                    )
                    return {
                        key: iqtree_holder.get(key)
                        for key in (
                            "command",
                            "working_directory",
                            "runtime_seconds",
                            "exit_code",
                            "version",
                            "stdout",
                            "stderr",
                        )
                    }

                engine.run_stage(
                    "iqtree",
                    iqtree_action,
                    {"alignment_sha256": _sha(aligned_path)},
                    [iqtree_tree_path, iqtree_report_path],
                    dependencies=["alignment_qc"],
                    parameters={
                        "model": phy_cfg.get("model", "MFP"),
                        "bootstrap": phy_cfg.get("bootstrap", 1000),
                        "threads": threads,
                    },
                    software_version=_tool_version("iqtree2")
                    if shutil.which("iqtree2")
                    else _tool_version("iqtree"),
                )
                engine.state["stages"]["iqtree"]["command"] = iqtree_holder.get("command")
                engine._save()
                tree_path = iqtree_tree_path
                phylo_info = {
                    k: v for k, v in iqtree_holder.items() if k not in ("tree", "stdout", "stderr")
                }
                if not phylo_info:
                    import re

                    report_text = iqtree_report_path.read_text(errors="replace")
                    ll_match = re.search(
                        r"Log-likelihood of the tree:\s*(-?\d+(?:\.\d+)?)", report_text
                    )
                    model_match = re.search(
                        r"Best-fit model according to .*?:\s*(\S+)", report_text
                    )
                    parsed_tree = Phylo.read(str(iqtree_tree_path), "newick")
                    phylo_info = {
                        "log_likelihood": float(ll_match.group(1)) if ll_match else None,
                        "selected_model": model_match.group(1)
                        if model_match
                        else phy_cfg.get("model", "MFP"),
                        "support_values_present": any(
                            t.confidence is not None for t in parsed_tree.get_nonterminals()
                        ),
                    }
                phylo_info["method"] = "IQ-TREE maximum likelihood"
            except ToolUnavailable as exc:
                engine.skip_stage("iqtree", str(exc), dependencies=["alignment"])
                warnings.append(str(exc) + " Falling back to explicitly labeled Neighbor Joining.")
            except Exception as exc:
                engine.state["stages"].setdefault("iqtree", {"status": "FAILED", "error": str(exc)})
                engine._save()
                warnings.append(
                    "IQ-TREE failed: {}. Neighbor Joining fallback is reported separately.".format(
                        exc
                    )
                )
        if tree_path is None:
            nj_path, svg_path = out / "phylogeny.nwk", out / "phylogeny.svg"

            def nj_action():
                tree = neighbor_joining(aligned)
                Phylo.write(tree, str(nj_path), "newick")
                _draw_tree_svg(tree, svg_path)

            engine.run_stage(
                "neighbor_joining",
                nj_action,
                {"alignment_sha256": _sha(aligned_path)},
                [nj_path, svg_path],
                dependencies=["alignment_qc"],
                parameters={"distance": "pairwise uncorrected p-distance", "support": "none"},
                software_version=__import__("Bio").__version__,
            )
            tree_path = nj_path
            phylo_info = {
                "method": "Neighbor Joining on uncorrected p-distance",
                "support": "No bootstrap or other branch support calculated",
                "tree_path": str(nj_path),
            }
        stats["phylogeny_method"] = phylo_info.get("method", "not run")
        stats["phylogeny_support"] = (
            "IQ-TREE support values present"
            if phylo_info.get("support_values_present")
            else phylo_info.get("support", "No support values calculated")
        )
        stats["phylogeny_log_likelihood"] = phylo_info.get("log_likelihood")

    # Domain scanning does not depend on the alignment/tree; overlap it with later work.
    domain_cfg = cfg.get("domains_config", {})
    domain_executor = None
    domain_future = None
    domain_db = _configured_path(domain_cfg.get("database"), config_path)
    if domain_cfg.get("enabled") and alphabet == "protein" and domain_db and not resume:
        domain_executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="evotrace-hmmer")
        domain_future = domain_executor.submit(
            run_hmmscan, inp, Path(domain_db), out / "domains", threads
        )

    # Alignment-aware motifs are mapped through the shared coordinate abstraction.
    motif_config = cfg.get("motifs", [])
    motif_hits = scan_motifs(records, motif_config, alignment=aligned)
    selection_results, selection_hits = [], []
    selection_cfg = cfg.get("selection_config", {})
    selection_tree_path = tree_path
    selection_tree_value = selection_cfg.get("tree")
    if selection_tree_value:
        selection_tree_path = Path(selection_tree_value).expanduser()
        if not selection_tree_path.is_absolute() and config_path:
            selection_tree_path = Path(config_path).resolve().parent / selection_tree_path
        selection_tree_path = selection_tree_path.resolve()
        if not selection_tree_path.is_file():
            raise ValueError("selection.tree does not exist: {}".format(selection_tree_path))
    selection_alignment = None
    selection_alignment_path = None
    codon_source = None
    codon_path_value = selection_cfg.get("codon_alignment")
    cds_path_value = selection_cfg.get("coding_sequences")
    if codon_path_value:
        codon_path = Path(codon_path_value).expanduser()
        if not codon_path.is_absolute() and config_path:
            codon_path = Path(config_path).resolve().parent / codon_path
        selection_alignment = read_fasta(codon_path)
        if {r.id for r in selection_alignment} != {r.id for r in aligned}:
            raise ValueError(
                "selection.codon_alignment IDs must match the analysis alignment IDs exactly."
            )
        ok, reason = _is_codon_alignment(selection_alignment)
        if not ok:
            raise ValueError("Invalid selection.codon_alignment: " + reason)
        selection_alignment_path = codon_path.resolve()
        codon_source = "user-provided codon-aware alignment"
        if alphabet == "protein":
            try:
                validate_codon_alignment(
                    selection_alignment, aligned, int(selection_cfg.get("genetic_code", 1))
                )
            except ValueError as exc:
                raise ValueError("Invalid codon/protein correspondence: {}".format(exc)) from exc
    elif alphabet == "protein" and cds_path_value:
        from evotrace.codon.backtranslate import backtranslate

        cds_path = Path(cds_path_value).expanduser()
        if not cds_path.is_absolute() and config_path:
            cds_path = Path(config_path).resolve().parent / cds_path
        selection_alignment, codon_meta = backtranslate(
            aligned, read_fasta(cds_path), int(selection_cfg.get("genetic_code", 1))
        )
        selection_alignment_path = out / "codon_alignment.fasta"
        write_fasta(selection_alignment, selection_alignment_path)
        _write_json(out / "codon_alignment_validation.json", codon_meta)
        codon_source = (
            "back-translated from the protein alignment after exact translation validation"
        )

    if selection_alignment is not None:
        codon_ok, codon_reason = _is_codon_alignment(selection_alignment)
        codon_report = validate_codon_alignment(
            selection_alignment,
            aligned if alphabet == "protein" else None,
            int(selection_cfg.get("genetic_code", 1)),
        ) if codon_ok else {}
        codon_length = len(selection_alignment[0].sequence) if selection_alignment else 0
        _write_json(
            out / "codon_alignment_qc.json",
            {
                "status": "PASS" if codon_ok else "FAIL",
                "sequence_count": len(selection_alignment),
                "codon_columns": codon_length // 3,
                "gap_containing_codons": sum(
                    1
                    for record in selection_alignment
                    for offset in range(0, len(record.sequence), 3)
                    if "-" in record.sequence[offset : offset + 3]
                ),
                "invalid_codons": 0 if codon_ok else 1,
                "premature_stops": 0 if codon_ok else 1,
                "frame_violations": 0 if codon_ok else 1,
                "reason": codon_reason or None,
                "source": codon_source,
                **codon_report,
            },
        )

    if cfg.get("pairwise_dnds"):
        if selection_alignment is None:
            selection_results.append(
                {
                    "status": "not_run",
                    "reason": "Pairwise dN/dS requires a codon-aware alignment. Provide selection.codon_alignment or provide protein input plus selection.coding_sequences.",
                }
            )
        elif not _is_codon_alignment(selection_alignment)[0]:
            selection_results.append(
                {"status": "not_run", "reason": _is_codon_alignment(selection_alignment)[1]}
            )
        else:
            selection_results.extend(pairwise_dnds(selection_alignment))
    if cfg.get("selection"):
        requested_engine = selection_cfg.get("engine", "hyphy")
        if selection_alignment is None:
            engine.skip_stage(
                "selection",
                "No validated codon alignment was configured.",
                dependencies=["alignment_qc"],
            )
            selection_results.append(
                {
                    "status": "skipped",
                    "reason": "HyPhy/PAML selection requires a validated codon alignment. Configure selection.codon_alignment or provide selection.coding_sequences with a protein alignment.",
                }
            )
        else:
            ok, reason = _is_codon_alignment(selection_alignment)
            if not ok:
                engine.skip_stage("selection", reason, dependencies=["alignment_qc"])
                selection_results.append(
                    {
                        "status": "skipped",
                        "reason": reason,
                    }
                )
            else:
                analyses = selection_cfg.get("analyses", ["FEL"])
                if not isinstance(analyses, list) or not analyses:
                    raise ValueError("selection.analyses must be a non-empty list.")
                if requested_engine == "hyphy":
                    from evotrace.tools.hyphy import METHODS, _method_name

                    unknown_methods = [x for x in analyses if _method_name(x) not in METHODS]
                    if unknown_methods:
                        raise ValueError(
                            "Unsupported HyPhy method(s): {}.".format(", ".join(map(str, unknown_methods)))
                        )
                    method_names = [_method_name(x) for x in analyses]
                    if len(set(method_names)) != len(method_names):
                        raise ValueError("selection.analyses contains duplicate methods.")
                if requested_engine == "hyphy":
                    try:
                        selection_folder = out / "selection"

                        def selection_action():
                            run_records = []
                            with ThreadPoolExecutor(
                                max_workers=min(workers, len(analyses)),
                                thread_name_prefix="evotrace-hyphy",
                            ) as pool:
                                futures = [
                                    pool.submit(
                                        run_hyphy,
                                        method,
                                        selection_alignment_path,
                                        selection_folder / str(method).lower(),
                                        selection_tree_path,
                                    )
                                    for method in analyses
                                ]
                                for method, future in zip(analyses, futures):
                                    item = future.result()
                                    run_records.append(
                                        {
                                            key: item.get(key)
                                            for key in (
                                                "method",
                                                "command",
                                                "working_directory",
                                                "exit_code",
                                                "runtime_seconds",
                                                "version",
                                                "stdout",
                                                "stderr",
                                            )
                                        }
                                    )
                                    _write_json(
                                        selection_folder / (str(method).lower() + ".run.json"),
                                        run_records[-1],
                                    )
                            return {
                                "command": [record["command"] for record in run_records],
                                "working_directory": str(out / "selection"),
                                "runtime_seconds": sum(
                                    record.get("runtime_seconds") or 0 for record in run_records
                                ),
                                "stdout": "\n".join(
                                    record.get("stdout") or "" for record in run_records
                                ),
                                "stderr": "\n".join(
                                    record.get("stderr") or "" for record in run_records
                                ),
                                "version": _tool_version("hyphy"),
                            }

                        engine.run_stage(
                            "selection",
                            selection_action,
                            {"alignment_sha256": _sha(selection_alignment_path)},
                            [
                                path
                                for m in analyses
                                for path in (
                                    selection_folder / str(m).lower() / (str(m).lower() + ".json"),
                                    selection_folder / (str(m).lower() + ".run.json"),
                                )
                            ],
                            dependencies=["alignment_qc"],
                            parameters={
                                "engine": "hyphy",
                                "analyses": analyses,
                                "codon_source": codon_source,
                            },
                            software_version=_tool_version("hyphy"),
                        )
                        engine.state["stages"]["selection"]["command"] = [
                            "hyphy {} --alignment ...".format(str(m).lower()) for m in analyses
                        ]
                        engine._save()
                        from evotrace.tools.hyphy import parse_hyphy

                        for method in analyses:
                            method_name = str(method).lower()
                            path = selection_folder / method_name / (method_name + ".json")
                            parsed = parse_hyphy(path, str(method))
                            selection_results.append(
                                {k: v for k, v in parsed.items() if k != "raw"}
                                | {"path": str(path)}
                            )
                            if alphabet == "protein":
                                selection_hits.extend(
                                    {
                                        **site,
                                        "alignment_position": site.get("site"),
                                        "method": method,
                                        "test": method,
                                    }
                                    for site in parsed.get("sites", [])
                                )
                    except ToolUnavailable as exc:
                        engine.skip_stage("selection", str(exc), dependencies=["alignment_qc"])
                        selection_results.append(
                            {"status": "unavailable", "engine": "hyphy", "reason": str(exc)}
                        )
                        warnings.append(str(exc))
                    except Exception as exc:
                        selection_results.append(
                            {"status": "failed", "engine": "hyphy", "reason": str(exc)}
                        )
                        warnings.append("HyPhy analysis failed: {}".format(exc))
                elif requested_engine == "paml":
                    if not selection_tree_path:
                        engine.skip_stage(
                            "selection",
                            "PAML/codeml selection requires a tree.",
                            dependencies=["alignment_qc"],
                        )
                        selection_results.append(
                            {"status": "skipped", "reason": "PAML/codeml requires a tree."}
                        )
                    else:
                        try:
                            model = selection_cfg.get("model", "branch_site")
                            paml_json = out / "selection" / "paml_result.json"

                            def paml_action():
                                item = run_codeml(
                                    selection_alignment_path,
                                    selection_tree_path,
                                    out / "selection" / "paml",
                                    model,
                                    genetic_code=int(selection_cfg.get("genetic_code", 1)),
                                )
                                compact = {
                                    k: v
                                    for k, v in item.items()
                                    if k not in ("raw_output", "null_raw_output")
                                }
                                _write_json(paml_json, compact)
                                return {
                                    key: compact.get(key)
                                    for key in (
                                        "command", "working_directory", "runtime_seconds",
                                        "exit_code", "version", "stdout", "stderr",
                                    )
                                }

                            engine.run_stage(
                                "selection",
                                paml_action,
                                {
                                    "codon_alignment_sha256": _sha(selection_alignment_path),
                                    "tree_sha256": _sha(selection_tree_path),
                                },
                                [
                                    paml_json,
                                    out / "selection" / "paml" / "codeml_alignment.phy",
                                    out / "selection" / "paml" / "codeml_tree.trees",
                                    out / "selection" / "paml" / "codeml.ctl",
                                    out / "selection" / "paml" / "codeml.out",
                                    out / "selection" / "paml" / "codeml_null.ctl",
                                    out / "selection" / "paml" / "codeml_null.out",
                                ],
                                dependencies=["alignment_qc"],
                                parameters={"engine": "paml", "model": model, "codon_source": codon_source, "genetic_code": int(selection_cfg.get("genetic_code", 1))},
                                software_version=_tool_version("codeml"),
                            )
                            paml_result = json.loads(paml_json.read_text(encoding="utf-8"))
                            selection_results.append(paml_result)
                            if alphabet == "protein":
                                selection_hits.extend(
                                    {
                                        **site,
                                        "alignment_position": site.get("site"),
                                        "method": "PAML/codeml",
                                        "test": "PAML/codeml",
                                        "selection_class": "codeml_BEB_posterior",
                                    }
                                    for site in paml_result.get("site_posterior_results", [])
                                )
                        except Exception as exc:
                            selection_results.append(
                                {
                                    "status": "unavailable_or_failed",
                                    "engine": "paml",
                                    "reason": str(exc),
                                }
                            )
                            warnings.append("PAML was not run successfully: {}".format(exc))

    domain_hits = []
    if domain_cfg.get("enabled"):
        db = _configured_path(domain_cfg.get("database"), config_path)
        if alphabet != "protein":
            reason = "HMMER domain search requires protein sequences; translate/select protein input first."
            engine.skip_stage("domains", reason, dependencies=["validation"])
            warnings.append(reason)
        elif not db:
            reason = "HMMER domain scan was requested but domains.database was not configured."
            engine.skip_stage("domains", reason, dependencies=["validation"])
            warnings.append(reason)
        else:
            try:
                domain_path = out / "domains" / "domains.domtblout"
                domain_holder = {}

                def domain_action():
                    if domain_future is not None:
                        domain_holder.update(domain_future.result())
                    else:
                        domain_holder.update(run_hmmscan(inp, Path(db), out / "domains", threads))
                    _write_json(
                        out / "domains.json",
                        {k: v for k, v in domain_holder.items() if k not in ("stdout", "stderr")},
                    )
                    return {
                        key: domain_holder.get(key)
                        for key in (
                            "command",
                            "working_directory",
                            "runtime_seconds",
                            "exit_code",
                            "version",
                            "stdout",
                            "stderr",
                        )
                    }

                engine.run_stage(
                    "domains",
                    domain_action,
                    {"sequence_sha256": input_hash, "database_sha256": _sha(db)},
                    [domain_path, out / "domains.json"],
                    dependencies=["validation"],
                    parameters={"database": str(db), "threads": threads},
                    software_version=_tool_version("hmmscan"),
                )
                engine.state["stages"]["domains"]["command"] = domain_holder.get("command")
                engine._save()
                if domain_holder:
                    domain_hits = domain_holder["hits"]
                else:
                    from evotrace.tools.hmmer import parse_domtblout

                    domain_hits = parse_domtblout(domain_path)
            except ToolUnavailable as exc:
                reason = str(exc)
                engine.skip_stage("domains", reason, dependencies=["validation"])
                warnings.append(reason)
            except Exception as exc:
                reason = str(exc)
                if (
                    "domains" not in engine.state["stages"]
                    or engine.state["stages"]["domains"].get("status") != "FAILED"
                ):
                    engine.state["stages"]["domains"] = {
                        "name": "domains",
                        "status": "FAILED",
                        "error": reason,
                    }
                    engine._save()
                warnings.append(reason)
            finally:
                if domain_executor is not None:
                    domain_executor.shutdown(wait=True)
    else:
        engine.skip_stage("domains", "Disabled in configuration.", dependencies=["validation"])

    structure_cfg = cfg.get("structure_config", {})
    structure_result = {
        "status": "structure_mapping_disabled",
        "reason": "Disabled in configuration.",
    }
    structure_path_json = out / "structure_mapping.json"
    if structure_cfg.get("enabled"):
        reference = next(r for r in aligned if r.id == (cfg.get("reference") or aligned[0].id))
        query_sequence = reference.sequence.replace("-", "").replace("?", "")
        structure_holder = {}
        try:
            if alphabet != "protein":
                raise ValueError("Protein structure mapping requires protein sequence input.")
            pdb_path = _configured_path(structure_cfg.get("pdb_path"), config_path)
            accession = structure_cfg.get("alphafold_accession")
            if pdb_path:
                provider = PDBProvider(Path(pdb_path))
                source_key = (
                    {"local_pdb_sha256": _sha(pdb_path)}
                    if Path(pdb_path).is_file()
                    else {"path": str(pdb_path)}
                )
                tool_version = "local PDB/mmCIF"
            elif accession:
                provider = AlphaFoldProvider(
                    str(accession),
                    out / ".cache" / "alphafold",
                    float(structure_cfg.get("cache_ttl_days", 30)),
                    int(structure_cfg.get("timeout", 30)),
                )
                source_key = {"alphafold_accession": str(accession)}
                tool_version = "AlphaFold DB API"
            else:
                raise ValueError("Set structure.pdb_path or structure.alphafold_accession.")

            def structure_action():
                structure_holder["source"] = provider.get_structure()
                source = structure_holder["source"]
                if source.get("status") == "structure_found":
                    structure_holder["mapping"] = map_sequence_to_structure(
                        query_sequence, Path(source["path"]), structure_cfg.get("chain")
                    )
                    structure_holder["mapping"]["source_provider"] = source.get(
                        "metadata", {}
                    ).get("provider", "local PDB/mmCIF")
                    structure_holder["mapping"]["structure_id"] = (
                        structure_cfg.get("alphafold_accession") or Path(source["path"]).stem
                    )
                else:
                    structure_holder["mapping"] = source
                structure_holder["result"] = {
                    "source": source,
                    "mapping": structure_holder["mapping"],
                }
                _write_json(structure_path_json, structure_holder["result"])

            engine.run_stage(
                "structure",
                structure_action,
                {
                    "reference_sequence_sha256": hashlib.sha256(
                        query_sequence.encode()
                    ).hexdigest(),
                    **source_key,
                },
                [structure_path_json],
                dependencies=["alignment_qc"],
                parameters={"chain": structure_cfg.get("chain")},
                software_version=tool_version,
            )
            if not structure_holder:
                structure_holder["result"] = json.loads(
                    structure_path_json.read_text(encoding="utf-8")
                )
                engine.state["stages"]["structure"]["cache_hit"] = True
                engine._save()
            structure_result = structure_holder["result"]
            if structure_result.get("mapping", {}).get("status") != "mapping_successful":
                warnings.append(
                    "Structure mapping status: {}".format(
                        structure_result.get("mapping", {}).get(
                            "reason", structure_result.get("mapping", {}).get("status")
                        )
                    )
                )
        except Exception as exc:
            structure_result = {"status": "mapping_failed", "reason": str(exc)}
            _write_json(structure_path_json, structure_result)
            if engine.state["stages"].get("structure", {}).get("status") != "FAILED":
                engine.state["stages"]["structure"] = {
                    "name": "structure",
                    "status": "FAILED",
                    "error": str(exc),
                }
                engine._save()
            warnings.append("Structure mapping failed: {}".format(exc))
    else:
        engine.skip_stage("structure", "Disabled in configuration.", dependencies=["alignment_qc"])
        _write_json(structure_path_json, structure_result)

    scores = evolutionary_scores(col, selection_hits)
    stats.update({"alignment_length": len(aligned[0].sequence), "sequence_count": len(aligned)})
    metrics_path = out / "alignment_metrics.csv"

    def write_metrics():
        fields = [
            "position",
            "dominant_residue",
            "dominant_frequency",
            "entropy",
            "normalized_entropy",
            "information_content",
            "effective_alphabet_size",
            "gap_fraction",
            "ambiguous_fraction",
            "occupancy",
            "conservation_score",
            "residue_frequencies",
        ]
        with metrics_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in col:
                row_copy = dict(row)
                row_copy["residue_frequencies"] = json.dumps(
                    row_copy["residue_frequencies"], sort_keys=True
                )
                writer.writerow(row_copy)
        _write_json(out / "alignment_stats.json", stats)
        _write_json(out / "identity_matrix.json", matrix)
        _write_json(out / "evolutionary_trace_scores.json", scores)
        _write_json(out / "motif_hits.json", motif_hits)
        _write_json(out / "pairwise_dnds.json", selection_results)
        _write_json(out / "selection.json", selection_results)

    analysis_outputs = [
        metrics_path,
        out / "alignment_stats.json",
        out / "identity_matrix.json",
        out / "evolutionary_trace_scores.json",
        out / "motif_hits.json",
        out / "pairwise_dnds.json",
        out / "selection.json",
    ]
    engine.run_stage(
        "conservation_and_evidence",
        write_metrics,
        {"alignment_sha256": _sha(aligned_path)},
        analysis_outputs,
        dependencies=["alignment_qc"],
        parameters={"gap_policy": cfg.get("gap_policy"), "reference": cfg.get("reference")},
    )
    modality_statuses = {
        "selection": "not_configured"
        if not cfg.get("selection")
        else "available"
        if any(x.get("result_level") for x in selection_results)
        else "unavailable",
        "domains": "not_configured"
        if not domain_cfg.get("enabled")
        else "available"
        if engine.state.get("stages", {}).get("domains", {}).get("status") == "SUCCESS"
        else "unavailable",
        "motifs": "available" if motif_config else "not_configured",
        "structure": "not_configured"
        if not structure_cfg.get("enabled")
        else structure_result.get("mapping", structure_result).get("status", "unavailable"),
    }
    candidate_rows = candidate_sites(
        aligned,
        col,
        scores,
        cfg.get("reference"),
        motif_hits,
        domain_hits,
        selection_hits,
        structure_result.get("mapping", structure_result),
        modality_statuses=modality_statuses,
    )
    write_candidate_sites(candidate_rows, out)
    selection_levels = {}
    for result in selection_results:
        level = result.get("result_level")
        if level:
            selection_levels[level] = selection_levels.get(level, 0) + 1
    evidence_integration = {
        "status": "complete",
        "candidate_site_count": len(candidate_rows),
        "selection_result_levels": selection_levels,
        "evidence_families": {
            "alignment": {
                "status": "available",
                "measures": ["conservation", "entropy", "occupancy"],
                "aggregation": "one correlated descriptive family",
            },
            "selection": {
                "status": modality_statuses["selection"],
                "methods": sorted({str(x.get("method", "")) for x in selection_hits}),
                "aggregation": "methods retained individually and grouped into one dependent family",
                "result_levels": ["site", "branch", "gene", "gene_branch_set"],
            },
            "annotations": {
                "domains": modality_statuses["domains"],
                "motifs": modality_statuses["motifs"],
                "structure": modality_statuses["structure"],
                "aggregation": "context only; not added as independent adaptation evidence",
            },
        },
        "ets": {
            "status": "descriptive_alignment_proxy",
            "definition": "normalized_entropy * occupancy",
        },
    }
    _write_json(out / "evidence_integration.json", evidence_integration)
    profile_svg(col, out / "conservation.svg")
    if not (out / "phylogeny.svg").exists() and tree_path and str(tree_path).endswith(".treefile"):
        try:
            _draw_tree_svg(Phylo.read(str(tree_path), "newick"), out / "phylogeny.svg")
        except Exception as exc:
            warnings.append("Tree SVG export failed: {}".format(exc))
    _write_json(out / "phylogeny_summary.json", phylo_info)
    warning_append = []
    if alphabet == "protein":
        warning_append.append("Protein input: codon analyses were not run.")
    if cfg.get("selection") and not selection_hits:
        warning_append.append(
            "No inferential selection sites are available; absence of parsed hits is not evidence for neutrality."
        )
    if domain_cfg.get("enabled") and not domain_hits:
        warning_append.append(
            "Domain annotations are unavailable or yielded no hits; inspect domain stage status and database settings."
        )
    warnings.extend(warning_append)
    provenance = {
        "evotrace_version": __version__,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "input_path": str(inp),
        "input_sha256": input_hash,
        "input_accession_retrieval": retrieval_info,
        "configuration_path": str(Path(config_path).resolve()) if config_path else None,
        "configuration": cfg,
        "configuration_sha256": hashlib.sha256(
            json.dumps(cfg, sort_keys=True, default=str).encode()
        ).hexdigest(),
        "python_version": sys.version,
        "platform": platform.platform(),
        "cpu_count": __import__("os").cpu_count(),
        "software_versions": {
            "mafft": mafft_version,
            "biopython": __import__("Bio").__version__,
            "iqtree": _tool_version("iqtree2")
            if shutil.which("iqtree2")
            else _tool_version("iqtree"),
            "hyphy": _tool_version("hyphy"),
            "paml_codeml": _tool_version("codeml"),
            "hmmer": _tool_version("hmmscan"),
        },
        "commands": {
            name: stage.get("result") or stage.get("command")
            for name, stage in engine.state["stages"].items()
        },
        "stage_status": engine.summarize(),
        "runtime_seconds": round(time.time() - started, 3),
    }
    metadata_path = out / "run_metadata.json"
    _write_json(metadata_path, provenance)
    _write_json(out / "provenance.json", provenance)
    _write_json(out / "residue_mapping.json", AlignmentMap(aligned, cfg.get("reference")).rows())
    report_path = out / "report.html"
    render_report(
        stats,
        warnings,
        provenance,
        motif_hits,
        selection_results,
        report_path,
        candidate_rows=candidate_rows,
    )
    summary = {
        "status": "complete_with_skips"
        if any(s.get("status") == "SKIPPED" for s in engine.state["stages"].values())
        else "complete",
        "warnings": warnings,
        "outputs": sorted(p.name for p in out.iterdir()),
        "stage_status": engine.summarize(),
    }
    _write_json(out / "run_summary.json", summary)
    log.info("Analysis complete: %s", out)
    return summary
