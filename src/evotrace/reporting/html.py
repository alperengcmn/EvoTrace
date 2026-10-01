"""Self-contained HTML report renderer."""

import html
import json
from pathlib import Path


def _selection_status(result):
    if result.get("status"):
        return str(result["status"])
    for key, label in (
        ("sites", "site"),
        ("branches", "branch"),
        ("gene_tests", "gene/branch-set"),
    ):
        values = result.get(key)
        if isinstance(values, list):
            return "{} {} result(s)".format(len(values), label)
    if "p_value" in result:
        return "gene-level p={}".format(result["p_value"])
    return "No inferential result fields parsed"


def _selection_lrt(result):
    if result.get("lrt_status") == "computed":
        return "LRT {} (df {}, p={})".format(
            result.get("lrt", "—"), result.get("lrt_degrees_of_freedom", "—"),
            result.get("lrt_p_value", "—"),
        )
    return "LRT {}".format(result.get("lrt_status", "—"))


def render_report(stats, warnings, tools, motif_hits, selection, out: Path, candidate_rows=None):
    esc = html.escape
    warning_html = "".join("<li>{}</li>".format(esc(str(w))) for w in warnings) or "<li>None</li>"
    motif_html = (
        "".join(
            "<tr><td>{}</td><td>{}</td><td>{}-{} </td><td>{}</td><td><code>{}</code></td></tr>".format(
                esc(str(h["motif"])),
                esc(str(h["sequence_id"])),
                h["start"],
                h["end"],
                esc(str(h.get("alignment_position", "—"))),
                esc(str(h["matched"])),
            )
            for h in motif_hits
        )
        or '<tr><td colspan="5">No user-defined motif hits</td></tr>'
    )
    sel_html = (
        "".join(
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
                esc(str(x.get("method", "pairwise NG" if x.get("sequence_a") else "unknown"))),
                esc(str(x.get("sequence_a", x.get("result_level", x.get("status", "—"))))),
                esc(
                    str(
                        x.get(
                            "sequence_b",
                            x.get(
                                "status",
                                _selection_status(x),
                            ),
                        )
                    )
                ),
                esc(
                    str(
                        x.get(
                            "dN",
                            x.get("p_value", (x.get("gene_tests") or [{}])[0].get("p_value", "—")),
                        )
                    )
                ),
                esc(str(x.get("dS", x.get("dN_dS", x.get("omega", "—"))))),
                esc(_selection_lrt(x)),
            )
            for x in selection
        )
        or '<tr><td colspan="6">No codon selection analysis configured or executed.</td></tr>'
    )
    candidate_rows = candidate_rows or []
    ranked = sorted(
        candidate_rows, key=lambda x: float(x.get("divergence_score") or 0), reverse=True
    )[:20]
    candidate_html = "".join(
        "<tr><td>{}</td><td>{}</td><td>{}</td><td>{:.3f}</td><td>{:.3f}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
            esc(str(x.get("alignment_position"))),
            esc(str(x.get("reference_position"))),
            esc(str(x.get("reference_residue"))),
            float(x.get("entropy") or 0),
            float(x.get("conservation") or 0),
            esc(str(x.get("domain") or "—")),
            esc(str(x.get("motif") or "—")),
            esc(
                "; ".join(
                    e.get("description", e.get("evidence_type", ""))
                    for e in json.loads(x.get("evidence_summary", "[]"))
                )
            ),
        )
        for x in ranked
    )
    body = """<!doctype html><html><head><meta charset="utf-8"><title>EvoTrace report</title><style>body{font:16px/1.55 system-ui;max-width:1100px;margin:36px auto;padding:0 24px;color:#18313b}h1,h2{color:#126777}table{border-collapse:collapse;width:100%}td,th{border:1px solid #d4e0e3;padding:8px;text-align:left}code{background:#edf4f5;padding:2px 5px}.note{background:#fff7df;padding:12px;border-left:4px solid #dba62f}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px}.card{padding:12px;border:1px solid #d4e0e3;border-radius:8px}</style></head><body><h1>EvoTrace evolutionary analysis</h1><p>Generated from computed input data. Interpret results within the documented methods and limitations.</p><div class="grid">"""
    for k in (
        "sequence_count",
        "alignment_length",
        "gap_fraction",
        "conserved_columns",
        "variable_columns",
        "mean_pairwise_identity",
    ):
        val = stats.get(k)
        body += '<div class="card"><b>{}</b><br>{}</div>'.format(
            esc(k.replace("_", " ").title()),
            esc(
                "NA"
                if val is None
                else ("{:.3f}".format(val) if isinstance(val, float) else str(val))
            ),
        )
    body += (
        """</div><h2>Executive summary</h2><p>Analyzed {} {} sequences across {} aligned positions. {} columns are invariant among observed canonical residues; mean pairwise identity is {}. These summaries describe this input set and do not establish adaptation.</p><h2>Dataset</h2><p>Input: {}. Alphabet: {}. Lengths: {}</p><h2>Conservation profile</h2><img src="conservation.svg" alt="Conservation profile" style="width:100%"><h2>Phylogeny</h2><p>{}</p><p>{}</p><p>See <code>phylogeny.svg</code> and <code>phylogeny.nwk</code> when generated.</p><h2>Candidate positions for review</h2><p>Sorted by descriptive divergence score; these are candidates for further study, not adaptive calls.</p><table><tr><th>Alignment pos.</th><th>Reference pos.</th><th>Reference residue</th><th>Entropy</th><th>Conservation</th><th>Domains</th><th>Motifs</th><th>Evidence</th></tr>""".format(
            esc(str(stats.get("sequence_count", 0))),
            esc(str(stats.get("alphabet", "unknown"))),
            esc(str(stats.get("alignment_length", 0))),
            esc(str(stats.get("conserved_columns", 0))),
            esc(
                "NA"
                if stats.get("mean_pairwise_identity") is None
                else "{:.3f}".format(stats.get("mean_pairwise_identity"))
            ),
            esc(Path(str(tools.get("input_path", "unknown"))).name),
            esc(str(stats.get("alphabet", "unknown"))),
            esc(json.dumps(stats.get("sequence_lengths", {}), sort_keys=True)),
            esc(str(stats.get("phylogeny_method", "not run"))),
            esc(str(stats.get("phylogeny_support", ""))),
        )
        + candidate_html
        + """</table><h2>Warnings and scope</h2><ul>"""
        + warning_html
        + """</ul><div class="note"><b>Scientific interpretation:</b> conservation, entropy and pairwise dN/dS are descriptive. This report does not infer positive selection or prove adaptation. Domain/motif and structure information is included only when actually computed; absent annotations are reported as unavailable.</div><h2>User-defined motif hits</h2><table><tr><th>Motif</th><th>Sequence</th><th>Sequence positions</th><th>Alignment position</th><th>Match</th></tr>"""
        + motif_html
        + """</table><h2>Selection analyses</h2><p>Pairwise descriptive dN/dS and model-based selection tests are shown as separate result records. P-values and posterior probabilities retain their own meanings.</p><table><tr><th>Method</th><th>Result level / pair</th><th>Status / partner</th><th>dN / p-value</th><th>dS / omega</th><th>Model LRT</th></tr>"""
        + sel_html
        + """</table><h2>Reproducibility</h2><p>Tool details and checksums are in <code>run_metadata.json</code>; machine-readable exports include CSV, TSV and JSON.</p><pre>"""
        + esc(json.dumps(tools, indent=2))
        + """</pre></body></html>"""
    )
    stage_status = tools.get("stage_status", {})
    codon_path = out.parent / "codon_alignment_qc.json"
    codon_qc = json.loads(codon_path.read_text(encoding="utf-8")) if codon_path.exists() else None
    structure_path = out.parent / "structure_mapping.json"
    structure_data = (
        json.loads(structure_path.read_text(encoding="utf-8")) if structure_path.exists() else {}
    )
    structure_status = structure_data.get("mapping", structure_data).get("status", "unavailable")
    evidence_path = out.parent / "evidence_integration.json"
    evidence_data = json.loads(evidence_path.read_text(encoding="utf-8")) if evidence_path.exists() else {}
    domain_count = sum(bool(row.get("domain")) for row in candidate_rows)
    sections = [
        (
            "Input dataset",
            "<p>Input file: <code>{}</code>; SHA-256: <code>{}</code>.</p>".format(
                esc(Path(str(tools.get("input_path", "unknown"))).name),
                esc(str(tools.get("input_sha256", "unavailable"))),
            ),
        ),
        (
            "Validation",
            "<p>Stage status: <b>{}</b>.</p>".format(
                esc(str(stage_status.get("validation", "unavailable")))
            ),
        ),
        (
            "Alignment QC",
            "<p>Alignment QC: <b>{}</b>; {} sequences, {} columns.</p>".format(
                esc(str(stage_status.get("alignment_qc", "unavailable"))),
                esc(str(stats.get("sequence_count"))),
                esc(str(stats.get("alignment_length"))),
            ),
        ),
        (
            "Codon QC",
            "<p>{}</p>".format(
                esc(
                    json.dumps(codon_qc, sort_keys=True)
                    if codon_qc
                    else "Not run: no codon alignment was supplied or generated."
                )
            ),
        ),
        (
            "Pairwise dN/dS",
            "<p>{}</p>".format(
                esc(
                    "{} descriptive pairwise estimate(s) in the result table.".format(
                        sum("sequence_a" in x for x in selection)
                    )
                    if any("sequence_a" in x for x in selection)
                    else "Not run: no validated codon alignment was available."
                )
            ),
        ),
        (
            "Selection inference",
            "<p>{}</p>".format(
                esc(
                    "Method-specific results are listed above; see machine-readable selection.json for full statistics."
                    if any(x.get("result_level") or x.get("method") for x in selection)
                    else "Not run or unavailable. This is not evidence for neutrality."
                )
            ),
        ),
        (
            "Conservation and entropy",
            "<p>Invariant columns: {}; variable columns: {}; mean pairwise identity: {}. Entropy is descriptive and sampling-dependent.</p>".format(
                esc(str(stats.get("conserved_columns"))),
                esc(str(stats.get("variable_columns"))),
                esc(str(stats.get("mean_pairwise_identity"))),
            ),
        ),
        (
            "Domains",
            "<p>Stage: {}; candidate positions overlapping at least one annotated domain: {}. No hits/unavailable does not establish biological absence.</p>".format(
                esc(str(stage_status.get("domains", "unavailable"))), esc(str(domain_count))
            ),
        ),
        (
            "Motifs",
            "<p>Observed user-defined motif hits: {}. Coordinates are mapped through the alignment when available.</p>".format(
                esc(str(len(motif_hits)))
            ),
        ),
        (
            "Structure",
            "<p>Structure mapping status: <b>{}</b>. Missing mapping is not structural evidence of absence.</p>".format(
                esc(str(structure_status))
            ),
        ),
        (
            "Evolutionary Trace Score",
            "<p>ETS is a descriptive normalized-entropy × occupancy proxy. The integration ledger groups conservation and entropy as one correlated family, retains selection methods and site/branch/gene result levels, and displays domain, motif, and structure as context. Status: <b>{}</b>.</p><pre>{}</pre>".format(
                esc(str(evidence_data.get("ets", {}).get("status", "unavailable"))),
                esc(json.dumps(evidence_data.get("evidence_families", {}), indent=2, sort_keys=True)),
            ),
        ),
        (
            "Candidate sites",
            "<p>{} alignment positions are tabulated; displayed candidates are ranked by descriptive divergence and include available evidence annotations.</p>".format(
                esc(str(len(candidate_rows)))
            ),
        ),
        (
            "Interpretation and limitations",
            "<p>Sites are candidates for further biological investigation. This analysis does not prove adaptation, function, or causality. Orthology is not independently validated.</p>",
        ),
        (
            "Reproducibility",
            "<p>Configuration SHA-256: <code>{}</code>. Per-stage inputs, outputs, parameters, tool versions, commands, and cache records are in <code>pipeline_stages.json</code>.</p>".format(
                esc(str(tools.get("configuration_sha256", "unavailable")))
            ),
        ),
    ]
    appendix = "".join("<h2>{}</h2>{}".format(esc(title), content) for title, content in sections)
    body = body.replace("</body>", appendix + "</body>")
    out.write_text(body, encoding="utf-8")
