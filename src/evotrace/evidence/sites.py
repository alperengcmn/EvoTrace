"""Per-alignment-position evidence ledger and portable candidate-site exports."""

import csv
import json
from pathlib import Path

from evotrace.evidence.models import Evidence
from evotrace.mapping.residues import AlignmentMap


def _position(value):
    try:
        number = int(value)
    except (TypeError, ValueError):
        return None
    return number if number > 0 else None


def _status(name, has_hit, statuses):
    if has_hit:
        return "hit"
    value = statuses.get(name, "not_run")
    if value in ("SUCCESS", "available", "mapping_successful"):
        return "no_hit" if name != "structure" else "unmapped"
    return value


def candidate_sites(
    alignment,
    conservation,
    scores,
    reference_id=None,
    motifs=None,
    domains=None,
    selection=None,
    structure_mapping=None,
    modality_statuses=None,
):
    """Join annotations on alignment columns without treating missing data as zero evidence.

    HyPhy site methods and PAML BEB are mapped through codon-site numbers to the
    matching protein alignment columns. Branch/gene-level tests remain intact in
    selection.json and are deliberately not projected onto residue rows.
    """
    mapping = AlignmentMap(alignment, reference_id)
    if len(conservation) != mapping.length or len(scores) != mapping.length:
        raise ValueError("Conservation metrics and scores must match alignment columns exactly.")
    motifs, domains, selection = motifs or [], domains or [], selection or []
    statuses = modality_statuses or {}
    motif_by_pos, domain_by_pos, selection_by_pos, structure_by_pos = {}, {}, {}, {}
    domain_coordinates = {}

    for hit in motifs:
        try:
            seq_id = hit["sequence_id"]
            start, end = _position(hit.get("start")), _position(hit.get("end"))
            if not start or not end or end < start:
                continue
            for seq_pos in range(start, end + 1):
                aln_pos = mapping.sequence_to_alignment(seq_id, seq_pos)
                motif_by_pos.setdefault(aln_pos, []).append(hit)
        except (KeyError, ValueError, TypeError):
            continue

    for hit in domains:
        try:
            seq_id = hit["sequence_id"]
            start, end = _position(hit.get("start")), _position(hit.get("end"))
            if not start or not end or end < start:
                continue
            ref_positions = []
            for seq_pos in range(start, end + 1):
                aln_pos = mapping.sequence_to_alignment(seq_id, seq_pos)
                domain_by_pos.setdefault(aln_pos, []).append(hit)
                ref_pos = mapping.alignment_to_reference(aln_pos)
                if ref_pos is not None:
                    ref_positions.append(ref_pos)
            domain_coordinates[id(hit)] = {
                "domain": hit.get("domain_name", hit.get("domain_id", "unknown")),
                "sequence_id": seq_id,
                "query_start": start,
                "query_end": end,
                "reference_start": min(ref_positions) if ref_positions else None,
                "reference_end": max(ref_positions) if ref_positions else None,
                "evalue": hit.get("evalue"),
                "score": hit.get("score"),
            }
        except (KeyError, ValueError, TypeError):
            continue

    for hit in selection:
        # HyPhy and codeml sites are 1-based codon positions; one codon maps to
        # one protein alignment column after the codon/protein correspondence QC.
        pos = _position(hit.get("alignment_position", hit.get("site", hit.get("position"))))
        if pos is not None and pos <= mapping.length and hit.get("result_level", "site") == "site":
            selection_by_pos.setdefault(pos, []).append(hit)

    for mapped in (structure_mapping or {}).get("mapped_residues", []):
        try:
            seq_pos = _position(mapped.get("sequence_position"))
            if seq_pos:
                aln_pos = mapping.sequence_to_alignment(mapping.reference_id, seq_pos)
                structure_by_pos[aln_pos] = mapped
        except (KeyError, ValueError, TypeError):
            continue

    ref = next(record for record in alignment if record.id == mapping.reference_id)
    rows = []
    for col, score in zip(conservation, scores):
        pos = int(col["position"])
        if pos < 1 or pos > mapping.length:
            raise ValueError("Conservation position lies outside the alignment: {}".format(pos))
        residue = ref.sequence[pos - 1]
        site_selection = selection_by_pos.get(pos, [])
        site_motifs = motif_by_pos.get(pos, [])
        site_domains = domain_by_pos.get(pos, [])
        structure_hit = structure_by_pos.get(pos)
        evidence = [
            Evidence(
                "alignment",
                pos,
                "conservation",
                "descriptive",
                "Conservation {:.3f}; entropy {:.3f} bits. These are correlated summaries of one alignment evidence family.".format(
                    col["conservation_score"], col["entropy"]
                ),
                evidence_family="alignment",
                result_level="site",
            ).to_dict()
        ]
        for test in site_selection:
            pv = test.get("p_value")
            qv = test.get("q_value", test.get("adjusted_p_value"))
            posterior = test.get("posterior_probability", test.get("posterior_positive"))
            evidence.append(
                Evidence(
                    str(test.get("method", test.get("test", "selection"))),
                    pos,
                    str(test.get("selection_class", "selection_candidate")),
                    "p_value_reported" if pv is not None else "posterior" if posterior is not None else "unclassified",
                    str(
                        test.get(
                            "evidence_status",
                            "Posterior probability {:.3f}; not a p-value.".format(float(posterior))
                            if posterior is not None
                            else "See the method-specific selection result and statistics.",
                        )
                    ),
                    float(pv) if pv is not None else None,
                    float(qv) if qv is not None else None,
                    str(test.get("method", test.get("test", "selection"))),
                    "selection",
                    "site",
                ).to_dict()
            )
        for hit in site_motifs:
            evidence.append(
                Evidence(
                    "motif",
                    pos,
                    "motif",
                    "annotation",
                    "Overlaps user-defined motif {} in {} positions {}-{}.".format(
                        hit.get("motif", "unknown"), hit.get("sequence_id"), hit.get("start"), hit.get("end")
                    ),
                    method="user_regex",
                    evidence_family="annotation",
                    result_level="site",
                ).to_dict()
            )
        for hit in site_domains:
            evidence.append(
                Evidence(
                    "HMMER",
                    pos,
                    "domain",
                    "annotation",
                    "Overlaps {} (E-value {}).".format(
                        hit.get("domain_name", hit.get("domain_id", "domain")), hit.get("evalue")
                    ),
                    method="hmmscan",
                    evidence_family="annotation",
                    result_level="site",
                ).to_dict()
            )
        if structure_hit:
            evidence.append(
                Evidence(
                    "structure",
                    pos,
                    "structure",
                    "mapped",
                    "Mapped to chain {} residue {}{}.".format(
                        structure_hit.get("structure_chain"),
                        structure_hit.get("structure_residue_number"),
                        structure_hit.get("structure_insertion_code", ""),
                    ),
                    method="sequence_structure_alignment",
                    evidence_family="annotation",
                    result_level="site",
                ).to_dict()
            )

        qvalues = [float(x["q_value"]) for x in site_selection if x.get("q_value") is not None]
        qvalues += [float(x["adjusted_p_value"]) for x in site_selection if x.get("adjusted_p_value") is not None]
        pvalues = [float(x["p_value"]) for x in site_selection if x.get("p_value") is not None]
        posteriors = [
            float(x[key])
            for x in site_selection
            for key in ("posterior_probability", "posterior_positive")
            if x.get(key) is not None
        ]
        coords = [domain_coordinates[id(hit)] for hit in site_domains if id(hit) in domain_coordinates]
        unique_motifs = sorted({str(x.get("motif", "unknown")) for x in site_motifs})
        unique_domains = sorted(
            {str(x.get("domain_name", x.get("domain_id", "unknown"))) for x in site_domains}
        )
        modality_rows = {
            "selection": _status("selection", bool(site_selection), statuses),
            "domains": _status("domains", bool(site_domains), statuses),
            "motifs": _status("motifs", bool(site_motifs), statuses),
            "structure": "mapped"
            if structure_hit
            else _status(
                "structure",
                False,
                {
                    **statuses,
                    "structure": statuses.get(
                        "structure", (structure_mapping or {}).get("status", "not_run")
                    ),
                },
            ),
        }
        selection_methods = sorted(
            {str(x.get("method", x.get("test", "selection"))) for x in site_selection}
        )
        evidence_families = {
            "alignment": {
                "status": "available",
                "measures": ["conservation", "entropy", "occupancy"],
                "counted_as_one_correlated_family": True,
            },
            "selection": {
                "status": modality_rows["selection"],
                "methods": selection_methods,
                "result_levels": sorted(
                    {str(x.get("result_level", "site")) for x in site_selection}
                ),
                "counted_as_one_family": True,
            },
            "annotations": {
                "domains": modality_rows["domains"],
                "motifs": modality_rows["motifs"],
                "structure": modality_rows["structure"],
                "numeric_adaptation_weight": None,
            },
        }
        rows.append(
            {
                "alignment_position": pos,
                "reference_position": mapping.alignment_to_reference(pos),
                "reference_residue": residue,
                "dominant_residue": col["dominant_residue"],
                "entropy": col["entropy"],
                "conservation": col["conservation_score"],
                "gap_fraction": col["gap_fraction"],
                "occupancy": col["occupancy"],
                "identity": col["dominant_frequency"],
                "selection_status": modality_rows["selection"],
                "selection_method": ";".join(sorted({str(x.get("method", x.get("test", ""))) for x in site_selection})),
                "selection_pvalue": min(pvalues) if pvalues else None,
                "selection_qvalue": min(qvalues) if qvalues else None,
                "selection_posterior": max(posteriors) if posteriors else None,
                "selection_result_levels": sorted({str(x.get("result_level", "site")) for x in site_selection}),
                "domain_status": _status("domains", bool(site_domains), statuses),
                "domain": ";".join(unique_domains),
                "domain_start": min([int(x["start"]) for x in site_domains], default=None),
                "domain_end": max([int(x["end"]) for x in site_domains], default=None),
                "domain_reference_start": min([x["reference_start"] for x in coords if x["reference_start"] is not None], default=None),
                "domain_reference_end": max([x["reference_end"] for x in coords if x["reference_end"] is not None], default=None),
                "domain_score": max([float(x["score"]) for x in site_domains if x.get("score") is not None], default=None),
                "domain_evalue": min([float(x["evalue"]) for x in site_domains if x.get("evalue") is not None], default=None),
                "domain_coordinates": json.dumps(coords, sort_keys=True),
                "motif_status": modality_rows["motifs"],
                "motif": ";".join(unique_motifs),
                "structure_status": modality_rows["structure"],
                "structure_id": (structure_mapping or {}).get("structure_id"),
                "structure_chain": structure_hit.get("structure_chain") if structure_hit else None,
                "structure_residue": "{}:{}{}".format(
                    structure_hit.get("structure_chain"),
                    structure_hit.get("structure_residue_number"),
                    structure_hit.get("structure_insertion_code", ""),
                ) if structure_hit else None,
                "structure_confidence": structure_hit.get("mean_b_factor")
                if structure_hit and (structure_mapping or {}).get("source_provider") == "AlphaFold DB"
                else None,
                "structure_b_factor": structure_hit.get("mean_b_factor") if structure_hit else None,
                "constraint_score": score["constraint_score"],
                "divergence_score": score["divergence_score"],
                "ets_descriptive_score": score.get("ets_descriptive_score"),
                "ets_status": score.get("ets_status", "descriptive_alignment_proxy"),
                "evidence_family_summary": json.dumps(evidence_families, sort_keys=True),
                "evidence_summary": json.dumps(evidence, sort_keys=True),
            }
        )
    return rows


def write_candidate_sites(rows, output_dir):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["alignment_position", "reference_position", "reference_residue"]
    for name, delimiter in (("candidate_sites.csv", ","), ("candidate_sites.tsv", "\t")):
        with (out / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(rows)
    (out / "candidate_sites.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
