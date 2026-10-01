"""Transparent descriptive constraint/divergence summaries."""

from typing import Any, Dict, List


def evolutionary_scores(
    columns: List[Dict[str, Any]], selection_sites: List[Dict[str, Any]] = None
) -> List[Dict[str, Any]]:
    """Return alignment trace summaries and a non-inferential site review proxy.

    Conservation and entropy are kept together as one correlated alignment family.
    Selection methods and annotations are joined in the site evidence ledger, not
    numerically summed: their scales and dependence structures are not calibrated.
    """
    methods_by_position = {}
    for site in selection_sites or []:
        position = site.get("alignment_position", site.get("site", site.get("position")))
        if isinstance(position, int) and position > 0 and site.get("result_level", "site") == "site":
            methods_by_position.setdefault(position, set()).add(
                str(site.get("method", site.get("test", "selection")))
            )
    scores = []
    for row in columns:
        # Constraint and divergence are complementary descriptions of the same column.
        constraint = float(row["conservation_score"]) * float(row["occupancy"])
        divergence = float(row["normalized_entropy"]) * float(row["occupancy"])
        scores.append(
            {
                "position": row["position"],
                "constraint_score": constraint,
                "divergence_score": divergence,
                "ets_descriptive_score": divergence,
                "ets_status": "descriptive_alignment_proxy",
                "selection_methods": sorted(methods_by_position.get(int(row["position"]), set())),
                "selection_family_count": int(bool(methods_by_position.get(int(row["position"])))),
                "interpretation": "ETS is a descriptive entropy/occupancy proxy. Selection methods are one dependent family; domain, motif, and structure are annotations.",
            }
        )
    return scores
