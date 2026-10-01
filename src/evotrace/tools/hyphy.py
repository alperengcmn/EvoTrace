"""HyPhy runners and method-specific parsers for documented JSON results."""

import json
import math
import re
from pathlib import Path
from typing import Dict, Optional

from evotrace.tools.common import locate, run

METHODS = {"FEL", "MEME", "FUBAR", "aBSREL", "BUSTED"}


def _method_name(method):
    return {name.lower(): name for name in METHODS}.get(str(method).lower())


def detect_hyphy():
    exe = locate(("hyphy",), "Install HyPhy from https://hyphy.org/.")
    res = run(exe, [exe, "--version"], timeout=30)
    return exe, ((res.stdout + res.stderr).strip().splitlines() or ["unknown"])[0]


def _norm(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def _number(value):
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def _pick(row, patterns):
    for key, value in row.items():
        if any(re.search(pattern, _norm(key)) for pattern in patterns):
            number = _number(value)
            if number is not None:
                return number
    return None


def benjamini_hochberg(p_values):
    """BH-adjust valid p-values in order; missing values remain missing."""
    ranked = sorted((float(p), i) for i, p in enumerate(p_values) if p is not None)
    result = [None] * len(p_values)
    minimum = 1.0
    for j in range(len(ranked) - 1, -1, -1):
        p, original = ranked[j]
        minimum = min(minimum, p * len(ranked) / (j + 1))
        result[original] = min(1.0, minimum)
    return result


def _mle_rows(data):
    """Read HyPhy FEL/MEME/FUBAR partitioned MLE headers/content tables."""
    mle = data.get("MLE", {})
    if not isinstance(mle, dict):
        return []
    headers, content = mle.get("headers", []), mle.get("content", {})
    # HyPhy serializes table headers as [name, description] pairs.
    headers = [
        item[0] if isinstance(item, (list, tuple)) and item else item for item in headers
    ]
    if isinstance(content, list):
        content = {str(i): part for i, part in enumerate(content)}
    partitions = data.get("data partitions", {})
    rows = []
    if not isinstance(content, dict):
        return rows
    for partition, values in content.items():
        if not isinstance(values, list):
            continue
        part_info = partitions.get(str(partition), {}) if isinstance(partitions, dict) else {}
        coverage = part_info.get("coverage", []) if isinstance(part_info, dict) else []
        if coverage and isinstance(coverage[0], (list, tuple)):
            coverage = coverage[0]
        coverage_offset = 1 if coverage and min(coverage) == 0 else 0
        for i, value in enumerate(values):
            row = (
                value
                if isinstance(value, dict)
                else dict(zip(headers, value))
                if isinstance(value, list) and len(value) == len(headers)
                else None
            )
            if row is None:
                continue
            site = next(
                (
                    _number(v)
                    for k, v in row.items()
                    if _norm(k) in {"site", "codon", "siteindex", "codonindex"}
                ),
                None,
            )
            if site is None:
                site = (
                    _number(coverage[i]) + coverage_offset
                    if i < len(coverage) and _number(coverage[i]) is not None
                    else i + 1
                )
            rows.append((str(partition), int(site), row))
    return rows


def _site_parser(data, method):
    sites = []
    for partition, site, row in _mle_rows(data):
        alpha = _pick(row, [r"^alpha$", r"^ds$"])
        beta = _pick(row, [r"^beta$", r"^dn$", r"betaplus"])
        p = _pick(row, [r"pvalue", r"^p$"])
        q = _pick(row, [r"qvalue", r"adjustedpvalue"])
        item = {
            "level": "site",
            "partition": partition,
            "site": site,
            "result_level": "site",
            "method": method,
            "alpha": alpha,
            "beta": beta,
            "dN": beta,
            "dS": alpha,
            "omega": beta / alpha if beta is not None and alpha not in (None, 0) else None,
            "lrt": _pick(row, [r"^lrt$"]),
            "p_value": p,
            "q_value": q,
            "statistics": row,
        }
        if method == "FEL":
            item["selection_class"] = (
                "positive/diversifying"
                if beta is not None and alpha is not None and beta > alpha
                else "negative/purifying"
                if beta is not None and alpha is not None and beta < alpha
                else "neutral_or_unclassified"
            )
        elif method == "MEME":
            item["selection_class"] = (
                "episodic_positive_candidate"
                if beta is not None and alpha is not None and beta > alpha
                else "episodic_positive_not_detected"
            )
        else:
            item["posterior_positive"] = _pick(
                row, [r"prob.*omega.*1", r"prob.*alpha.*beta", r"posterior.*positive"]
            )
            item["posterior_negative"] = _pick(
                row, [r"prob.*omega.*lt1", r"prob.*beta.*alpha", r"posterior.*negative"]
            )
            item["p_value"] = None
            item["q_value"] = None
            item["evidence_semantics"] = "posterior probability; not a p-value"
        sites.append(item)
    if method == "FUBAR":
        grid = data.get("grid", [])
        posterior = data.get("posterior", {})
        partitions = data.get("data partitions", {})
        if isinstance(posterior, dict):
            for partition, per_site in posterior.items():
                part = partitions.get(str(partition), {}) if isinstance(partitions, dict) else {}
                coverage = part.get("coverage", []) if isinstance(part, dict) else []
                if coverage and isinstance(coverage[0], (list, tuple)):
                    coverage = coverage[0]
                coverage_offset = 1 if coverage and min(coverage) == 0 else 0
                if not isinstance(per_site, list):
                    continue
                for i, probabilities in enumerate(per_site):
                    if not isinstance(probabilities, list) or len(probabilities) != len(grid):
                        continue
                    positive = 0.0
                    for cell, probability in zip(grid, probabilities):
                        if isinstance(cell, (list, tuple)) and len(cell) >= 2:
                            dn, ds, prob = _number(cell[0]), _number(cell[1]), _number(probability)
                            if (
                                dn is not None
                                and ds not in (None, 0)
                                and dn / ds > 1
                                and prob is not None
                            ):
                                positive += prob
                    site_number = (
                        int(_number(coverage[i]) + coverage_offset)
                        if i < len(coverage) and _number(coverage[i]) is not None
                        else i + 1
                    )
                    existing = next(
                        (
                            x
                            for x in sites
                            if str(x["partition"]) == str(partition) and x["site"] == site_number
                        ),
                        None,
                    )
                    if existing is None:
                        existing = {
                            "level": "site",
                            "result_level": "site",
                            "method": method,
                            "partition": str(partition),
                            "site": site_number,
                            "alpha": None,
                            "beta": None,
                            "dN": None,
                            "dS": None,
                            "omega": None,
                            "lrt": None,
                            "p_value": None,
                            "q_value": None,
                            "statistics": {},
                            "evidence_semantics": "posterior probability; not a p-value",
                        }
                        sites.append(existing)
                    existing["posterior_positive"] = positive
                    existing["selection_class"] = (
                        "positive_posterior" if positive >= 0.9 else "unclassified"
                    )
    if method in {"FEL", "MEME"}:
        for item, q in zip(sites, benjamini_hochberg([x["p_value"] for x in sites])):
            item["q_value"] = q
            if q is not None and q <= 0.05:
                item["selection_class"] += "_BH_q<=0.05"
    return sites


def _branch_parser(data):
    branches = []
    attrs = data.get("branch attributes", {})
    if not isinstance(attrs, dict):
        return branches
    for partition, nodes in attrs.items():
        if partition == "attributes" or not isinstance(nodes, dict):
            continue
        for name, values in nodes.items():
            if not isinstance(values, dict):
                continue
            corrected = _number(
                next((v for k, v in values.items() if _norm(k) == "correctedpvalue"), None)
            )
            raw_p = _number(
                next((v for k, v in values.items() if _norm(k) == "uncorrectedpvalue"), None)
            )
            branches.append(
                {
                    "level": "branch",
                    "result_level": "branch",
                    "partition": partition,
                    "branch": values.get("original name", name),
                    "rate_classes": values.get("Rate classes"),
                    "rate_distributions": values.get("Rate Distributions"),
                    "baseline_omega": values.get("Baseline MG94xREV omega"),
                    "lrt": _pick(values, [r"^lrt$"]),
                    "p_value": raw_p,
                    "adjusted_p_value": corrected,
                    "statistics": values,
                    "selection_class": "branch_positive_selection_candidate"
                    if corrected is not None and corrected <= 0.05
                    else "not_supported",
                }
            )
    return branches


def _busted_parser(data):
    tests, fits = data.get("test results", {}), data.get("fits", {})
    return [
        {
            "level": "gene_branch_set",
            "result_level": "gene_branch_set",
            "method": "BUSTED",
            "background": data.get("background"),
            "tested_branch_sets": data.get("tested"),
            "lrt": _pick(tests, [r"^lrt$"]),
            "p_value": _pick(tests, [r"pvalue"]),
            "log_likelihoods": {
                k: v.get("Log Likelihood") for k, v in fits.items() if isinstance(v, dict)
            },
            "rate_distributions": {
                k: v.get("Rate Distributions") for k, v in fits.items() if isinstance(v, dict)
            },
            "statistics": tests,
            "interpretation": "Gene/branch-set test; it does not identify selected sites.",
        }
    ]


def parse_hyphy(path: Path, method: str) -> Dict[str, object]:
    method = _method_name(method)
    if method is None:
        raise ValueError("Unsupported HyPhy method: {}".format(method))
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("HyPhy JSON output must be an object.")
    if method in {"FEL", "MEME", "FUBAR"}:
        return {
            "method": method,
            "result_level": "site",
            "sites": _site_parser(data, method),
            "raw": data,
        }
    if method == "aBSREL":
        return {
            "method": method,
            "result_level": "branch",
            "branches": _branch_parser(data),
            "raw": data,
        }
    return {
        "method": method,
        "result_level": "gene_branch_set",
        "gene_tests": _busted_parser(data),
        "raw": data,
    }


def run_hyphy(
    method: str, alignment: Path, output_dir: Path, tree: Optional[Path] = None, timeout=86400
) -> Dict[str, object]:
    method = _method_name(method)
    if method is None:
        raise ValueError(
            "Unsupported HyPhy analysis {!r}; choose {}.".format(method, ", ".join(sorted(METHODS)))
        )
    alignment, output_dir = Path(alignment).resolve(), Path(output_dir).resolve()
    if not alignment.is_file():
        raise ValueError("Alignment does not exist: {}".format(alignment))
    output_dir.mkdir(parents=True, exist_ok=True)
    exe, version = detect_hyphy()
    out = output_dir / (method.lower() + ".json")
    cmd = [exe, method.lower(), "--alignment", str(alignment), "--output", str(out)]
    if tree:
        if not Path(tree).is_file():
            raise ValueError("Tree does not exist: {}".format(tree))
        cmd += ["--tree", str(Path(tree).resolve())]
    result = run(exe, cmd, cwd=str(output_dir), timeout=timeout)
    if not out.is_file():
        raise RuntimeError("HyPhy completed without producing expected JSON: {}".format(out))
    parsed = parse_hyphy(out, method)
    parsed.update(
        {
            "path": str(out),
            "version": version,
            "command": cmd,
            "working_directory": result.working_directory,
            "exit_code": result.returncode,
            "runtime_seconds": result.runtime_seconds,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
    )
    return parsed
