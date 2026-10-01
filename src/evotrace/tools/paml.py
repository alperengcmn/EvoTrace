"""PAML/codeml model runners, nested likelihood-ratio tests, and labeled parsers."""

import math
import re
from pathlib import Path
from io import StringIO

from Bio import AlignIO, Phylo
from Bio.Data import CodonTable
from evotrace.tools.common import locate, run

MODELS = {"branch": (2, 0), "site": (0, 2), "branch_site": (2, 2)}


def parse_codeml_output(text, model):
    """Parse model likelihood, omega estimates, and BEB posterior site results."""
    likelihoods = [
        float(x)
        for x in re.findall(
            r"lnL\([^\n]*?\)\s*:\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)", text
        )
    ]
    omegas = [
        float(x)
        for x in re.findall(
            r"omega\s*\(\s*dN\s*/\s*dS\s*\)\s*[:=]\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
            text,
            re.I,
        )
    ]
    branch_omega_classes = []
    branch_rate_table = []
    class_match = re.search(r"w\s*\(dN/dS\) for branches:\s*([^\n]+)", text, re.I)
    if class_match:
        branch_omega_classes = [
            float(value)
            for value in re.findall(r"\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?", class_match.group(1))
        ]
    in_branch_table = False
    for line in text.splitlines():
        if "dN & dS for each branch" in line:
            in_branch_table = True
            continue
        if not in_branch_table:
            continue
        values = line.split()
        if len(values) >= 7 and re.fullmatch(r"\d+\.\.\d+", values[0]):
            try:
                branch_rate_table.append(
                    {
                        "paml_branch_label": values[0],
                        "omega": float(values[4]),
                        "dN": float(values[5]),
                        "dS": float(values[6]),
                    }
                )
            except ValueError:
                continue
        elif branch_rate_table:
            in_branch_table = False
    beb_sites = []
    if model in ("site", "branch_site"):
        in_beb = False
        for line in text.splitlines():
            if "Bayes Empirical Bayes" in line or "Naive Empirical Bayes" in line:
                in_beb = True
                continue
            if in_beb:
                match = re.match(r"\s*(\d+)\s+([A-Z*])\s+([0-9]*\.?[0-9]+)(\*{0,3})\s*(.*)$", line)
                if match:
                    beb_sites.append(
                        {
                            "level": "site",
                            "result_level": "site",
                            "method": "PAML/codeml",
                            "site": int(match.group(1)),
                            "amino_acid": match.group(2),
                            "posterior_probability": float(match.group(3)),
                            "star_annotation": match.group(4),
                            "rate_class": match.group(5).strip() or None,
                            "evidence_semantics": "Bayesian posterior probability; not a p-value",
                        }
                    )
                elif beb_sites and not line.strip():
                    in_beb = False
    return {
        "model": model,
        "result_level": "gene" if model in ("site", "branch_site") else "branch",
        "site_result_level": "site" if beb_sites else "unavailable",
        "log_likelihood": likelihoods[-1] if likelihoods else None,
        "likelihood_records": likelihoods,
        "omega_estimates": omegas,
        "branch_omega_classes": branch_omega_classes,
        "branch_rate_table": branch_rate_table,
        "site_posterior_results": beb_sites,
    }


def _control_text(alignment, tree, outfile, model, nssites, fix_omega, omega):
    return "\n".join(
        (
            "seqfile = {}".format(alignment),
            "treefile = {}".format(tree),
            "outfile = {}".format(outfile),
            "noisy = 0",
            "verbose = 1",
            "runmode = 0",
            "seqtype = 1",
            # EvoTrace validates the standard genetic code before codeml. Keep
            # the PAML control explicit instead of relying on its implicit default.
            "icode = 0",
            "CodonFreq = 2",
            "model = {}".format(model),
            "NSsites = {}".format(nssites),
            "fix_omega = {}".format(fix_omega),
            "omega = {}".format(omega),
            "cleandata = 0",
            "",
        )
    )


def _likelihood_ratio(null_lnl, alt_lnl, model):
    if null_lnl is None or alt_lnl is None:
        return {"statistic": None, "p_value": None, "status": "likelihood_unavailable"}
    statistic = max(0.0, 2.0 * (alt_lnl - null_lnl))
    if model == "site":
        df = 2  # M1a (NSsites=1) nested in M2a (NSsites=2).
        p_value = math.exp(-statistic / 2.0)
        comparison = {"null": "M1a", "alternative": "M2a", "distribution": "chi_square_df2"}
    elif model == "branch_site":
        df = 1
        p_value = 0.5 * math.erfc(math.sqrt(statistic / 2.0)) if statistic else 1.0
        comparison = {
            "null": "branch-site omega2 fixed at 1",
            "alternative": "branch-site omega2 free",
            "distribution": "50:50 point mass at zero and chi_square_df1",
        }
    else:
        df = 1
        p_value = math.erfc(math.sqrt(statistic / 2.0))
        comparison = {
            "null": "one-ratio model",
            "alternative": "two-ratio branch model",
            "distribution": "chi_square_df1",
        }
    return {
        "statistic": statistic,
        "degrees_of_freedom": df,
        "p_value": min(1.0, max(0.0, p_value)),
        "status": "computed",
        "comparison": comparison,
    }


def _prepare_codeml_inputs(alignment_path, tree_path, output_dir):
    """Convert supported aligned inputs to PAML PHYLIP and matching short tree labels."""
    alignment = None
    errors = []
    for fmt in ("fasta", "phylip-relaxed", "phylip"):
        try:
            alignment = AlignIO.read(str(alignment_path), fmt)
            break
        except Exception as exc:
            errors.append("{}: {}".format(fmt, exc))
    if alignment is None:
        raise ValueError(
            "Cannot parse codon alignment as FASTA or PHYLIP: {}".format("; ".join(errors))
        )
    if len(alignment) < 3 or len({len(record.seq) for record in alignment}) != 1:
        raise ValueError("PAML needs at least three equal-length aligned sequences.")
    if len(alignment[0].seq) % 3:
        raise ValueError("PAML codon alignment length must be divisible by three.")
    ids = [record.id for record in alignment]
    if len(ids) != len(set(ids)):
        raise ValueError("PAML sequence IDs must be unique.")
    normalized = []
    for record in alignment:
        sequence = str(record.seq).upper().replace("U", "T").replace(".", "-")
        for i in range(0, len(sequence), 3):
            codon = sequence[i : i + 3]
            if "-" in codon and codon != "---":
                raise ValueError("{} has a partial codon gap at alignment position {}.".format(record.id, i // 3 + 1))
            if codon != "---" and set(codon) - set("ACGT"):
                raise ValueError("{} has an ambiguous codon at alignment position {}.".format(record.id, i // 3 + 1))
            if codon != "---" and codon in CodonTable.unambiguous_dna_by_name["Standard"].stop_codons:
                raise ValueError("{} has a stop codon at alignment position {}.".format(record.id, i // 3 + 1))
        normalized.append((record.id, sequence))

    tree_text = Path(tree_path).read_text(encoding="utf-8")
    lines = [line.strip() for line in tree_text.splitlines() if line.strip() and not line.lstrip().startswith("//")]
    if lines and re.fullmatch(r"\d+\s+\d+", lines[0]):
        lines = lines[1:]
    newick = "\n".join(lines)
    try:
        tree = Phylo.read(StringIO(newick), "newick")
    except Exception as exc:
        raise ValueError("Cannot parse selection tree as one Newick tree: {}".format(exc)) from exc
    tips = [terminal.name for terminal in tree.get_terminals()]
    if len(tips) != len(set(tips)) or set(tips) != set(ids):
        raise ValueError("PAML tree tip names must match codon alignment IDs exactly.")
    aliases = {seq_id: "ET{:06d}".format(i + 1) for i, (seq_id, _) in enumerate(normalized)}
    for terminal in tree.get_terminals():
        terminal.name = aliases[terminal.name]

    phylip_path = Path(output_dir) / "codeml_alignment.phy"
    tree_path_out = Path(output_dir) / "codeml_tree.trees"
    with phylip_path.open("w", encoding="utf-8") as handle:
        handle.write("{} {}\n".format(len(normalized), len(normalized[0][1])))
        for seq_id, sequence in normalized:
            handle.write("{:<10}{}\n".format(aliases[seq_id], sequence))
    newick_out = StringIO()
    Phylo.write(tree, newick_out, "newick")
    tree_path_out.write_text(
        "{} 1\n{}".format(len(normalized), newick_out.getvalue()), encoding="utf-8"
    )
    return phylip_path, tree_path_out, aliases


def run_codeml(
    alignment: Path, tree: Path, output_dir: Path, model="branch_site", omega=1.0, timeout=86400,
    genetic_code=1,
):
    if model not in MODELS:
        raise ValueError("model must be branch, site, or branch_site")
    alignment, tree, output_dir = (
        Path(alignment).resolve(),
        Path(tree).resolve(),
        Path(output_dir).resolve(),
    )
    if not alignment.is_file() or not tree.is_file():
        raise ValueError("codeml requires existing codon alignment and tree files.")
    if omega <= 0:
        raise ValueError("omega must be positive")
    if int(genetic_code) != 1:
        raise ValueError("PAML codeml currently supports only genetic_code=1 (standard code).")
    output_dir.mkdir(parents=True, exist_ok=True)
    paml_alignment, paml_tree, aliases = _prepare_codeml_inputs(alignment, tree, output_dir)
    tree_text = paml_tree.read_text(encoding="utf-8")
    if model in ("branch", "branch_site") and not re.search(r"#1(?!\d)", tree_text):
        raise ValueError(
            "PAML {} requires a user-specified foreground branch tagged #1 in the Newick tree.".format(model)
        )
    exe = locate(("codeml",), "Install PAML and ensure codeml is on PATH.")
    alt_model, alt_sites = MODELS[model]
    if model == "site":
        null_model, null_sites, null_fix = 0, 1, 0
    elif model == "branch_site":
        null_model, null_sites, null_fix = alt_model, alt_sites, 1
    else:
        null_model, null_sites, null_fix = 0, 0, 0

    alt_out, null_out = output_dir / "codeml.out", output_dir / "codeml_null.out"
    alt_ctl, null_ctl = output_dir / "codeml.ctl", output_dir / "codeml_null.ctl"
    alt_ctl.write_text(
        _control_text(paml_alignment, paml_tree, alt_out, alt_model, alt_sites, 0, omega), encoding="utf-8"
    )
    null_ctl.write_text(
        _control_text(paml_alignment, paml_tree, null_out, null_model, null_sites, null_fix, omega),
        encoding="utf-8",
    )
    alternative_run = run(exe, [exe, str(alt_ctl)], cwd=str(output_dir), timeout=timeout)
    if not alt_out.is_file():
        raise RuntimeError("codeml alternative model returned without creating its output file.")
    alternative_text = alt_out.read_text(errors="replace")
    alternative = parse_codeml_output(alternative_text, model)
    null_run = run(exe, [exe, str(null_ctl)], cwd=str(output_dir), timeout=timeout)
    if not null_out.is_file():
        raise RuntimeError("codeml null model returned without creating its output file.")
    null_text = null_out.read_text(errors="replace")
    null = parse_codeml_output(null_text, model)
    lrt = _likelihood_ratio(null["log_likelihood"], alternative["log_likelihood"], model)
    version_match = re.search(r"paml version\s+([0-9A-Za-z.]+)", alternative_text, re.I)
    return {
        **alternative,
        "control_file": str(alt_ctl),
        "output_file": str(alt_out),
        "null_control_file": str(null_ctl),
        "null_output_file": str(null_out),
        "paml_alignment": str(paml_alignment),
        "paml_tree": str(paml_tree),
        "taxon_aliases": aliases,
        "null_log_likelihood": null["log_likelihood"],
        "lrt": lrt.get("statistic"),
        "lrt_p_value": lrt.get("p_value"),
        "lrt_degrees_of_freedom": lrt.get("degrees_of_freedom"),
        "lrt_status": lrt["status"],
        "lrt_comparison": lrt.get("comparison"),
        "version": "PAML {}".format(version_match.group(1)) if version_match else "unknown",
        "command": [[exe, str(alt_ctl)], [exe, str(null_ctl)]],
        "working_directory": str(output_dir),
        "exit_code": alternative_run.returncode if alternative_run.returncode else null_run.returncode,
        "runtime_seconds": round(alternative_run.runtime_seconds + null_run.runtime_seconds, 4),
        "stdout": alternative_run.stdout + "\n" + null_run.stdout,
        "stderr": alternative_run.stderr + "\n" + null_run.stderr,
        "raw_output": alternative_text,
        "null_raw_output": null_text,
    }
