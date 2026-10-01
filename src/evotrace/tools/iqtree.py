"""IQ-TREE 2 maximum-likelihood integration and result parsing."""

import re
from pathlib import Path
from typing import Dict
from Bio import Phylo
from Bio import AlignIO
from evotrace.tools.common import locate, run


def detect_iqtree():
    exe = locate(("iqtree2", "iqtree"), "Install IQ-TREE 2 from https://www.iqtree.org/ and retry.")
    proc = run(exe, [exe, "-version"], timeout=30)
    version = (proc.stdout + "\n" + proc.stderr).strip().splitlines()
    return exe, (version[0] if version else "unknown")


def run_iqtree(
    alignment: Path, output_dir: Path, model="MFP", bootstrap=1000, threads=1, timeout=86400
) -> Dict[str, object]:
    alignment, output_dir = Path(alignment).resolve(), Path(output_dir).resolve()
    if not alignment.is_file():
        raise ValueError("Alignment file does not exist: {}".format(alignment))
    records = AlignIO.read(str(alignment), "fasta")
    if len(records) < 3 or len({len(r.seq) for r in records}) != 1:
        raise ValueError("IQ-TREE needs at least three equal-length aligned sequences.")
    if not re.fullmatch(r"[A-Za-z0-9+*_.-]+", model):
        raise ValueError("Invalid IQ-TREE model string.")
    output_dir.mkdir(parents=True, exist_ok=True)
    exe, version = detect_iqtree()
    prefix = output_dir / "iqtree"
    cmd = [
        exe,
        "-s",
        str(alignment),
        "-m",
        model,
        "-pre",
        str(prefix),
        "-nt",
        str(max(1, int(threads))),
    ]
    if bootstrap:
        if int(bootstrap) < 0:
            raise ValueError("bootstrap must be >=0")
        cmd += ["-B", str(int(bootstrap))]
    result = run(exe, cmd, cwd=str(output_dir), timeout=timeout)
    tree_path = Path(str(prefix) + ".treefile")
    report_path = Path(str(prefix) + ".iqtree")
    if not tree_path.is_file():
        raise RuntimeError(
            "IQ-TREE exited successfully but did not create {}".format(tree_path.name)
        )
    tree = Phylo.read(str(tree_path), "newick")
    text = report_path.read_text(errors="replace") if report_path.exists() else ""
    likelihood_match = re.search(r"Log-likelihood of the tree:\s*(-?\d+(?:\.\d+)?)", text)
    model_match = re.search(r"Best-fit model according to .*?:\s*(\S+)", text)
    return {
        "tree_path": str(tree_path),
        "report_path": str(report_path),
        "tree": tree,
        "log_likelihood": float(likelihood_match.group(1)) if likelihood_match else None,
        "selected_model": model_match.group(1) if model_match else model,
        "requested_model": model,
        "version": version,
        "command": cmd,
        "working_directory": result.working_directory,
        "exit_code": result.returncode,
        "runtime_seconds": result.runtime_seconds,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "support_values_present": any(t.confidence is not None for t in tree.get_nonterminals()),
    }
