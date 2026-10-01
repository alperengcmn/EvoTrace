"""Local dependency diagnostics with actionable optional-tool guidance."""

import importlib
import platform
import shutil
import subprocess
import sys


TOOLS = {
    "MAFFT": (("mafft",), "MAFFT: https://mafft.cbrc.jp/alignment/software/", True),
    "IQ-TREE": (("iqtree2", "iqtree", "iqtree3"), "IQ-TREE: https://www.iqtree.org/", True),
    "HyPhy": (("hyphy",), "HyPhy: https://hyphy.org/", False),
    "PAML/codeml": (("codeml",), "PAML: http://abacus.gene.ucl.ac.uk/software/paml.html", False),
    "HMMER/hmmscan": (("hmmscan",), "HMMER: http://hmmer.org/", True),
    "HMMER/hmmsearch": (("hmmsearch",), "HMMER: http://hmmer.org/", False),
    "BLAST+": (
        ("blastp", "blastn"),
        "BLAST+: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download",
        False,
    ),
}
PYTHON_MODULES = {
    "Biopython": "Bio",
    "NumPy": "numpy",
    "PyYAML": "yaml",
    "pandas": "pandas",
    "Plotly": "plotly",
    "Streamlit": "streamlit",
    "pytest": "pytest",
    "Ruff": "ruff",
    "mypy": "mypy",
}


def _version(command):
    try:
        p = subprocess.run(command, capture_output=True, text=True, timeout=4, check=False)
        return ((p.stdout or p.stderr).strip().splitlines() or ["available"])[0]
    except Exception:
        return "available (version query failed)"


def diagnostics():
    rows = [
        {
            "name": "Operating system",
            "available": True,
            "required": False,
            "version": platform.platform(),
            "path": "",
            "guidance": "",
            "architecture": platform.machine(),
        },
        {
            "name": "Python",
            "available": True,
            "required": True,
            "version": sys.version.split()[0],
            "path": sys.executable,
            "guidance": "",
            "architecture": platform.machine(),
        }
    ]
    for name, command in (("Conda", "conda"), ("Homebrew", "brew"), ("Docker", "docker"), ("Git", "git")):
        path = shutil.which(command)
        rows.append({
            "name": name, "available": bool(path), "required": False,
            "version": _version([path, "--version"]) if path else "not found",
            "path": path, "guidance": "", "architecture": "",
        })
    for name, (executables, guidance, required) in TOOLS.items():
        path = next((shutil.which(x) for x in executables if shutil.which(x)), None)
        rows.append(
            {
                "name": name,
                "available": bool(path),
                "required": required,
                "version": (
                    "available (no version flag)"
                    if path and name == "PAML/codeml"
                    else _version([path, "-h"] if path and name.startswith("HMMER/") else [path, "-version"] if path and name == "BLAST+" else [path, "--version"])
                    if path
                    else "not found"
                ),
                "path": path,
                "guidance": "" if path else guidance,
                "architecture": _version(["file", "-b", path]) if path and shutil.which("file") else "unknown",
            }
        )
    for name, module in PYTHON_MODULES.items():
        try:
            m = importlib.import_module(module)
            rows.append(
                {
                    "name": name,
                    "available": True,
                    "version": getattr(m, "__version__", "installed"),
                    "guidance": "",
                    "path": getattr(m, "__file__", ""),
                    "required": name in ("Biopython", "NumPy", "PyYAML", "pandas", "Plotly", "Streamlit"),
                    "architecture": "",
                }
            )
        except ImportError:
            required = name in ("Biopython", "NumPy", "PyYAML", "pandas", "Plotly", "Streamlit")
            rows.append(
                {
                    "name": name,
                    "available": False,
                    "version": "not installed",
                    "guidance": "pip install evotrace[web]"
                    if name in ("Plotly", "Streamlit")
                    else "pip install evotrace[dev]"
                    if name in ("Ruff", "mypy", "pytest")
                    else "pip install pyyaml",
                    "required": required,
                    "path": "",
                    "architecture": "",
                }
            )
    rows.append({
        "name": "EvoTrace package", "available": importlib.util.find_spec("evotrace") is not None,
        "required": True, "version": "installed" if importlib.util.find_spec("evotrace") else "not installed",
        "path": str(importlib.util.find_spec("evotrace").origin) if importlib.util.find_spec("evotrace") else "",
        "guidance": "Install with: python -m pip install -e '.[all]'" if importlib.util.find_spec("evotrace") is None else "",
        "architecture": "",
    })
    return rows


def print_doctor(strict=False):
    print(
        "EvoTrace dependency check (✓ available, ✗ missing; missing optional tools do not block the core workflow)"
    )
    print(
        "{:<18} {:<10} {:<9} {:<30} {}".format(
            "Tool", "Status", "Required", "Version", "Architecture / path / guidance"
        )
    )
    rows = diagnostics()
    for row in rows:
        mark = "✓" if row["available"] else "✗"
        location = "{} {}".format(row.get("architecture", ""), row.get("path") or row.get("guidance", "")).strip()
        print(
            "{:<18} {:<10} {:<9} {:<30} {}".format(
                row["name"],
                mark,
                "yes" if row.get("required") else "no",
                row.get("version", ""),
                location,
            )
        )
    missing = [row["name"] for row in rows if row.get("required") and not row["available"]]
    if strict and missing:
        print("Strict check failed; missing required components: {}".format(", ".join(missing)))
        return 1
    if missing:
        print("Required components missing: {}".format(", ".join(missing)))
    return 0
