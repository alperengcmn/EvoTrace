"""HMMER hmmscan integration and domtblout parser."""

from pathlib import Path
from evotrace.tools.common import locate, run


def parse_domtblout(path: Path):
    hits = []
    for line in Path(path).read_text(errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split(maxsplit=22)
        if len(fields) < 22:
            continue
        try:
            hits.append(
                {
                    "domain_id": fields[0],
                    "domain_name": fields[0],
                    "sequence_id": fields[3],
                    "evalue": float(fields[12]),
                    "score": float(fields[13]),
                    "start": int(fields[17]),
                    "end": int(fields[18]),
                    "database": str(path),
                    "description": fields[22] if len(fields) > 22 else "",
                }
            )
        except (ValueError, IndexError):
            continue
    return hits


def run_hmmscan(fasta: Path, database: Path, output_dir: Path, threads=1, timeout=86400):
    fasta, database, output_dir = (
        Path(fasta).resolve(),
        Path(database).resolve(),
        Path(output_dir).resolve(),
    )
    if not fasta.is_file() or not database.is_file():
        raise ValueError("hmmscan requires existing sequence FASTA and HMM database paths.")
    output_dir.mkdir(parents=True, exist_ok=True)
    exe = locate(("hmmscan",), "Install HMMER3 and make hmmscan available on PATH.")
    domtbl = output_dir / "domains.domtblout"
    cmd = [
        exe,
        "--cpu",
        str(max(1, int(threads))),
        "--domtblout",
        str(domtbl),
        "--noali",
        str(database),
        str(fasta),
    ]
    result = run(exe, cmd, cwd=str(output_dir), timeout=timeout)
    if not domtbl.is_file():
        raise RuntimeError("hmmscan completed without a domain table.")
    ver = run(exe, [exe, "-h"], timeout=30).stdout.splitlines()
    return {
        "hits": parse_domtblout(domtbl),
        "domtblout": str(domtbl),
        "version": next((x.strip() for x in ver if "HMMER" in x), "unknown"),
        "command": cmd,
        "working_directory": result.working_directory,
        "exit_code": result.returncode,
        "runtime_seconds": result.runtime_seconds,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }
