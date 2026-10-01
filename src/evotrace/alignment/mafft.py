"""Safe MAFFT integration."""

import shutil
import subprocess
import time
import os
from pathlib import Path
from typing import List
from evotrace.io.fasta import read_fasta, write_fasta
from evotrace.models import SequenceRecord


class AlignmentError(RuntimeError):
    pass


def run_mafft(
    records: List[SequenceRecord],
    output: Path,
    strategy: str = "auto",
    threads: int = 1,
    timeout: int = 86400,
) -> dict:
    exe = shutil.which("mafft")
    if not exe:
        raise AlignmentError(
            "MAFFT was not found on PATH. Install MAFFT and retry; no alignment was fabricated."
        )
    if strategy not in ("auto", "linsi", "einsi", "ginsi"):
        raise AlignmentError("MAFFT strategy must be auto, linsi, einsi, or ginsi.")
    output.parent.mkdir(parents=True, exist_ok=True)
    input_path = output.parent / (output.stem + ".input.fasta")
    write_fasta(records, input_path)
    cmd = [exe, "--" + strategy, "--thread", str(max(1, int(threads))), "--quiet", str(input_path)]
    started = time.monotonic()
    cwd = os.getcwd()
    try:
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise AlignmentError("MAFFT could not complete: {}".format(exc))
    if proc.returncode:
        raise AlignmentError(
            "MAFFT exited with code {}: {}".format(proc.returncode, proc.stderr[-3000:])
        )
    output.write_text(proc.stdout, encoding="utf-8")
    aligned = read_fasta(output)
    lengths = {len(r.sequence) for r in aligned}
    if len(aligned) != len(records) or len(lengths) != 1:
        raise AlignmentError(
            "MAFFT output failed validation (record count or aligned lengths differ)."
        )
    version = "unknown"
    try:
        version_result = subprocess.run(
            [exe, "--version"], capture_output=True, text=True, timeout=10, check=False
        )
        version = (
            (version_result.stdout or version_result.stderr).strip().splitlines() or ["unknown"]
        )[0]
    except (OSError, subprocess.TimeoutExpired):
        pass
    return {
        "command": cmd,
        "working_directory": cwd,
        "exit_code": proc.returncode,
        "runtime_seconds": round(time.monotonic() - started, 4),
        "version": version,
        # MAFFT's stdout is the alignment artifact itself, not a diagnostic log.
        "stdout": "",
        "stderr": proc.stderr,
    }
