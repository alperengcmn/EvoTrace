"""Shared safe executable discovery and subprocess records."""

import shutil
import subprocess
import time
from dataclasses import dataclass
from typing import Sequence


@dataclass
class ToolResult:
    executable: str
    version: str
    command: list
    returncode: int
    stdout: str
    stderr: str
    working_directory: str
    runtime_seconds: float


class ToolUnavailable(RuntimeError):
    pass


def locate(names: Sequence[str], guidance: str) -> str:
    for name in names:
        path = shutil.which(name)
        if path:
            return path
    raise ToolUnavailable("{} executable was not found on PATH. {}".format(names[0], guidance))


def run(executable: str, command: Sequence[str], cwd=None, timeout=86400) -> ToolResult:
    started = time.monotonic()
    try:
        proc = subprocess.run(
            list(command), cwd=cwd, capture_output=True, text=True, timeout=timeout, check=False
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            "Command timed out after {} seconds: {}".format(timeout, " ".join(command))
        ) from exc
    except OSError as exc:
        raise RuntimeError("Could not execute {}: {}".format(executable, exc)) from exc
    result = ToolResult(
        executable,
        "unknown",
        list(command),
        proc.returncode,
        proc.stdout or "",
        proc.stderr or "",
        str(cwd) if cwd is not None else str(__import__("os").getcwd()),
        round(time.monotonic() - started, 4),
    )
    if proc.returncode:
        raise RuntimeError(
            "{} exited with code {}. stderr: {}".format(
                executable, proc.returncode, result.stderr[-4000:]
            )
        )
    return result
