"""Dependency-aware pipeline stages with content-addressed artifact reuse."""

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Callable, Dict, Iterable

STATUSES = {"PENDING", "RUNNING", "SUCCESS", "FAILED", "SKIPPED"}


class PipelineEngine:
    def __init__(self, output_dir: Path, resume=False):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.manifest = self.output_dir / "pipeline_stages.json"
        self.resume = resume
        try:
            self.state = (
                json.loads(self.manifest.read_text()) if self.manifest.exists() else {"stages": {}}
            )
        except (ValueError, OSError):
            self.state = {"stages": {}}

    def _save(self):
        tmp = self.manifest.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(self.state, indent=2, sort_keys=True), encoding="utf-8")
        os.replace(str(tmp), str(self.manifest))

    def run_stage(
        self,
        name: str,
        action: Callable[[], object],
        inputs: Dict[str, object],
        outputs: Iterable[Path],
        dependencies=(),
        parameters=None,
        software_version="internal",
        command=None,
    ):
        output_paths = [Path(p).resolve() for p in outputs]
        for dependency in dependencies:
            dep = self.state["stages"].get(dependency, {})
            if dep.get("status") not in ("SUCCESS", "SKIPPED"):
                raise RuntimeError(
                    "Stage {!r} requires successful stage {!r}.".format(name, dependency)
                )
        key_data = {
            "inputs": inputs,
            "parameters": parameters or {},
            "software_version": software_version,
        }
        try:
            fingerprint = hashlib.sha256(
                json.dumps(key_data, sort_keys=True, default=str).encode()
            ).hexdigest()
        except (TypeError, ValueError) as exc:
            raise ValueError("Stage cache parameters must be JSON serializable: {}".format(exc))
        previous = self.state["stages"].get(name, {})
        previous_hashes = previous.get("output_sha256", {})
        if (
            self.resume
            and previous.get("status") == "SUCCESS"
            and previous.get("fingerprint") == fingerprint
            and all(p.is_file() for p in output_paths)
            and all(previous_hashes.get(str(p)) == self._sha256(p) for p in output_paths)
        ):
            previous["status"] = "SUCCESS"
            previous["cache_hit"] = True
            self._save()
            return previous.get("result")
        stage = {
            "name": name,
            "status": "RUNNING",
            "inputs": inputs,
            "outputs": [str(p) for p in output_paths],
            "parameters": parameters or {},
            "software_version": software_version,
            "command": command,
            "started_unix": time.time(),
            "fingerprint": fingerprint,
            "cache_hit": False,
        }
        self.state["stages"][name] = stage
        self._save()
        start = time.monotonic()
        try:
            result = action()
            missing = [str(p) for p in output_paths if not p.is_file()]
            if missing:
                raise RuntimeError(
                    "Stage {!r} did not create declared outputs: {}".format(
                        name, ", ".join(missing)
                    )
                )
            stage["status"] = "SUCCESS"
            stage["output_sha256"] = {str(p): self._sha256(p) for p in output_paths}
            if isinstance(result, dict):
                for key in (
                    "command",
                    "working_directory",
                    "exit_code",
                    "runtime_seconds",
                    "version",
                ):
                    if key in result:
                        stage[key] = result[key]
                logs = {}
                for key in ("stdout", "stderr"):
                    if isinstance(result.get(key), str):
                        logs[key] = result[key][-64000:]
                        logs[key + "_truncated"] = len(result[key]) > 64000
                if logs:
                    stage["logs"] = logs
            if isinstance(result, dict):
                stage["result"] = {
                    key: value
                    for key, value in result.items()
                    if key not in ("stdout", "stderr", "raw_output")
                }
            else:
                stage["result"] = (
                    result
                    if isinstance(result, (str, int, float, bool, type(None), list))
                    else None
                )
            stage["error"] = None
        except Exception as exc:
            stage["status"] = "FAILED"
            stage["error"] = "{}: {}".format(type(exc).__name__, exc)
            raise
        finally:
            stage["runtime_seconds"] = round(time.monotonic() - start, 4)
            stage["finished_unix"] = time.time()
            self._save()
        return stage.get("result")

    @staticmethod
    def _sha256(path):
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    def skip_stage(self, name: str, reason: str, inputs=None, dependencies=()):
        for dependency in dependencies:
            if self.state["stages"].get(dependency, {}).get("status") not in ("SUCCESS", "SKIPPED"):
                raise RuntimeError("Stage {!r} requires {!r}.".format(name, dependency))
        self.state["stages"][name] = {
            "name": name,
            "status": "SKIPPED",
            "inputs": inputs or {},
            "outputs": [],
            "parameters": {},
            "runtime_seconds": 0,
            "error": reason,
        }
        self._save()

    def summarize(self):
        return {name: stage.get("status") for name, stage in self.state["stages"].items()}
