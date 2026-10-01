"""Configuration loading with safe JSON/YAML support."""

import copy
import json
import re
from pathlib import Path
from typing import Any, Dict, Optional, cast

DEFAULTS: Dict[str, Any] = {
    "alphabet": "auto",
    "aligner": "mafft",
    "mafft_strategy": "auto",
    "max_ambiguity": 0.5,
    "selection": False,
    "motifs": [],
    "reference": None,
    "gap_policy": "ignore",
    "pairwise_dnds": False,
    "phylogeny": {"tool": "nj", "model": "MFP", "bootstrap": 1000},
    "pipeline": {"workers": 1},
    "selection_config": {},
    "domains_config": {},
    "input_config": {},
    "structure_config": {},
}

SECTIONS = {
    "project",
    "input",
    "alignment",
    "conservation",
    "selection",
    "domains",
    "structure",
    "report",
    "pipeline",
    "phylogeny",
}


def load_config(path: Optional[str]) -> Dict[str, Any]:
    cfg: Dict[str, Any] = copy.deepcopy(DEFAULTS)
    if not path:
        return cfg
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    if p.suffix.lower() == ".json":
        data = json.loads(text)
    else:
        try:
            import yaml
        except ImportError:
            raise ValueError("YAML config requires PyYAML. Install with: pip install evotrace[all]")
        data = yaml.safe_load(text)
    if not isinstance(data, dict):
        raise ValueError("Configuration root must be a mapping.")
    data = cast(Dict[str, Any], data)
    unknown = sorted(set(data) - set(DEFAULTS) - SECTIONS)
    if unknown:
        raise ValueError("Unknown configuration keys: {}".format(", ".join(unknown)))
    cfg.update(
        {
            k: v
            for k, v in data.items()
            if k in DEFAULTS
            and k
            not in (
                "phylogeny",
                "pipeline",
                "selection_config",
                "domains_config",
                "input_config",
                "structure_config",
            )
        }
    )
    # Accept the nested v1 configuration while retaining the original flat keys.
    alignment = data.get("alignment", {})
    if isinstance(alignment, dict):
        cfg["aligner"] = alignment.get("tool", cfg["aligner"])
        cfg["mafft_strategy"] = alignment.get("strategy", cfg["mafft_strategy"])
    conservation = data.get("conservation", {})
    if isinstance(conservation, dict):
        cfg["gap_policy"] = conservation.get("gap_policy", cfg["gap_policy"])
    if isinstance(data.get("input"), dict):
        cfg["input_config"] = data["input"]
        cfg["alphabet"] = data["input"].get("alphabet", cfg["alphabet"])
    if isinstance(data.get("phylogeny"), dict):
        cfg["phylogeny"].update(data["phylogeny"])
    if isinstance(data.get("pipeline"), dict):
        cfg["pipeline"].update(data["pipeline"])
    if isinstance(data.get("selection"), dict):
        cfg["selection_config"] = data["selection"]
        cfg["selection"] = bool(data["selection"].get("enabled", False))
        cfg["pairwise_dnds"] = bool(data["selection"].get("pairwise_dnds", False))
    elif isinstance(data.get("selection"), bool):
        cfg["selection"] = data["selection"]
        cfg["pairwise_dnds"] = data["selection"]
    if isinstance(data.get("domains"), dict):
        cfg["domains_config"] = data["domains"]
    if isinstance(data.get("structure"), dict):
        cfg["structure_config"] = data["structure"]
    motifs = data.get("motifs")
    if isinstance(motifs, dict):
        cfg["motifs"] = motifs.get("custom", [])
    if cfg["aligner"] != "mafft":
        raise ValueError(
            "Only MAFFT alignment is currently implemented; requested {!r}.".format(cfg["aligner"])
        )
    if cfg["gap_policy"] not in ("ignore", "include"):
        raise ValueError("conservation.gap_policy must be 'ignore' or 'include'.")
    if cfg["phylogeny"].get("tool", "nj") not in ("nj", "iqtree"):
        raise ValueError("phylogeny.tool must be 'nj' or 'iqtree'.")
    if int(cfg["phylogeny"].get("bootstrap", 1000)) < 0:
        raise ValueError("phylogeny.bootstrap must be non-negative.")
    if not isinstance(cfg.get("pipeline"), dict):
        raise ValueError("pipeline configuration must be a mapping.")
    if int(cfg["pipeline"].get("workers", 1)) < 1:
        raise ValueError("pipeline.workers must be a positive integer.")
    for key, section in (("phylogeny", cfg.get("phylogeny")), ("selection", cfg.get("selection_config")),
                         ("domains", cfg.get("domains_config")), ("structure", cfg.get("structure_config"))):
        if not isinstance(section, dict):
            raise ValueError("{} configuration must be a mapping.".format(key))
    selection = cfg.get("selection_config", {})
    if cfg.get("selection"):
        engine = selection.get("engine", "hyphy")
        if engine not in ("hyphy", "paml"):
            raise ValueError("selection.engine must be 'hyphy' or 'paml'.")
        if engine == "hyphy":
            methods = selection.get("analyses", ["FEL"])
            allowed = {"FEL", "MEME", "FUBAR", "ABSREL", "BUSTED"}
            if not isinstance(methods, list) or not methods:
                raise ValueError("selection.analyses must be a non-empty list.")
            normalized = [str(method).upper() for method in methods]
            if any(method not in allowed for method in normalized):
                raise ValueError("selection.analyses contains an unsupported HyPhy method.")
            if len(set(normalized)) != len(normalized):
                raise ValueError("selection.analyses must not contain duplicate methods.")
        elif selection.get("model", "branch_site") not in ("branch", "site", "branch_site"):
            raise ValueError("PAML selection.model must be branch, site, or branch_site.")
        elif int(selection.get("genetic_code", 1)) != 1:
            raise ValueError("PAML codeml currently supports only genetic_code=1 (standard code).")
    if not isinstance(cfg.get("motifs", []), list):
        raise ValueError("motifs must be a list of named regular-expression patterns.")
    for motif in cfg.get("motifs", []):
        if not isinstance(motif, dict) or not motif.get("name") or not motif.get("pattern"):
            raise ValueError("Each motif requires a name and regex pattern.")
        try:
            if re.compile(str(motif["pattern"])).match("") is not None:
                raise ValueError("Motif {} can match an empty sequence.".format(motif["name"]))
        except re.error as exc:
            raise ValueError("Invalid regex for motif {}: {}".format(motif["name"], exc)) from exc
    if not 0 <= float(cfg["max_ambiguity"]) <= 1:
        raise ValueError("max_ambiguity must be between 0 and 1.")
    return cfg
