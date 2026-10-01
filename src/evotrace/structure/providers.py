"""Real PDB/AlphaFold retrieval and gap-aware residue-to-structure mapping."""

import json
import re
import time
import urllib.request
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict

from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1_extended


class StructureUnavailable(RuntimeError):
    pass


class StructureProvider(ABC):
    @abstractmethod
    def get_structure(self) -> Dict[str, object]:
        """Return a structure source status and, when available, a local coordinate path."""


class PDBProvider(StructureProvider):
    def __init__(self, path: Path):
        self.path = Path(path)

    def get_structure(self):
        if not self.path.is_file():
            return {
                "status": "structure_unavailable",
                "reason": "PDB/mmCIF file not found: {}".format(self.path),
            }
        return {"status": "structure_found", "path": str(self.path.resolve()), "retrieved": False}


class AlphaFoldProvider(StructureProvider):
    """Fetch public AlphaFold DB model metadata and mmCIF; no credentials are used."""

    def __init__(self, accession: str, cache_dir: Path, cache_ttl_days=30, timeout=30):
        if not re.fullmatch(r"[A-Za-z0-9-]+", accession):
            raise ValueError("Invalid UniProt accession for AlphaFold retrieval.")
        self.accession = accession
        self.cache_dir = Path(cache_dir)
        self.cache_ttl_days = cache_ttl_days
        self.timeout = timeout

    def get_structure(self):
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        cif = self.cache_dir / (self.accession + ".cif")
        meta = self.cache_dir / (self.accession + ".json")
        if cif.is_file() and meta.is_file():
            metadata = json.loads(meta.read_text(encoding="utf-8"))
            age = (time.time() - float(metadata["retrieved_unix"])) / 86400
            if age <= self.cache_ttl_days:
                return {
                    "status": "structure_found",
                    "path": str(cif),
                    "retrieved": False,
                    "cache_age_days": age,
                    "metadata": metadata,
                }
        url = "https://alphafold.ebi.ac.uk/api/prediction/{}".format(self.accession)
        request = urllib.request.Request(
            url, headers={"User-Agent": "EvoTrace/0.1 public-structure-client"}
        )
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                payload = json.loads(response.read().decode("utf-8"))
        except Exception as exc:
            raise StructureUnavailable(
                "AlphaFold DB lookup failed for {}: {}".format(self.accession, exc)
            ) from exc
        if not isinstance(payload, list) or not payload or not payload[0].get("cifUrl"):
            return {
                "status": "structure_unavailable",
                "accession": self.accession,
                "reason": "No AlphaFold DB model was returned for this accession.",
            }
        cif_url = payload[0]["cifUrl"]
        try:
            req = urllib.request.Request(
                cif_url, headers={"User-Agent": "EvoTrace/0.1 public-structure-client"}
            )
            with urllib.request.urlopen(req, timeout=self.timeout) as response:
                body = response.read()
        except Exception as exc:
            raise StructureUnavailable("AlphaFold mmCIF download failed: {}".format(exc)) from exc
        tmp = cif.with_suffix(".cif.tmp")
        tmp.write_bytes(body)
        tmp.replace(cif)
        metadata = {
            "provider": "AlphaFold DB",
            "accession": self.accession,
            "metadata_url": url,
            "cif_url": cif_url,
            "retrieved_unix": time.time(),
            "model_identifier": payload[0].get("entryId"),
        }
        meta.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        return {
            "status": "structure_found",
            "path": str(cif),
            "retrieved": True,
            "metadata": metadata,
        }


def map_sequence_to_structure(
    sequence: str, structure_path: Path, chain_id=None
) -> Dict[str, object]:
    """Map 1-based query residues to real chain/residue IDs using global alignment."""
    sequence = sequence.replace("-", "").replace(".", "").upper()
    if not sequence:
        return {"status": "mapping_failed", "reason": "Query protein sequence is empty."}
    path = Path(structure_path)
    parser = (
        MMCIFParser(QUIET=True)
        if path.suffix.lower() in (".cif", ".mmcif")
        else PDBParser(QUIET=True)
    )
    structure = parser.get_structure("evotrace", str(path))
    model = next(structure.get_models())
    chains = []
    for chain in model:
        if chain_id and chain.id != chain_id:
            continue
        letters, residues = [], []
        for residue in chain:
            if residue.id[0] != " ":
                continue
            aa = protein_letters_3to1_extended.get(residue.resname.upper(), "X")
            letters.append(aa)
            residues.append(residue)
        if letters:
            chains.append((chain.id, "".join(letters), residues))
    if not chains:
        return {
            "status": "mapping_failed",
            "reason": "No protein residues found in requested structure chain.",
        }
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -6
    aligner.extend_gap_score = -0.5
    candidates = []
    for cid, struct_seq, residues in chains:
        alignment = aligner.align(sequence, struct_seq)[0]
        aligned_query, aligned_structure = alignment.aligned
        matches = sum(
            sequence[q0 + i] == struct_seq[s0 + i]
            for (q0, q1), (s0, s1) in zip(aligned_query, aligned_structure)
            for i in range(min(q1 - q0, s1 - s0))
        )
        aligned_count = sum(
            min(q1 - q0, s1 - s0) for (q0, q1), (s0, s1) in zip(aligned_query, aligned_structure)
        )
        identity = matches / aligned_count if aligned_count else 0.0
        candidates.append(
            (alignment.score, identity, cid, struct_seq, residues, aligned_query, aligned_structure)
        )
    best = max(candidates, key=lambda row: row[0])
    _, identity, cid, struct_seq, residues, aq, ass = best
    mapping = []
    for (q0, q1), (s0, s1) in zip(aq, ass):
        for qpos, spos in zip(range(q0, q1), range(s0, s1)):
            residue = residues[spos]
            atoms = list(residue.get_atoms())
            ca = residue["CA"] if "CA" in residue else None
            mapping.append(
                {
                    "sequence_position": qpos + 1,
                    "structure_chain": cid,
                    "structure_residue_number": int(residue.id[1]),
                    "structure_insertion_code": str(residue.id[2]).strip(),
                    "structure_residue": struct_seq[spos],
                    "identity_match": sequence[qpos] == struct_seq[spos],
                    "ca_coordinates_angstrom": list(ca.coord) if ca is not None else None,
                    "mean_b_factor": sum(float(atom.bfactor) for atom in atoms) / len(atoms)
                    if atoms
                    else None,
                }
            )
    coverage = len(mapping) / max(1, len(sequence))
    return {
        "status": "mapping_successful",
        "structure_path": str(path.resolve()),
        "chain": cid,
        "sequence_identity": identity,
        "query_coverage": coverage,
        "mapped_residues": mapping,
        "unmapped_positions": sorted(
            set(range(1, len(sequence) + 1)) - {x["sequence_position"] for x in mapping}
        ),
    }
