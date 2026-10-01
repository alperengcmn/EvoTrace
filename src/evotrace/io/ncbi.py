"""NCBI E-utilities accession retrieval with explicit, timestamped local caching."""

import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List
from evotrace.io.fasta import FastaError, read_fasta, write_fasta

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


class AccessionFetchError(RuntimeError):
    pass


def _request(accessions: List[str], database: str, timeout: int) -> str:
    query = urllib.parse.urlencode(
        {"db": database, "id": ",".join(accessions), "rettype": "fasta", "retmode": "text"}
    )
    request = urllib.request.Request(
        BASE_URL + "?" + query, headers={"User-Agent": "EvoTrace/0.1 (NCBI E-utilities client)"}
    )
    last: Exception = RuntimeError("No request attempt has run.")
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                text = bytes(response.read()).decode("utf-8")
                if not text.lstrip().startswith(">"):
                    raise AccessionFetchError("NCBI response did not contain FASTA records.")
                return text
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code not in (429, 500, 502, 503, 504):
                break
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
        if attempt < 3:
            time.sleep(min(2**attempt, 8))
    raise AccessionFetchError("NCBI E-utilities request failed: {}".format(last))


def fetch_accessions(
    accessions: Iterable[str],
    output_fasta: Path,
    cache_dir: Path,
    database="protein",
    cache_ttl_days=7,
    timeout=30,
) -> Dict[str, object]:
    """Download one or more NCBI nucleotide/protein accessions or reuse fresh cache.

    Stale cache is never silently substituted after a network failure. A caller may
    explicitly choose offline operation by passing only cached results in a future API.
    """
    ids = list(dict.fromkeys(str(x).strip() for x in accessions))
    if not ids or any(not re.fullmatch(r"[A-Za-z0-9_.]+", x) for x in ids):
        raise ValueError(
            "Provide one or more valid NCBI accession identifiers (letters, digits, dot, underscore)."
        )
    if database not in ("protein", "nuccore"):
        raise ValueError("NCBI database must be 'protein' or 'nuccore'.")
    cache_dir, output_fasta = Path(cache_dir), Path(output_fasta)
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256((database + "\0" + "\0".join(ids)).encode()).hexdigest()
    cached_fasta, metadata_path = cache_dir / (key + ".fasta"), cache_dir / (key + ".json")
    now = time.time()
    if cached_fasta.is_file() and metadata_path.is_file():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        age_days = (now - float(metadata["retrieved_unix"])) / 86400
        if age_days <= cache_ttl_days:
            records = read_fasta(cached_fasta)
            write_fasta(records, output_fasta)
            return {
                "status": "cache_hit",
                "database": database,
                "accessions": ids,
                "retrieved_utc": metadata["retrieved_utc"],
                "cache_age_days": age_days,
                "cache_key": key,
            }
    text = _request(ids, database, timeout)
    temporary = cache_dir / (key + ".download.tmp")
    temporary.write_text(text, encoding="utf-8")
    try:
        records = read_fasta(temporary)
    except FastaError as exc:
        temporary.unlink(missing_ok=True)
        raise AccessionFetchError("Could not parse NCBI FASTA response: {}".format(exc))
    observed = {record.id.split("|")[-1] for record in records}
    missing = [
        acc
        for acc in ids
        if acc not in observed and not any(acc in rec.description for rec in records)
    ]
    if missing:
        temporary.unlink(missing_ok=True)
        raise AccessionFetchError(
            "NCBI returned no matching FASTA record for: {}".format(", ".join(missing))
        )
    temporary.replace(cached_fasta)
    retrieved = datetime.now(timezone.utc).isoformat()
    metadata = {
        "database": database,
        "accessions": ids,
        "retrieved_utc": retrieved,
        "retrieved_unix": now,
        "source": BASE_URL,
        "cache_key": key,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    write_fasta(records, output_fasta)
    return {"status": "downloaded", **metadata}
