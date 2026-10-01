#!/usr/bin/env python3
"""Independent external validation of EvoTrace alignment constraint summaries.

Downloads ProteinGym v1.3 DMS data temporarily, selects a deterministic panel
without reading assay outcomes, samples its supplied A2M homolog alignments,
and compares EvoTrace's per-site constraint summary to experimental site fitness.
No files larger than the compact results are retained in the project.
"""

import argparse
import csv
import hashlib
import io
import json
import math
import random
import re
import tempfile
import zipfile
from datetime import date
from collections import OrderedDict, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

from Bio import SeqIO

from evotrace.conservation.core import columns
from evotrace.mapping.residues import AlignmentMap
from evotrace.models import SequenceRecord
from evotrace.scoring.ets import evolutionary_scores


BASE_URL = "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/"
DMS_ARCHIVE = BASE_URL + "DMS_ProteinGym_substitutions.zip"
DMS_METADATA = BASE_URL + "DMS_substitutions.csv"
MSA_ARCHIVE = BASE_URL + "DMS_msa_files.zip"
VERSION = "ProteinGym v1.3"
CATEGORIES = ("Activity", "OrganismalFitness", "Stability")
AA = set("ACDEFGHIKLMNPQRSTVWY")
MUTATION = re.compile(r"^([A-Z])([1-9][0-9]*)([A-Z])$")


class HTTPRangeReader(io.RawIOBase):
    """Seekable, small-block reader for ZIP central directories and members."""

    block_size = 2 * 1024 * 1024

    def __init__(self, url, size):
        self.url, self.size, self.position = url, size, 0
        self.cache = OrderedDict()

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.position

    def seek(self, offset, whence=0):
        if whence == 0:
            self.position = offset
        elif whence == 1:
            self.position += offset
        elif whence == 2:
            self.position = self.size + offset
        else:
            raise ValueError("Invalid seek mode.")
        return self.position

    def _block(self, index):
        if index not in self.cache:
            start = index * self.block_size
            end = min(self.size - 1, start + self.block_size - 1)
            request = Request(self.url, headers={"Range": "bytes={}-{}".format(start, end)})
            with urlopen(request, timeout=120) as response:
                if response.status != 206:
                    raise RuntimeError("Benchmark server did not honor HTTP byte-range requests.")
                self.cache[index] = response.read()
            while len(self.cache) > 8:
                self.cache.popitem(last=False)
        self.cache.move_to_end(index)
        return self.cache[index]

    def read(self, count=-1):
        if count < 0:
            count = self.size - self.position
        count = min(count, self.size - self.position)
        chunks = []
        while count:
            index = self.position // self.block_size
            offset = self.position % self.block_size
            block = self._block(index)
            take = min(count, len(block) - offset)
            chunks.append(block[offset : offset + take])
            self.position += take
            count -= take
        return b"".join(chunks)

    def readinto(self, buffer):
        data = self.read(len(buffer))
        buffer[: len(data)] = data
        return len(data)


def _sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _download(url, path):
    request = Request(url, headers={"User-Agent": "EvoTrace external validation/1.0"})
    with urlopen(request, timeout=120) as response, Path(path).open("wb") as output:
        while True:
            block = response.read(1024 * 1024)
            if not block:
                break
            output.write(block)


def _num(row, name, fallback=0.0):
    try:
        value = float(row.get(name, ""))
        return value if math.isfinite(value) else fallback
    except (TypeError, ValueError):
        return fallback


def _select_panel(metadata, msa_names, per_category=8):
    eligible = []
    for row in metadata:
        if (
            row.get("MSA_filename") not in msa_names
            or row.get("coarse_selection_type") not in CATEGORIES
            or not 50 <= _num(row, "seq_len") <= 350
            or _num(row, "DMS_number_single_mutants") < 1000
            or _num(row, "MSA_num_seqs") < 100
            or _num(row, "MSA_num_seqs") > 200000
            or _num(row, "MSA_perc_cov") < 0.80
        ):
            continue
        eligible.append(row)
    selected, used_targets = [], set()
    for category in CATEGORIES:
        candidates = [r for r in eligible if r["coarse_selection_type"] == category]
        candidates.sort(
            key=lambda r: hashlib.sha256(("20261001:" + r["DMS_id"]).encode()).hexdigest()
        )
        count = 0
        for row in candidates:
            if row["UniProt_ID"] in used_targets:
                continue
            selected.append(row)
            used_targets.add(row["UniProt_ID"])
            count += 1
            if count == per_category:
                break
    return selected


def _normalize_a2m(sequence):
    return "".join("-" if char in ".-" else char.upper() for char in sequence)


def _sample_msa(zip_file, member, accession, target, sample_size, seed):
    rng = random.Random(seed)
    reservoir, target_aligned, target_header = [], None, None
    width = None
    seen_records = 0
    with zip_file.open(member) as binary:
        text = io.TextIOWrapper(binary, encoding="utf-8")
        for record in SeqIO.parse(text, "fasta"):
            aligned = _normalize_a2m(str(record.seq))
            if width is None:
                width = len(aligned)
            if len(aligned) != width:
                raise ValueError("A2M records have inconsistent alignment widths.")
            header_id = record.id.split("/", 1)[0]
            if header_id == accession:
                target_aligned, target_header = aligned, record.id
                continue
            seen_records += 1
            if len(reservoir) < sample_size:
                reservoir.append(aligned)
            else:
                replace = rng.randrange(seen_records)
                if replace < sample_size:
                    reservoir[replace] = aligned
    if target_aligned is None:
        raise ValueError("No target sequence matching {} was found in its MSA.".format(accession))
    if target_aligned.replace("-", "") != target:
        raise ValueError("ProteinGym MSA target does not match its DMS target sequence.")
    if not reservoir:
        raise ValueError("No non-target homologs found in the MSA.")
    records = [SequenceRecord("target", target_aligned)]
    records.extend(
        SequenceRecord("homolog_{:04d}".format(i), sequence)
        for i, sequence in enumerate(reservoir, 1)
    )
    return records, target_header, seen_records


def _rank(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        average = (start + 1 + end) / 2.0
        for i in range(start, end):
            ranks[order[i]] = average
        start = end
    return ranks


def _spearman(xs, ys):
    if len(xs) < 3:
        return None
    rx, ry = _rank(xs), _rank(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    covariance = sum((x - mx) * (y - my) for x, y in zip(rx, ry))
    vx = sum((x - mx) ** 2 for x in rx)
    vy = sum((y - my) ** 2 for y in ry)
    if not vx or not vy:
        return None
    return covariance / math.sqrt(vx * vy)


def _bootstrap(values, seed=20261001, iterations=20000):
    rng = random.Random(seed)
    means = []
    for _ in range(iterations):
        means.append(sum(rng.choice(values) for _ in values) / len(values))
    means.sort()
    return [means[int(0.025 * iterations)], means[int(0.975 * iterations)]]


def _sign_test_p(negative, nonzero):
    if not nonzero:
        return None
    return sum(math.comb(nonzero, k) for k in range(negative, nonzero + 1)) / (2**nonzero)


def _analyze_assay(row, assay_csv, msa_zip, sample_size):
    accession = row["UniProt_ID"]
    target = row["target_seq"].strip().upper()
    member = "DMS_msa_files/" + row["MSA_filename"]
    info = msa_zip.getinfo(member)
    seed = int(hashlib.sha256(("20261001:" + row["DMS_id"]).encode()).hexdigest()[:8], 16)
    alignment, target_header, available = _sample_msa(
        msa_zip, member, accession, target, sample_size, seed
    )
    col_metrics = columns(alignment, alphabet="protein", gap_policy="ignore")
    evo_scores = evolutionary_scores(col_metrics)
    coordinate_map = AlignmentMap(alignment, "target")
    with io.TextIOWrapper(assay_csv, encoding="utf-8", newline="") as text:
        data = list(csv.DictReader(text))
    per_site = defaultdict(list)
    mapped_variants = 0
    for variant in data:
        match = MUTATION.fullmatch(variant.get("mutant", ""))
        if not match:
            continue
        wildtype, position, alternate = match.group(1), int(match.group(2)), match.group(3)
        if alternate not in AA or wildtype not in AA or position > len(target):
            continue
        if target[position - 1] != wildtype or alternate == wildtype:
            continue
        try:
            alignment_position = coordinate_map.sequence_to_alignment("target", position)
        except ValueError:
            continue
        try:
            score = float(variant.get("DMS_score", ""))
            fit = float(variant.get("DMS_score_bin", ""))
        except (TypeError, ValueError):
            continue
        if not math.isfinite(score) or fit not in (0.0, 1.0):
            continue
        per_site[position].append((score, fit, alignment_position))
        mapped_variants += 1

    x_constraint, x_frequency, y_fit, y_score = [], [], [], []
    site_records = []
    for position, variants in sorted(per_site.items()):
        if len(variants) < 5:
            continue
        alignment_position = variants[0][2]
        site_score = evo_scores[alignment_position - 1]
        column = col_metrics[alignment_position - 1]
        fit_fraction = sum(v[1] for v in variants) / len(variants)
        mean_dms = sum(v[0] for v in variants) / len(variants)
        x_constraint.append(site_score["constraint_score"])
        x_frequency.append(column["conservation_score"])
        y_fit.append(fit_fraction)
        y_score.append(mean_dms)
        site_records.append(
            {
                "target_position": position,
                "alignment_position": alignment_position,
                "reference_residue": target[position - 1],
                "variant_count": len(variants),
                "fit_fraction": fit_fraction,
                "mean_dms_score": mean_dms,
                "constraint": site_score["constraint_score"],
                "conservation_frequency": column["conservation_score"],
                "occupancy": column["occupancy"],
            }
        )
    if len(site_records) < 10:
        raise ValueError("Fewer than 10 adequately measured sites after coordinate/QC filtering.")
    return {
        "dms_id": row["DMS_id"],
        "uniprot_id": accession,
        "selection_type": row["coarse_selection_type"],
        "sequence_length": len(target),
        "msa_name": row["MSA_filename"],
        "msa_crc32": "{:08x}".format(info.CRC),
        "msa_uncompressed_bytes": info.file_size,
        "protein_gym_msa_depth": int(_num(row, "MSA_num_seqs")),
        "sampled_homologs": len(alignment) - 1,
        "available_non_target_records": available,
        "target_header": target_header,
        "mapped_single_missense_variants": mapped_variants,
        "evaluated_sites": len(site_records),
        "rho_constraint_vs_fit_fraction": _spearman(x_constraint, y_fit),
        "rho_conservation_vs_fit_fraction": _spearman(x_frequency, y_fit),
        "rho_constraint_vs_mean_dms_score": _spearman(x_constraint, y_score),
        "sites": site_records,
    }


def run_validation(output_dir, per_category=8, sample_size=500):
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="evotrace-proteingym-") as temp:
        temp = Path(temp)
        dms_path, metadata_path = temp / "dms.zip", temp / "metadata.csv"
        _download(DMS_ARCHIVE, dms_path)
        _download(DMS_METADATA, metadata_path)
        metadata = list(csv.DictReader(metadata_path.open(encoding="utf-8", newline="")))
        range_request = Request(MSA_ARCHIVE, method="HEAD")
        with urlopen(range_request, timeout=60) as response:
            msa_bytes = int(response.headers["Content-Length"])
        msa_zip = zipfile.ZipFile(HTTPRangeReader(MSA_ARCHIVE, msa_bytes))
        msa_names = {Path(name).name for name in msa_zip.namelist() if not name.endswith("/")}
        panel = _select_panel(metadata, msa_names, per_category)
        if not panel:
            raise RuntimeError("No target passed the pre-specified ProteinGym panel filters.")

        assays = []
        with zipfile.ZipFile(dms_path) as dms_zip:
            assay_names = set(dms_zip.namelist())
            for row in panel:
                member = "DMS_ProteinGym_substitutions/" + row["DMS_filename"]
                if member not in assay_names:
                    continue
                try:
                    with dms_zip.open(member) as assay_csv:
                        result = _analyze_assay(row, assay_csv, msa_zip, sample_size)
                    result["dms_member_crc32"] = "{:08x}".format(dms_zip.getinfo(member).CRC)
                    assays.append(result)
                    print("Validated {}: {} sites, {} variants, rho={:.3f}".format(
                        result["dms_id"], result["evaluated_sites"],
                        result["mapped_single_missense_variants"],
                        result["rho_constraint_vs_fit_fraction"],
                    ), flush=True)
                except Exception as exc:
                    assays.append(
                        {
                            "dms_id": row["DMS_id"],
                            "uniprot_id": row["UniProt_ID"],
                            "selection_type": row["coarse_selection_type"],
                            "status": "excluded_after_qc",
                            "reason": "{}: {}".format(type(exc).__name__, exc),
                            "msa_name": row["MSA_filename"],
                        }
                    )
                    print("Excluded {} after QC: {}".format(row["DMS_id"], exc), flush=True)

    valid = [r for r in assays if r.get("rho_constraint_vs_fit_fraction") is not None]
    if not valid:
        raise RuntimeError("No assay remained after independent data and coordinate QC.")
    correlations = [r["rho_constraint_vs_fit_fraction"] for r in valid]
    negative = sum(value < 0 for value in correlations)
    nonzero = sum(value != 0 for value in correlations)
    mean_rho = sum(correlations) / len(correlations)
    by_category = {}
    for category in CATEGORIES:
        vals = [r["rho_constraint_vs_fit_fraction"] for r in valid if r["selection_type"] == category]
        if vals:
            by_category[category] = {"n_targets": len(vals), "mean_rho": sum(vals) / len(vals)}
    summary = {
        "status": "computed",
        "benchmark": VERSION,
        "assays_selected": len(panel),
        "independent_targets_analyzed": len(valid),
        "excluded_after_qc": len(assays) - len(valid),
        "primary_metric": "Equal-target mean Spearman correlation between per-site EvoTrace constraint and ProteinGym fraction-fit of single missense substitutions.",
        "direction_expected_from_method": "negative: more conserved positions tolerate a smaller fraction of missense substitutions",
        "mean_target_spearman_rho": mean_rho,
        "target_bootstrap_95_percentile_ci": _bootstrap(correlations),
        "negative_target_correlations": negative,
        "nonzero_target_correlations": nonzero,
        "one_sided_exact_sign_test_p": _sign_test_p(negative, nonzero),
        "category_summary": by_category,
        "total_evaluated_sites": sum(r["evaluated_sites"] for r in valid),
        "total_mapped_single_missense_variants": sum(r["mapped_single_missense_variants"] for r in valid),
        "inference_limit": "This evaluates the alignment-derived conservation/constraint component against external protein DMS assays. It does not validate codon selection tests, structural/domain annotations, or establish causal biological adaptation.",
    }
    protocol = {
        "benchmark": VERSION,
        "retrieved_on": date.today().isoformat(),
        "local_file_sha256": {
            "dms_data": _sha256(temp / "dms.zip"),
            "metadata": _sha256(temp / "metadata.csv"),
        },
        "urls": {"dms_data": DMS_ARCHIVE, "metadata": DMS_METADATA, "msa_archive": MSA_ARCHIVE},
        "panel_selection_without_outcomes": {
            "coarse_selection_types": list(CATEGORIES),
            "sequence_length_min": 50,
            "sequence_length_max": 350,
            "minimum_single_mutants": 1000,
            "minimum_msa_records": 100,
            "maximum_msa_records": 200000,
            "minimum_msa_coverage": 0.8,
            "targets_per_category": per_category,
            "selection": "sort eligible DMS IDs by SHA-256('20261001:'+DMS_id), take first targets/category while avoiding duplicate UniProt IDs",
        },
        "msa_subsample": {"homologs_per_target": sample_size, "sampling": "uniform reservoir sample over ProteinGym A2M records excluding the target row", "seed": "SHA-256('20261001:'+DMS_id) first 32 bits"},
        "a2m_conversion": "uppercase all amino-acid symbols; map A2M dot and dash gaps to '-' while retaining aligned columns",
        "coordinate_qc": "MSA target ungapped sequence must exactly match ProteinGym target_seq; single substitutions must match the reference residue; positions are projected through AlignmentMap",
        "site_aggregation": "mean DMS_score_bin over at least five single missense measurements at that residue",
        "primary_metric": "per-target Spearman rho(constraint_score, site mean DMS_score_bin), equally aggregated over independent UniProt targets",
        "uncertainty": "20,000 bootstrap resamples of targets; one-sided exact sign test for negative per-target correlations",
        "no_tuning": True,
        "panel": [
            {"dms_id": row["DMS_id"], "uniprot_id": row["UniProt_ID"], "category": row["coarse_selection_type"], "msa": row["MSA_filename"]}
            for row in panel
        ],
    }
    dms_path = output_dir / "per_target.csv"
    fields = [k for k, v in valid[0].items() if k != "sites"]
    with dms_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in assays:
            writer.writerow({key: row.get(key, "") for key in fields})
    site_path = output_dir / "per_site.csv"
    with site_path.open("w", newline="", encoding="utf-8") as stream:
        fields = ["dms_id", "uniprot_id", "selection_type"] + list(valid[0]["sites"][0])
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in valid:
            for site in row["sites"]:
                writer.writerow({"dms_id": row["dms_id"], "uniprot_id": row["uniprot_id"], "selection_type": row["selection_type"], **site})
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (output_dir / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n", encoding="utf-8")
    (output_dir / "excluded_targets.json").write_text(
        json.dumps([r for r in assays if r.get("status") == "excluded_after_qc"], indent=2) + "\n",
        encoding="utf-8",
    )
    (output_dir / "report.md").write_text(_report(summary, protocol, valid), encoding="utf-8")
    return summary


def _report(summary, protocol, valid):
    lines = [
        "# EvoTrace independent ProteinGym validation",
        "",
        "## Scope and protocol",
        "",
        "This external check tests EvoTrace's alignment-derived constraint summary against ProteinGym v1.3 experimental deep-mutational-scanning (DMS) results. It does not evaluate selection-model inference or annotation modules.",
        "",
        "The target panel was selected from metadata only, before reading measured mutation outcomes: 50–350 aa proteins, at least 1,000 single mutants, at least 100 MSA records, at least 80% MSA coverage, and up to {} unique targets from each of Activity, OrganismalFitness, and Stability. Each UniProt target appears at most once. The protocol and exact panel are recorded in `protocol.json`.".format(protocol["panel_selection_without_outcomes"]["targets_per_category"]),
        "",
        "ProteinGym documents `DMS_score_bin=1` as fit and 0 as not fit, with higher continuous DMS scores indicating higher fitness. For each target residue, at least five single missense outcomes were averaged; the primary metric is Spearman correlation between EvoTrace constraint and the per-site fraction-fit. The pre-specified expected direction is negative.",
        "",
        "## Result",
        "",
        "- Independent UniProt targets analyzed: {} ({} selected; {} excluded by sequence/coordinate/QC checks).".format(summary["independent_targets_analyzed"], summary["assays_selected"], summary["excluded_after_qc"]),
        "- Evaluated sites: {}; mapped single-missense outcomes: {}.".format(summary["total_evaluated_sites"], summary["total_mapped_single_missense_variants"]),
        "- Equal-target mean Spearman rho: {:.3f}; target-bootstrap 95% interval [{:.3f}, {:.3f}].".format(summary["mean_target_spearman_rho"], *summary["target_bootstrap_95_percentile_ci"]),
        "- Negative target correlations: {}/{}; one-sided exact sign-test p={:.4g}.".format(summary["negative_target_correlations"], summary["nonzero_target_correlations"], summary["one_sided_exact_sign_test_p"]),
        "",
        "| Assay category | Targets | Mean Spearman rho |",
        "|---|---:|---:|",
    ]
    for category, result in summary["category_summary"].items():
        lines.append("| {} | {} | {:.3f} |".format(category, result["n_targets"], result["mean_rho"]))
    lines += ["", "Per-target and per-site values are in `per_target.csv` and `per_site.csv`. Failed target QC is listed in `excluded_targets.json`.", "", "## Interpretation and limits", "", summary["inference_limit"], "The external target bootstrap and sign test quantify consistency across this selected target panel; they are not a substitute for broader prospective validation, replicate assays, or a peer-reviewed benchmark comparison.", "", "## Data provenance", "", "- Dataset release: ProteinGym v1.3.", "- DMS benchmark archive and references: [ProteinGym repository](https://github.com/OATML-Markslab/ProteinGym) and [release data](https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/).", "- ProteinGym reports 217 DMS substitution assays and describes target sequences and assay processing in its reference file. See the repository README and benchmark paper.", "- MSA member CRC32 values, DMS member CRC32 values, exact target list, filters, sampling seed, and coordinate rules are retained in `per_target.csv` and `protocol.json`.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="results/independent_validation", help="Compact validation artifacts output directory")
    parser.add_argument("--per-category", type=int, default=8)
    parser.add_argument("--sample-size", type=int, default=500)
    args = parser.parse_args()
    print(json.dumps(run_validation(args.output, args.per_category, args.sample_size), indent=2))


if __name__ == "__main__":
    main()
