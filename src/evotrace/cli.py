"""Command-line interface for EvoTrace."""

import argparse
import logging
import sys
from pathlib import Path

from evotrace import __version__
from evotrace.doctor import print_doctor
from evotrace.io.fasta import read_fasta
from evotrace.pipeline.analyze import analyze
from evotrace.tools.hmmer import run_hmmscan
from evotrace.tools.hyphy import run_hyphy
from evotrace.tools.iqtree import run_iqtree
from evotrace.tools.paml import run_codeml
from evotrace.validation.fasta import validate_records


def parser():
    p = argparse.ArgumentParser(
        prog="evotrace", description="Reproducible evolutionary sequence analysis"
    )
    p.add_argument("--version", action="version", version="EvoTrace {}".format(__version__))
    sub = p.add_subparsers(dest="command", required=True)
    for name in ("analyze", "pipeline"):
        q = sub.add_parser(name, help="Run the integrated analysis workflow")
        inputs = q.add_mutually_exclusive_group()
        inputs.add_argument("--input", "-i", required=False, help="FASTA input")
        inputs.add_argument("--accession", nargs="+", help="Fetch NCBI protein accession(s)")
        q.add_argument(
            "--output", "-o", default="results", help="Output directory (default: results)"
        )
        q.add_argument("--config", "-c", help="JSON/YAML configuration")
        q.add_argument("--threads", "-t", type=int, default=1)
        q.add_argument(
            "--workers",
            type=int,
            default=1,
            help="Concurrent independent analyses (reserved; currently serial)",
        )
        q.add_argument("--dry-run", action="store_true")
        q.add_argument(
            "--resume", action="store_true", help="Reuse successful content-matched stages"
        )
        q.add_argument("--verbose", action="store_true")
    v = sub.add_parser("validate", help="Validate FASTA input")
    v.add_argument("--input", "-i", required=True)
    v.add_argument("--alphabet", choices=("auto", "dna", "protein"), default="auto")
    phy = sub.add_parser("phylogeny", help="Run a maximum-likelihood IQ-TREE analysis")
    phy.add_argument("--alignment", required=True)
    phy.add_argument("--output", default="results/phylogeny")
    phy.add_argument("--model", default="MFP")
    phy.add_argument("--bootstrap", type=int, default=1000)
    phy.add_argument("--threads", type=int, default=1)
    sel = sub.add_parser("selection", help="Run an explicitly selected HyPhy or PAML analysis")
    sel.add_argument("--alignment", required=True)
    sel.add_argument("--output", default="results/selection")
    sel.add_argument("--engine", choices=("hyphy", "paml"), default="hyphy")
    sel.add_argument(
        "--method", choices=("FEL", "MEME", "FUBAR", "aBSREL", "BUSTED"), default="FEL"
    )
    sel.add_argument("--tree")
    sel.add_argument("--model", choices=("branch", "site", "branch_site"), default="branch_site")
    dom = sub.add_parser("domains", help="Run HMMER hmmscan against a local profile database")
    dom.add_argument("--sequences", required=True)
    dom.add_argument("--database", required=True)
    dom.add_argument("--output", default="results/domains")
    dom.add_argument("--threads", type=int, default=1)
    doctor = sub.add_parser("doctor", help="Check installed executables and Python dependencies")
    doctor.add_argument("--strict", action="store_true", help="Return failure if a required dependency is missing")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if getattr(args, "verbose", False) else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    try:
        if args.command == "doctor":
            return print_doctor(strict=args.strict)
        if args.command == "validate":
            rec = read_fasta(args.input)
            alphabet, warnings = validate_records(rec, args.alphabet)
            print("Valid {} FASTA: {} sequences".format(alphabet, len(rec)))
            for warning in warnings:
                print("WARNING: {}".format(warning))
            return 0
        if args.command in ("analyze", "pipeline"):
            if not args.input and not args.accession and not args.config:
                raise ValueError("Provide --input, --accession, or --config with an input source.")
            analyze(
                args.input,
                args.output,
                args.config,
                args.dry_run,
                args.threads,
                args.resume,
                workers=args.workers,
                accessions=args.accession,
            )
            return 0
        if args.command == "phylogeny":
            result = run_iqtree(
                Path(args.alignment), Path(args.output), args.model, args.bootstrap, args.threads
            )
            print("IQ-TREE complete: {}".format(result["tree_path"]))
            print(
                "Model: {}; log-likelihood: {}".format(
                    result["selected_model"], result["log_likelihood"]
                )
            )
            return 0
        if args.command == "selection":
            if args.engine == "hyphy":
                result = run_hyphy(
                    args.method,
                    Path(args.alignment),
                    Path(args.output),
                    Path(args.tree) if args.tree else None,
                )
            else:
                if not args.tree:
                    raise ValueError("PAML/codeml requires --tree.")
                result = run_codeml(
                    Path(args.alignment), Path(args.tree), Path(args.output), args.model
                )
            print(
                "{} completed: {}".format(
                    args.engine, result.get("path", result.get("output_file"))
                )
            )
            return 0
        if args.command == "domains":
            result = run_hmmscan(
                Path(args.sequences), Path(args.database), Path(args.output), args.threads
            )
            print("hmmscan completed: {} domain hits".format(len(result["hits"])))
            return 0
        return 2
    except Exception as exc:
        logging.error("%s", exc)
        if getattr(args, "verbose", False):
            logging.exception("Detailed failure")
        return 2


if __name__ == "__main__":
    sys.exit(main())
