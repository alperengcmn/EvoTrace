#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${EVOTRACE_PYTHON:-python3}"
export PATH="$(dirname "$PY"):$PATH"
cd "$ROOT"

echo "== Runtime and imports =="
"$PY" --version
"$PY" -m pip --version
"$PY" - <<'PY'
import importlib
for name in ("evotrace", "Bio", "numpy", "yaml", "pandas", "plotly", "streamlit", "pytest"):
    module = importlib.import_module(name)
    print(f"{name}: {getattr(module, '__version__', 'available')} ({getattr(module, '__file__', '')})")
PY
"$PY" -m evotrace --version
"$PY" -m evotrace doctor

echo "== Required executable smoke tests =="
for tool in mafft hmmscan hmmbuild hmmpress; do
    command -v "$tool" >/dev/null || { echo "Required tool missing: $tool" >&2; exit 1; }
done
MAFFT="$(command -v mafft)"
IQTREE="$(command -v iqtree2 || command -v iqtree || command -v iqtree3 || true)"
[[ -n "$IQTREE" ]] || { echo "Required tool missing: IQ-TREE (iqtree2/iqtree/iqtree3)" >&2; exit 1; }

TMP="$(mktemp -d "${TMPDIR:-/tmp}/evotrace-smoke.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
cat > "$TMP/dna.fasta" <<'EOF'
>taxon_a
ACGTTGCAACGTTCGATCGATGCTAGCTAGCTACGATCGTACGTTAGC
>taxon_b
ACGTTGCAACGTTCGATCGATGCTAGCTAGCTACGATCGTACGCTAGC
>taxon_c
ACGTTGCAACGTTCGATCGATGCTAGCTAACTACGATCGTACGCTAGC
>taxon_d
ACGTTGCAACGTTCGATCGATGCTAGCTAACTACGATCGTACGTTAGC
EOF
"$MAFFT" --auto --quiet "$TMP/dna.fasta" > "$TMP/aligned.fasta"
"$PY" - "$TMP/aligned.fasta" <<'PY'
import sys
from Bio import AlignIO
a = AlignIO.read(sys.argv[1], "fasta")
assert len(a) == 4 and a.get_alignment_length() > 0
print(f"MAFFT smoke: {len(a)} sequences, {a.get_alignment_length()} aligned sites")
PY
"$IQTREE" -s "$TMP/aligned.fasta" --seqtype DNA -m GTR -nt 1 --prefix "$TMP/iqtree" > "$TMP/iqtree.log" 2>&1
[[ -s "$TMP/iqtree.treefile" ]] || { cat "$TMP/iqtree.log"; echo "IQ-TREE produced no tree" >&2; exit 1; }
echo "IQ-TREE smoke: tree produced"

cat > "$TMP/proteins.faa" <<'EOF'
>p1
MKWVTFISLLFLFSSAYS
>p2
MKWVTFISLLFLFSSAYN
>p3
MKWVTFISLLFLFSSAYS
EOF
cat > "$TMP/proteins.sto" <<'EOF'
# STOCKHOLM 1.0
p1 MKWVTFISLLFLFSSAYS
p2 MKWVTFISLLFLFSSAYN
p3 MKWVTFISLLFLFSSAYS
//
EOF
hmmbuild "$TMP/test.hmm" "$TMP/proteins.sto" > "$TMP/hmmbuild.log" 2>&1
hmmpress -f "$TMP/test.hmm" > "$TMP/hmmpress.log" 2>&1
hmmscan --noali --tblout "$TMP/hits.tbl" "$TMP/test.hmm" "$TMP/proteins.faa" > "$TMP/hmmscan.log" 2>&1
grep -q '^proteins[[:space:]]' "$TMP/hits.tbl"
echo "HMMER smoke: profile built and hmmscan completed"

echo "== Selection executable smoke tests =="
for tool in hyphy codeml; do
    command -v "$tool" >/dev/null || { echo "Required selection tool missing: $tool" >&2; exit 1; }
done
cat > "$TMP/codons.fasta" <<'EOF'
>taxA
ATGGCTGCCGAAGTTTTTGCTGACCCTAAGGGTTATCAAGCTGCAAATCGTTGC
>taxB
ATGGCCGCCGAAGTTTTTGCTGACCCTAAGGGTTATCAAGCTGCAAATCGTTGC
>taxC
ATGGCTGCTGAAGTTTTTGATGACCCTAAGGGTTATCAAGCTGCAAATCGTTGC
>taxD
ATGGCGGCCGAGGTTTTTGCTGACCCTAAGGGTTATCAAGCTGCAAATCGTTGC
EOF
cat > "$TMP/codons.nwk" <<'EOF'
((taxA,taxB),(taxC,taxD));
EOF
"$PY" - "$TMP/codons.fasta" "$TMP/codons.nwk" "$TMP/hyphy" <<'PY'
import sys
from pathlib import Path
from evotrace.tools.hyphy import run_hyphy
r = run_hyphy("FEL", Path(sys.argv[1]), Path(sys.argv[3]), Path(sys.argv[2]))
assert r["result_level"] == "site" and r["sites"]
print(f"HyPhy FEL smoke: {len(r['sites'])} parsed site records ({r['version']})")
PY

echo "== Test suite =="
"$PY" -m pytest -q
echo "== Strict dependency check =="
"$PY" -m evotrace doctor --strict
echo "Environment verification complete."
