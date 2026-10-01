import tempfile
import unittest
import shutil
import subprocess
from pathlib import Path
from evotrace.io.fasta import read_fasta, FastaError
from evotrace.validation.fasta import validate_records, ValidationError
from evotrace.models import SequenceRecord
from evotrace.conservation.core import columns, alignment_stats
from evotrace.selection.ng import pairwise_dnds
from evotrace.scoring.ets import evolutionary_scores
from evotrace.conservation.information import shannon_entropy, jensen_shannon_divergence
from evotrace.mapping.residues import AlignmentMap
from evotrace.codon.backtranslate import backtranslate, CodonAlignmentError, validate_codon_alignment
from evotrace.evidence.sites import candidate_sites
from evotrace.pipeline.engine import PipelineEngine
from evotrace.tools.hmmer import parse_domtblout
from evotrace.tools.hyphy import parse_hyphy
from evotrace.tools.paml import parse_codeml_output, _likelihood_ratio, run_codeml
from evotrace.reporting.html import render_report
from evotrace.structure.providers import map_sequence_to_structure
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from evotrace.pipeline.analyze import analyze


class CoreTests(unittest.TestCase):
    def test_fasta_and_duplicate(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.fa"
            p.write_text(">a\nACGT\n>b\nACGA\n")
            self.assertEqual(len(read_fasta(p)), 2)
            p.write_text(">a\nACGT\n>a again\nACGA\n")
            with self.assertRaises(FastaError):
                read_fasta(p)

    def test_validation_alphabet(self):
        kind, w = validate_records([SequenceRecord("a", "ACGT"), SequenceRecord("b", "ACGA")])
        self.assertEqual(kind, "dna")
        with self.assertRaises(ValidationError):
            validate_records([SequenceRecord("a", "AC!T")], "dna")

    def test_entropy_conservation(self):
        a = [SequenceRecord("a", "AA-"), SequenceRecord("b", "AT-")]
        c = columns(a)
        self.assertEqual(c[0]["entropy"], 0)
        self.assertEqual(c[1]["entropy"], 1)
        self.assertEqual(c[2]["gap_fraction"], 1)
        self.assertEqual(alignment_stats(a)["alignment_length"], 3)

    def test_descriptive_alignment_summaries(self):
        x = evolutionary_scores(columns([SequenceRecord("a", "AA"), SequenceRecord("b", "AT")]))
        self.assertFalse(any("adaptation" in key for key in x[0]))
        self.assertIn("divergence_score", x[0])

    def test_dnds_requires_codon_frame(self):
        x = pairwise_dnds([SequenceRecord("a", "ATG"), SequenceRecord("b", "ATA")])
        self.assertEqual(x[0]["status"], "descriptive_only")
        x = pairwise_dnds([SequenceRecord("a", "AT"), SequenceRecord("b", "AA")])
        self.assertEqual(x[0]["status"], "not_run")

    def test_gap_policy_and_information_metrics(self):
        self.assertEqual(shannon_entropy("AAAA")["entropy"], 0)
        self.assertAlmostEqual(shannon_entropy("AATT")["entropy"], 1.0)
        self.assertEqual(shannon_entropy("-A", "ignore")["entropy"], 0)
        self.assertAlmostEqual(shannon_entropy("-A", "include", alphabet_size=4)["entropy"], 1.0)
        dna_gap = shannon_entropy("-A", "include", alphabet_size=4, valid_symbols=set("ACGT"))
        self.assertLessEqual(dna_gap["normalized_entropy"], 1.0)
        protein_ambiguous = shannon_entropy(
            "AX", "ignore", alphabet_size=20, valid_symbols=set("ACDEFGHIKLMNPQRSTVWY")
        )
        self.assertEqual(protein_ambiguous["entropy"], 0.0)
        self.assertAlmostEqual(jensen_shannon_divergence({"A": 1}, {"B": 1}), 1.0)
        self.assertAlmostEqual(jensen_shannon_divergence({"A": 1}, {"A": 1}), 0.0)

    def test_gap_aware_coordinate_round_trip(self):
        aln = [SequenceRecord("ref", "A-CGT"), SequenceRecord("q", "ATCG-")]
        m = AlignmentMap(aln, "ref")
        self.assertEqual(m.sequence_to_alignment("ref", 2), 3)
        self.assertIsNone(m.alignment_to_sequence("ref", 2))
        self.assertEqual(m.alignment_to_reference(3), 2)
        self.assertIsNone(m.coordinate("q", 5).sequence_position)
        self.assertEqual(m.coordinate("q", 4).sequence_position, 4)
        for sequence_id, position in (("missing", 1), ("ref", 0), ("ref", 6)):
            with self.assertRaises(ValueError):
                m.coordinate(sequence_id, position)

    def test_codon_backtranslation_and_frame_errors(self):
        protein = [SequenceRecord("x", "M-A"), SequenceRecord("y", "M-A")]
        cds = [SequenceRecord("x", "ATGGCT"), SequenceRecord("y", "ATGGCT")]
        projected, meta = backtranslate(protein, cds)
        self.assertEqual(projected[0].sequence, "ATG---GCT")
        self.assertEqual(meta["sequences"]["x"]["codons_used"], 2)
        with self.assertRaises(CodonAlignmentError):
            backtranslate([SequenceRecord("x", "MA")], [SequenceRecord("x", "ATGA")])
        with self.assertRaises(CodonAlignmentError):
            backtranslate(
                [SequenceRecord("x", "MA")],
                [SequenceRecord("x", "ATGGCT"), SequenceRecord("unused", "ATGGCT")],
            )

    def test_codon_alignment_matches_protein_gap_pattern_and_translation(self):
        protein = [SequenceRecord("x", "M-A"), SequenceRecord("y", "M-A")]
        codons = [SequenceRecord("x", "ATG---GCT"), SequenceRecord("y", "ATG---GCC")]
        report = validate_codon_alignment(codons, protein)
        self.assertTrue(report["protein_correspondence_checked"])
        self.assertEqual(report["gap_containing_codons"], 2)
        with self.assertRaises(CodonAlignmentError):
            validate_codon_alignment(
                [SequenceRecord("x", "ATGGCT"), SequenceRecord("y", "ATGGCC")],
                [SequenceRecord("x", "MA"), SequenceRecord("y", "MM")],
            )

    def test_candidate_site_reference_positions(self):
        aln = [SequenceRecord("ref", "A-C"), SequenceRecord("b", "ATC")]
        cs = columns(aln, "protein")
        sites = candidate_sites(aln, cs, evolutionary_scores(cs), reference_id="ref")
        self.assertEqual([x["reference_position"] for x in sites], [1, None, 2])

    def test_candidate_site_integrates_motif_domain_selection_structure(self):
        aln = [SequenceRecord("ref", "A-C"), SequenceRecord("b", "ATC")]
        cs = columns(aln, "protein")
        motifs = [
            {
                "motif": "test_motif",
                "sequence_id": "ref",
                "start": 1,
                "end": 1,
                "alignment_position": 1,
                "alignment_end": 1,
            }
        ]
        domains = [
            {
                "sequence_id": "ref",
                "start": 1,
                "end": 1,
                "domain_name": "D1",
                "score": 30.0,
                "evalue": 1e-5,
            }
        ]
        selection = [
            {
                "alignment_position": 1,
                "method": "FEL",
                "selection_class": "negative/purifying",
                "p_value": 0.01,
                "q_value": 0.02,
            }
        ]
        structure = {
            "status": "mapping_successful",
            "structure_id": "demo",
            "mapped_residues": [
                {
                    "sequence_position": 1,
                    "structure_chain": "A",
                    "structure_residue_number": 4,
                    "structure_residue": "A",
                    "mean_b_factor": 91.0,
                }
            ],
        }
        sites = candidate_sites(
            aln, cs, evolutionary_scores(cs), "ref", motifs, domains, selection, structure
        )
        row = sites[0]
        self.assertEqual(row["domain"], "D1")
        self.assertEqual(row["motif"], "test_motif")
        self.assertEqual(row["structure_residue"], "A:4")
        self.assertEqual(row["domain_evalue"], 1e-5)
        evidence = __import__("json").loads(row["evidence_summary"])
        self.assertTrue(
            {"conservation", "negative/purifying", "motif", "domain", "structure"}.issubset(
                {item["evidence_type"] for item in evidence}
            )
        )

    def test_candidate_projection_from_query_coordinates_and_missing_statuses(self):
        aln = [SequenceRecord("ref", "A-CGT"), SequenceRecord("query", "ATCG-")]
        cs = columns(aln, "protein")
        rows = candidate_sites(
            aln, cs, evolutionary_scores(cs), "ref",
            domains=[{"sequence_id": "query", "start": 2, "end": 4, "domain_name": "D"}],
            selection=[{"method": "aBSREL", "result_level": "branch", "site": 2}],
            modality_statuses={"domains": "available", "motifs": "not_configured", "selection": "available", "structure": "not_configured"},
        )
        self.assertEqual(rows[2]["domain"], "D")  # query residue 2 -> alignment column 3
        self.assertEqual(rows[2]["domain_reference_start"], 2)
        self.assertEqual(rows[1]["selection_status"], "no_hit")
        self.assertEqual(rows[1]["motif_status"], "not_configured")
        self.assertFalse(any("adaptation" in key for key in rows[1]))

    def test_structure_mapping_preserves_insertion_codes(self):
        with tempfile.TemporaryDirectory() as d:
            pdb = Path(d) / "insertions.pdb"
            rows = []
            for i, (resname, number, insertion) in enumerate(
                (("ALA", 10, " "), ("CYS", 10, "A"), ("ASP", 12, " ")), 1
            ):
                rows.append(
                    "ATOM   {:4d}  CA  {} A{:4d}{}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}           C\n".format(
                        i, resname, number, insertion, float(i), 0.0, 0.0, 1.0, 80.0
                    )
                )
            pdb.write_text("".join(rows) + "TER\nEND\n")
            result = map_sequence_to_structure("ACD", pdb, "A")
            self.assertEqual(result["status"], "mapping_successful")
            self.assertEqual(result["mapped_residues"][1]["structure_insertion_code"], "A")

    def test_pipeline_cache_records_stages(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            artifact = root / "stage.txt"
            engine = PipelineEngine(root)
            engine.run_stage(
                "one", lambda: (artifact.write_text("ok"), None)[1], {"input": "abc"}, [artifact]
            )
            reused = PipelineEngine(root, resume=True)
            result = reused.run_stage(
                "one", lambda: self.fail("cache miss"), {"input": "abc"}, [artifact]
            )
            self.assertIsNone(result)
            self.assertTrue(reused.state["stages"]["one"]["cache_hit"])
            artifact.write_text("tampered")
            refreshed = PipelineEngine(root, resume=True)
            refreshed.run_stage(
                "one", lambda: artifact.write_text("recomputed"), {"input": "abc"}, [artifact]
            )
            self.assertEqual(artifact.read_text(), "recomputed")

    def test_hmmer_domtblout_parser(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.domtbl"
            p.write_text(
                "PF00001 PF00001.1 100 queryA - 120 1e-20 50 0 1 1 1e-10 1e-10 45 0 2 30 10 38 8 40 0.9 kinase domain\n"
            )
            hit = parse_domtblout(p)[0]
            self.assertEqual(hit["domain_id"], "PF00001")
            self.assertEqual((hit["start"], hit["end"]), (10, 38))

    def test_hyphy_method_specific_parsers(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "hyphy.json"
            common = {
                "MLE": {
                    "headers": [
                        ["alpha", "Synonymous rate"],
                        ["beta", "Non-synonymous rate"],
                        ["LRT", "Likelihood ratio"],
                        ["p-value", "Asymptotic p-value"],
                    ],
                    "content": {"0": [[0.2, 1.4, 5.0, 0.01], [0.8, 0.2, 0.1, 0.8]]},
                },
                "data partitions": {"0": {"coverage": [[0, 5]]}},
            }
            path.write_text(__import__("json").dumps(common))
            fel = parse_hyphy(path, "FEL")["sites"]
            self.assertEqual([x["site"] for x in fel], [1, 6])
            self.assertEqual(fel[0]["selection_class"], "positive/diversifying_BH_q<=0.05")
            meme = parse_hyphy(path, "MEME")["sites"]
            self.assertEqual(meme[0]["selection_class"], "episodic_positive_candidate_BH_q<=0.05")
            fubar_input = {
                "grid": [[0.2, 1.0, 0.5], [2.0, 1.0, 0.5]],
                "posterior": {"0": [[0.05, 0.95]]},
                "data partitions": {"0": {"coverage": [11]}},
            }
            path.write_text(__import__("json").dumps(fubar_input))
            fubar = parse_hyphy(path, "FUBAR")["sites"][0]
            self.assertAlmostEqual(fubar["posterior_positive"], 0.95)
            self.assertIsNone(fubar["p_value"])
            absrel_input = {
                "branch attributes": {
                    "0": {
                        "branchA": {
                            "Uncorrected P-value": 0.01,
                            "Corrected P-value": 0.03,
                            "LRT": 4.2,
                            "Rate classes": 2,
                            "Rate Distributions": [[2.0, 1.0]],
                        }
                    },
                    "attributes": {},
                }
            }
            path.write_text(__import__("json").dumps(absrel_input))
            branch = parse_hyphy(path, "aBSREL")["branches"][0]
            self.assertEqual(branch["level"], "branch")
            self.assertEqual(branch["adjusted_p_value"], 0.03)
            busted_input = {"test results": {"LRT": 7.0, "p-value": 0.004}, "fits": {}}
            path.write_text(__import__("json").dumps(busted_input))
            busted = parse_hyphy(path, "BUSTED")["gene_tests"][0]
            self.assertEqual(busted["level"], "gene_branch_set")
            self.assertEqual(busted["p_value"], 0.004)

    def test_codeml_labeled_parser_keeps_posterior_not_pvalue(self):
        text = """lnL(ntime: 8 np: 9): -123.456
omega (dN/dS) = 1.250
Bayes Empirical Bayes (BEB) analysis
  87 A 0.973** 3.41
"""
        result = parse_codeml_output(text, "branch_site")
        self.assertEqual(result["log_likelihood"], -123.456)
        self.assertEqual(result["omega_estimates"], [1.25])
        self.assertEqual(result["site_posterior_results"][0]["site"], 87)
        self.assertNotIn("p_value", result["site_posterior_results"][0])
        self.assertEqual(result["site_posterior_results"][0]["result_level"], "site")
        self.assertEqual(result["result_level"], "gene")
        branch_site_lrt = _likelihood_ratio(-100, -98, "branch_site")
        self.assertAlmostEqual(branch_site_lrt["statistic"], 4)
        self.assertAlmostEqual(branch_site_lrt["p_value"], 0.5 * __import__("math").erfc(2 ** 0.5))
        self.assertEqual(_likelihood_ratio(-100, -99, "site")["degrees_of_freedom"], 2)

    def test_codeml_standard_code_and_stop_codon_validation(self):
        from evotrace.tools.paml import _control_text, _prepare_codeml_inputs
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            alignment = root / "codons.fa"
            tree = root / "tree.nwk"
            alignment.write_text(">a\nATGTAA\n>b\nATGTAA\n>c\nATGTAA\n")
            tree.write_text("((a,b),c);\n")
            with self.assertRaisesRegex(ValueError, "stop codon"):
                _prepare_codeml_inputs(alignment, tree, root / "out")
        control = _control_text("a.phy", "a.tree", "out", 2, 2, 0, 1)
        self.assertIn("icode = 0", control)

    def test_report_shows_paml_lrt_without_conflating_beb_posterior(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "report.html"
            render_report({}, [], {}, [], [{
                "method": "PAML/codeml", "result_level": "gene", "status": "ok",
                "site_posterior_results": [{"posterior_probability": 0.98}],
                "lrt": 3.2, "lrt_p_value": 0.04, "lrt_degrees_of_freedom": 1,
                "lrt_status": "computed",
            }], out)
            report = out.read_text()
            self.assertIn("LRT 3.2 (df 1, p=0.04)", report)
            self.assertIn("posterior probabilities retain their own meanings", report)

    @unittest.skipUnless(shutil.which("codeml"), "PAML codeml is required for real model smoke test")
    def test_real_codeml_site_branch_and_branch_site_models(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            alignment = root / "codons.fasta"
            tree = root / "tree.nwk"
            alignment.write_text(
                ">a\nATGGCTTTTGAACCTGATCAAGGTCGTAAACGT\n"
                ">b\nATGGCCTTTGAACCCGATCAGGGTCGTAAGCGT\n"
                ">c\nATGGCTTTCGAGCCTGACCAAGGCCGTAAACGC\n"
                ">d\nATGGCCTTCGAGCCCGACCAGGGCCGCAAGCGC\n"
            )
            tree.write_text("((a,b),(c,d)#1);\n")
            for model in ("site", "branch", "branch_site"):
                result = run_codeml(alignment, tree, root / model, model=model)
                self.assertEqual(result["lrt_status"], "computed", model)
                self.assertIsNotNone(result["lrt_p_value"], model)

    @unittest.skipUnless(
        all(shutil.which(x) for x in ("mafft", "hmmbuild", "hmmpress", "hmmscan", "hyphy")),
        "MAFFT, HMMER and HyPhy are required for integrated modality workflow",
    )
    def test_real_multimodality_candidate_site_integration(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            codons = [
                "ATGGCTGCCGAAGTTTTTGCTGACCCTAAGGGTTATCAAGCTGCAAATCGTTGC" * 4,
                "ATGGCCGCCGAAGTTTTTGCTGACCCTAAGGGTTATCAAGCTGCAAATCGTTGC" * 4,
                "ATGGCTGCTGAAGTTTTTGATGACCCTAAGGGTTATCAAGCTGCAAATCGTTGC" * 4,
                "ATGGCGGCCGAGGTTTTTGCTGACCCTAAGGGTTATCAAGCTGCAAATCGTTGC" * 4,
            ]
            ids = ("taxA", "taxB", "taxC", "taxD")
            cds = root / "cds.fasta"
            proteins = root / "proteins.fasta"
            cds.write_text("".join(f">{name}\n{seq}\n" for name, seq in zip(ids, codons)))
            protein_records = [(name, str(Seq(seq).translate())) for name, seq in zip(ids, codons)]
            proteins.write_text("".join(f">{name}\n{seq}\n" for name, seq in protein_records))
            profile = root / "profile.sto"
            profile.write_text(
                "# STOCKHOLM 1.0\n" + "".join(f"{name} {seq}\n" for name, seq in protein_records) + "//\n"
            )
            hmm = root / "profile.hmm"
            subprocess.run([shutil.which("hmmbuild"), str(hmm), str(profile)], check=True, capture_output=True)
            subprocess.run([shutil.which("hmmpress"), "-f", str(hmm)], check=True, capture_output=True)
            pdb = root / "structure.pdb"
            pdb_lines = []
            for i, aa in enumerate(protein_records[0][1], 1):
                pdb_lines.append(
                    "ATOM   {:4d}  CA  {:3s} A{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}           C\n".format(
                        i, seq3(aa).upper(), i, float(i), 0.0, 0.0, 1.0, 90.0
                    )
                )
            pdb.write_text("".join(pdb_lines) + "TER\nEND\n")
            config = root / "config.yaml"
            config.write_text(
                "input:\n  alphabet: protein\n"
                "phylogeny:\n  tool: nj\n"
                "selection:\n  enabled: true\n  engine: hyphy\n  analyses: [FEL]\n  coding_sequences: cds.fasta\n"
                "domains:\n  enabled: true\n  database: profile.hmm\n"
                "structure:\n  enabled: true\n  pdb_path: structure.pdb\n  chain: A\n"
                "motifs:\n  custom:\n    - name: start_pair\n      pattern: 'M.'\n"
            )
            output = root / "out"
            summary = analyze(str(proteins), str(output), str(config))
            self.assertEqual(summary["stage_status"]["selection"], "SUCCESS")
            rows = __import__("json").loads((output / "candidate_sites.json").read_text())
            self.assertTrue(any(r["selection_status"] == "hit" for r in rows))
            self.assertTrue(any(r["domain_status"] == "hit" for r in rows))
            self.assertTrue(any(r["motif_status"] == "hit" for r in rows))
            self.assertTrue(any(r["structure_status"] == "mapped" for r in rows))
            self.assertFalse(any("adaptation" in key for key in rows[0]))
            integration = __import__("json").loads((output / "evidence_integration.json").read_text())
            self.assertEqual(integration["status"], "complete")

    @unittest.skipUnless(
        shutil.which("mafft"), "MAFFT is required for the real workflow integration test"
    )
    def test_real_example_end_to_end(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "analysis"
            summary = analyze(config_path=root / "configs" / "example.yaml", output_dir=output)
            self.assertEqual(summary["stage_status"]["alignment_qc"], "SUCCESS")
            self.assertTrue((output / "alignment.fasta").is_file())
            self.assertTrue((output / "candidate_sites.csv").is_file())
            report = (output / "report.html").read_text(encoding="utf-8")
            self.assertIn("tp53_mammals.fasta", report)
            self.assertTrue((output / "provenance.json").is_file())
            resumed = analyze(
                config_path=root / "configs" / "example.yaml", output_dir=output, resume=True
            )
            stages = __import__("json").loads((output / "pipeline_stages.json").read_text())[
                "stages"
            ]
            self.assertTrue(stages["alignment"]["cache_hit"])
            self.assertEqual(resumed["stage_status"]["alignment_qc"], "SUCCESS")


if __name__ == "__main__":
    unittest.main()
