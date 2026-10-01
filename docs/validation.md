# Independent validation

EvoTrace's alignment-derived constraint summary was compared with ProteinGym v1.3 experimental deep-mutational-scanning (DMS) substitution assays. The benchmark is external to EvoTrace: measured DMS outcomes were not used to fit or tune the method. This validates one component only; it does not validate HyPhy/PAML selection inference, domain/motif/structure integration, or biological adaptation claims.

## Observed result

The metadata-only panel selected 24 unique targets (up to 8 each in Activity, OrganismalFitness, and Stability). Five targets failed exact MSA target-sequence QC and were excluded. For each remaining target, at least five single-missense measurements were aggregated per residue. Across 19 targets, 2,690 sites, and 34,974 mapped variants, per-site EvoTrace constraint was negatively correlated with the experimental fraction-fit:

- Equal-target mean Spearman rho: **-0.433**.
- 95% percentile interval from 20,000 target-level bootstrap resamples: **[-0.515, -0.346]**.
- All 19 target-level correlations were negative; one-sided exact sign-test p = **1.91 × 10⁻⁶**.

The association is consistent with conserved residues tolerating fewer substitutions. It is not a causal test and does not show that EvoTrace predicts adaptation. The stability stratum had only three passing targets, so its estimate is particularly preliminary. Targets were selected by deterministic hash ordering from metadata under pre-specified coverage/depth/length filters; inference applies to this selected panel rather than all proteins or assay types.

## Reproduce

Install EvoTrace with its development dependencies, then run:

```bash
python scripts/validate_proteingym.py --output results/independent_validation
```

The script downloads the benchmark temporarily, performs deterministic reservoir sampling of up to 500 homologs per target, checks target sequence and residue mapping, and writes compact artifacts. Results include the exact panel and protocol, source URLs and SHA-256 checksums for the downloaded DMS archive and metadata, assay/MSA member CRC32 values, exclusions, per-target correlations, and per-site data. A full ProteinGym MSA archive is accessed with HTTP byte ranges and is not retained locally.

## Sources

- [ProteinGym repository and benchmark documentation](https://github.com/OATML-Markslab/ProteinGym)
- [ProteinGym v1.3 release data](https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/)
- [ProteinGym benchmark paper](https://papers.nips.cc/paper_files/paper/2023/file/cac723e5ff29f65e3fcbb0739ae91bee-Paper-Datasets_and_Benchmarks.pdf)

Machine-readable protocol and results are in `results/independent_validation/` (ignored by Git under the default results policy).
