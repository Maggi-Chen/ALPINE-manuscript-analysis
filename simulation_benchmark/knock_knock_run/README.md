# Knock-knock benchmark run scripts

Scripts for running [knock-knock](https://github.com/jeffhussmann/knock-knock) v0.8.2 on the simulated PacBio HiFi and Oxford Nanopore datasets used in the ALPINE manuscript benchmark, plus utilities for parsing runtime logs and computing per-class precision/recall/F1 metrics.

These scripts produce the `outcome_list.txt` and `outcome_counts.csv` files that `../check_accuracy_knock.py` consumes for accuracy comparison against the ground truth.

## Pipeline overview

```
simulation FASTQ
      │
      ▼
trim_smrtbell.py  ────►  trimmed FASTQ
                              │
                              ▼
            bench_kk.sh / bench_kk_ont_only.sh
                              │  (knock-knock process-sample for HiFi;
                              │   run_nanopore.py for ONT)
                              ▼
                      outcome_list.txt
                      outcome_counts.csv
                              │
                              ├──►  parse_timing.py  ──►  benchmark_summary.csv
                              │
                              └──►  compute_benchmark_metrics.py
                                            │
                                            ▼
                                  per_class_metrics.csv
                                  (precision, recall, F1)
                                            │
                                            ▼
                              ../check_accuracy_knock.py
                                  (final accuracy comparison)
```

## Files

| Script | Purpose |
|---|---|
| `trim_smrtbell.py` | Exhaustive k-mer match SMRTbell adapter trimmer. CLI: `--input FASTQ.gz --output FASTQ.gz --min-k INT`. Use `--min-k 13` for PacBio HiFi, `--min-k 11` for ONT. |
| `bench_kk.sh` | 4-quadrant benchmark wrapper. Runs HiFi via knock-knock CLI (`process-sample`) and ONT via `run_nanopore.py`. Wraps each call in GNU `/usr/bin/time -v` for runtime / RAM / CPU capture. |
| `bench_kk_ont_only.sh` | ONT-only resume script (used when the HiFi half had already completed). |
| `run_nanopore.py` | Python driver that instantiates `NanoporeExperiment` directly. Required because knock-knock's CLI does not dispatch the nanopore platform (its `experiment.get_exp_class()` only wires illumina/pacbio). The driver also patches `HDR.Architecture.__init__` to set `max_indel_allowed_in_donor=5` for nanopore, which knock-knock 0.8.2 does not initialize for that platform. |
| `parse_timing.py` | Parses GNU `/usr/bin/time -v` reports under `timing_logs/` into `benchmark_summary.csv`. Handles the `Elapsed (wall clock) time (h:mm:ss or m:ss): VALUE` line whose label contains internal colons. |
| `compute_benchmark_metrics.py` | Computes per-class TP/FP/FN/precision/recall/F1 from `outcome_list.txt`, mapping knock-knock's hierarchical `(category, subcategory)` labels onto ALPINE's 13-class schema. CLI: `--base-dir <run_dir>`. Truth labels are parsed from PBSIM3 read names (`simulated_<gene>_<class>.<index>`). |

## Why the SMRTbell trimming step

Preliminary runs of knock-knock on the raw simulation FASTQs placed 11.7–34.5% of reads in a single `malformed layout: extra copy of primer` category, with a 3× TRBC vs. TRAC asymmetry that was platform-flat across HiFi and ONT. This asymmetry was traced to incidental sequence similarity between the 24-nt PacBio SMRTbell adapter motif and the interior of the TRBC forward primer, which passed knock-knock's ≤5-edit threshold in `realign_edge_to_primer` (`knock_knock/architecture/__init__.py:1582`) and triggered `len(als) > 1` at `HDR.py:600`. Trimming SMRTbell adapter motifs from read ends with `trim_smrtbell.py` reduces the `extra_copy_of_primer` rate by 7×–241× and recovers knock-knock's per-class F1 on classes it can natively express (HDR, simple indels, deletions) to within a few points of ALPINE on the same input. Full forensic chain documenting this finding is in the manuscript supplementary methods and in the Zenodo deposit's `docs/METHODS.md`.

## Conda environment

The benchmark runs were executed in a Python 3.14 conda environment (knock-knock requires Python ≥ 3.12):

```bash
conda create -n py313 python=3.14
conda activate py313
pip install knock-knock==0.8.2
conda install -c conda-forge time blast minimap2 samtools
```

`compute_benchmark_metrics.py` and `parse_timing.py` are pure Python (stdlib only) and have no external dependencies beyond Python ≥ 3.8.

## To reproduce the manuscript benchmark

The complete benchmark archive — including the simulation FASTQs (trimmed and untrimmed), per-quadrant `outcome_list.txt` and `outcome_counts.csv` from knock-knock, sample sheets, references, timing logs, and the merged per-class metrics workbook — is deposited on Zenodo (DOI: see manuscript Data Availability statement). To reproduce the four-quadrant benchmark from scratch:

1. Pull the original simulation FASTQs from the Zenodo deposit (`simulation_<platform>_<gene>.fastq.gz`, four files).
2. Trim SMRTbell adapter motifs:
   ```bash
   python trim_smrtbell.py --input simulation_hifi_trbc.fastq.gz --output simulation_hifi_trbc_trim.fastq.gz --min-k 13
   ```
   Use `--min-k 11` for ONT inputs.
3. Build knock-knock strategies from the references (`samtools faidx` on each FASTA first), then run the four quadrants:
   ```bash
   ./bench_kk.sh
   ```
4. Compute per-class metrics and parse timing:
   ```bash
   python parse_timing.py timing_logs/ > benchmark_summary.csv
   python compute_benchmark_metrics.py --base-dir results/
   ```
5. Compare against the ground-truth labels using `../check_accuracy_knock.py` (existing in this repository).

## Reference

These scripts were used to generate the knock-knock benchmark numbers reported in the ALPINE manuscript revision response to Reviewer 2 Comment #4 (head-to-head with knock-knock on simulated data).
