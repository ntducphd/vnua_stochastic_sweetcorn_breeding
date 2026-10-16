# Stochastic sweet corn breeding

Simulation engine, founder populations and outputs for the paper

> Nguyen Trung Duc, Nguyen Thi Hoai Linh, Ta Viet Ton. Stochastic sweet corn breeding: a reproducible
> simulation workflow for programme design. Manuscript in preparation.

## Running

From the repository root:

| Command | Output | Time |
|---|---|---|
| `bash run_all.sh` | every table, figure and quoted number, from the stored outputs | about 1 minute |
| `bash run_all.sh --full` | the engine tests, then the application rerun into `out_full/` from the stored founder objects and compared record for record with `out/`; Experiments 1 and 2 rerun with new coalescent draws | about 13 CPU-hours |
| `Rscript sim/tests/run_tests.R` | the engine tests (`sim/tests/test_results.txt`) | about 3 minutes |

`run_all.sh` ends with three checks: the rerun of replicate 4001 (`ALL RECORDS IDENTICAL`), every number quoted
in the paper (`scripts/verify_numbers.py`, `ALL NUMBERS REPRODUCE`) and every regenerated CSV, JSON and PNG file
against `checksums/outputs.sha256`, ignoring line endings. PNG files match only with the Arial font and cairo
build of the deposited run; PDF figures carry a creation date and are not checksummed.

## Requirements

R 4.5.1 with AlphaSimR 2.1.0, rrBLUP 4.6.3 and jsonlite; Python 3.13 with the packages in `requirements.txt`.
All runs use one thread.

    Rscript -e 'options(repos = "https://cloud.r-project.org"); install.packages("remotes"); remotes::install_version("AlphaSimR", "2.1.0"); remotes::install_version("rrBLUP", "4.6.3"); install.packages("jsonlite")'
    pip install -r requirements.txt

AlphaSimR 2.1.0 is required. Its `runMacs` does not return the same founder genome twice under one seed
(Experiment 1), and the fix in the development version changes the order of the random draws, so the founder
objects in `out/founders/`, not the seeds, define the application. A saved founder object is reloaded, with its
seed, configuration and SHA-256 checksum checked. `run_all.sh` stops under another version unless
`ALLOW_VERSION=1` is set.

## Contents

| Path | Contents | In the paper |
|---|---|---|
| `sim/R/` | the simulation engine: configuration (Table 1), founder objects, training archive and ridge prediction, programmes as strategy objects, the cycle of Algorithm 1, outcomes, input and output | Sections 2.1, 2.2; Table 1; Algorithm 1 |
| `sim/run_application.R`, `sim/run_exp1.R`, `sim/run_exp2.R` | drivers for the application (all scenarios, programmes and count-matched controls of one replicate per call) and Experiments 1 and 2; each skips finished work | Sections 2.3, 2.4, 2.6 |
| `sim/tests/` | engine tests: Mendelian ratios, count-matched replay identities, shared start of paired arms, visual-sorting identity, reproducibility from a saved founder object in a new session, plot accounting | Section 3.1 |
| `sim/SPEC_NOTES.md` | implementation choices, each with the sentence of the paper it implements | |
| `R/01_pool_application.py` to `R/04_analyse_matched_controls.R` | pooled records, analyses, prediction bias in t ha⁻¹, decomposition | Sections 3.3, 3.4; Tables 4–6 |
| `R/05_bootstrap.py`, `R/06_inference.py` | bootstrap intervals; inference and `out/numbers.json` | Sections 3.2–3.5; Tables 3, 4, 6 |
| `R/07_figures_1_2.R`, `R/08_figures_3_5.R`, `R/fig_style.R` | figures | Figures 1–5 |
| `R/exp2_pool.py`, `R/exp2_summary.py` | Experiment 2 tables | Table 3 |
| `scripts/` | number check, rerun comparison, extraction of the published tables, session record | |
| `out/founders/` | the ten founder objects of the application (seeds 4001–4010) with checksums | |
| `out/application/` | per-cycle records of the five programmes in the four scenarios (`raw_w<seed>.csv`), the pooled records and their summaries | Sections 3.3, 3.4 |
| `out/matched_controls/` | per-cycle records of the count-matched controls (`raw_m<seed>.csv`) and their analysis; `estimator_distinct.csv` (Table 4) and `bootstrap_intervals.csv` (Table 4, Fig. 1d) | Sections 3.2, 3.3, 3.4.4 |
| `out/exp1/` | Experiment 1 records, Table 2 and the saved programme-layer founder object | Section 3.1; Table 2 |
| `out/exp2/` | Experiment 2 shards and tables; its founder populations were not saved, so a rerun repeats the experiment with new founders | Table 3 |
| `out/inference/`, `out/numbers.json` | inference tables; every number the paper quotes from them | Sections 3.2–3.5; Tables 3, 4, 6; Figures 1, 2 |
| `out/rerun_4001/` | replicate 4001 rerun from its stored founder object in a new session, and the comparison | Section 3.1 |
| `data/` | Tables 2 and 3 of Peixoto et al. (2024), one row per cell | Section 3.5 |
| `figures/` | Figures 1–5 as 600 dpi PNG and PDF | |
| `environment/` | R session (packages, BLAS, LAPACK, platform) and Python versions of the deposited run | Section 2.8 |
| `checksums/` | SHA-256 of the founder objects, raw records and regenerated outputs | |

`DATA_DICTIONARY.md` describes every column. `scripts/extract_published_tables.py` needs PyMuPDF and the
open-access PDF of Peixoto et al. (2024), which is not included.

## Seeds

| Use | Seeds |
|---|---|
| Application, one founder object and programme seed per replicate | 4001–4010 |
| Experiment 1 | founder seed 20260925; programme seed 7 |
| Experiment 2, four shards of ten replicates | base seeds 500000, 510000, 520000, 530000 (paired replicate r: base + r; unpaired arm k: base + 100000k + r) |
| Bootstraps: costs, pairing factors, inference | 20260925, 20261013, 20261014 |

## Licence

Code: MIT (`LICENSE`). Outputs (`out/`) and figures: CC BY 4.0. `data/peixoto2024_tables2_3.csv`: CC BY 4.0,
copyright the authors of Peixoto et al. (2024), G3 14: jkae128, doi:10.1093/g3journal/jkae128; values unchanged,
one row per cell.
