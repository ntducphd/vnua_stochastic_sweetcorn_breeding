# Data dictionary

All files are comma-separated with a header row. Yield is on the t ha⁻¹ scale set by the elite
founders (mean 14 t ha⁻¹, additive variance 1 (t ha⁻¹)², both assumed) unless a column says it is
on the simulator's internal genetic scale; `out/application/founder_scale.csv` gives each replicate's
factor between the two.

## Per-cycle records: `out/application/raw_w<seed>.csv`, `out/application/pooled_records.csv`, `out/matched_controls/raw_m<seed>.csv`, `out/matched_controls/pooled_records.csv`, `out/rerun_4001/raw_[wm]4001.csv`

One row per replicate, scenario, arm and cycle, written by `run_programme()` in `sim/R/programme.R`
(columns in the order of `RECORD_COLUMNS` in `sim/R/outcomes.R`). `sim/run_application.R
--pipeline-layout` runs one replicate (seed) per call, with all four scenarios in one run, and
splits its records into `out/application/raw_w<seed>.csv` (CLOSED, ELITE2, DIRECT, BRIDGE and MAS:
4 scenarios × 5 arms × 20 cycles = 400 rows) and `out/matched_controls/raw_m<seed>.csv` (the
ELITE2M arms, run in the baseline and resource-matched scenarios only: 2 × 3 × 20 = 120 rows).
`R/01_pool_application.py` stacks the ten `raw_w<seed>.csv` files in seed order into
`out/application/pooled_records.csv` (4,000 rows), and `R/04_analyse_matched_controls.R` stacks the
`raw_m<seed>.csv` files into `out/matched_controls/pooled_records.csv` (1,200 rows).
`out/rerun_4001/` holds the rerun of replicate 4001 from its saved founder object in a new R session
(baseline and resource-matched scenarios) and `comparison.txt`, written by `scripts/compare_rerun.py`.

The population columns (`meanG` to `pctSweet`, `gapCycle`, `nProgeny`) are measured at step 3 of
Algorithm 1 on the unselected nursery of the cycle, **before** that cycle's selection and admission,
so the admissions of cycle 20 do not affect cycle 20's outcomes. The admission columns (`nAdmitted`
to `admitTrueGap`, `nF2screened`, `nSweetUsable`) describe what steps 4 and 5 of the same cycle
offered and admitted, and `plotsGrown` counts the plots of the same cycle.

| Column | Meaning | Unit or scale |
|---|---|---|
| `cycle`, `year` | cycle 1 to 20; year at its end (three years per cycle) | |
| `meanG` | mean genetic value of the elite nursery | t ha⁻¹ |
| `gain` | `meanG` minus the post-burn-in base (`baseMean`); zero at cycle 1 | t ha⁻¹ |
| `genicVar` | genic variance, equation (11) | (t ha⁻¹)² |
| `genicPct` | genic variance as a percentage of the base | % |
| `totalVarG` | realised genetic variance of the nursery | (t ha⁻¹)² |
| `accuracy` | correlation of GEBV and true genetic value among all the elite candidates of the nursery | |
| `donorPct` | donor-allele content of the nursery at the diagnostic markers (definition below) | % |
| `pctSweet` | share of the nursery homozygous sweet at the endosperm locus | % |
| `nAdmitted` | lines admitted this cycle, at most 2 | count |
| `nCand` | candidates offered for admission (donor-derived lines; for ELITE2 and ELITE2M the whole elite nursery; 0 for CLOSED) | count |
| `nClearBar` | candidates clearing the threshold on GEBV, before the cap, equation (8) | count |
| `nClearTrue` | candidates clearing the true-value threshold, equation (9) | count |
| `accDonor` | correlation of GEBV and true value inside the donor-derived pool; NA for CLOSED, ELITE2 and ELITE2M, or with fewer than two candidates | |
| `biasDonor` | mean GEBV minus true value inside the donor-derived pool, equation (10); NA as for `accDonor` | **internal genetic scale** |
| `admitTrueGap` | admitted lines' mean true value minus the true threshold; NA when nothing is admitted | **internal genetic scale** |
| `gapCycle` | elite nursery mean minus the frozen donor pool mean | t ha⁻¹ |
| `nRedraw` | rejected draws of the QTL and marker sets before neither contained the endosperm locus (0 when the first draw was valid); founder-level | count |
| `nF2screened`, `nSweetUsable` | F2 (and BRIDGE second-step) progeny or MAS seedlings screened, and sweet homozygotes among them | count |
| `nProgeny` | progeny per cross of the elite nursery measured in the row: 15, except from cycle 2 on in the inflow arms and their ELITE2M controls under resource matching, where it is 10 (the cycle-1 nursery is common to every arm) | count |
| `plotsGrown` | plots grown this cycle: plants of the elite nursery plus plants of that cycle's introgression nurseries (every F2 and second-step line under field discard, only sweet homozygotes under visual sorting, only sweet seedlings in MAS); F1 plants are not counted | count |
| `strategy`, `scenario`, `visualSort`, `resourceMatch` | arm and scenario | |
| `rep`, `seed` | replicate (seed − 4000) and its seed, which is both the founder seed and the programme seed of the replicate | |
| `RVersion`, `AlphaSimR`, `ompThreads` | software and thread count of the run (4.5.1, 2.1.0, 1) | |
| `locusA`, `h2`, `nAdmit`, `nDonorCross` | parameters of the run: the additive effect of the sugar trait at the endosperm locus (5.82; recorded, not simulated), the heritability in the combined founder population (0.30), the admission cap κ (2) and the donor crosses per cycle (6) | |
| `realisedGap`, `donorMean`, `baseMean`, `locusChr`, `maxFreqDiff`, `cullDonor` | founder-level quantities of the replicate: `baseMean`, mean of the cycle-1 elite nursery (after burn-in); `donorMean`, mean of the kept donor founders; `realisedGap` = `baseMean` − `donorMean`; `locusChr`, chromosome of the endosperm locus; `maxFreqDiff`, absolute difference in sweet-allele frequency between the pools at the locus; `cullDonor`, share of the simulated donor-pool individuals culled for not being homozygous non-sweet | t ha⁻¹ where applicable |

`donorPct` is computed by `donor_allele_pct()` in `sim/R/genetics.R`. A marker is *diagnostic* when
the post-burn-in elite nursery of cycle 1 is fixed for one allele and the frozen donor pool carries
the other. With p_j the frequency in the population of the allele that the cycle-1 nursery lacks at
diagnostic marker j, and q_j its frequency in the donor pool,
`donorPct` = 100 × Σ_j p_j / Σ_j q_j. It is exactly 0 for any population descending from the cycle-1
nursery alone (CLOSED, ELITE2, ELITE2M: there is no mutation), 100 for the donor pool itself and
about 50 in expectation for unselected donor × elite F2. It tracks marker alleles, not identity by
descent, and is not used in the paper.

## Bootstrap intervals: `out/matched_controls/bootstrap_intervals.csv` (`R/05_bootstrap.py`)

Paired percentile intervals, 4,000 resamples of replicates, seed 20261013: the paired-only gain
factors of the application (Fig. 1d, Section 3.2) and of Experiment 2, and the Table 4 median
relative standard deviations with the slope-to-endpoint ratio of medians (rows `table4 ...`).
Columns: `quantity`, `estimate`, `lo`, `hi`, `resamples`. Table 3 uses the Experiment 2 intervals
in `out/inference/exp2_robust.csv` (seed 20261014).

## Per-run summaries

`out/application/per_run_summary.csv` (`R/02_analyse_application.R`), one row per replicate, scenario and arm:
`estEndpoint` (cumulative gain at cycle 20, t ha⁻¹), `estLateMean` (mean gain over cycles 16 to
20), `estSlope` (OLS slope of gain on year over cycles 16 to 20, t ha⁻¹ yr⁻¹), `genicEnd`,
`genicLate`, `accEnd`, `accLate`, `donorEnd`, the per-cycle means `candPerCyc`, `clearBarPerCyc` and
`clearTruePerCyc` (of `nCand`, `nClearBar`, `nClearTrue`), the run means of `accDonor`, `biasDonor` and
`admitTrueGap` (internal genetic scale), total `admitted`, the per-cycle means `plotsPerCyc` (of `plotsGrown`),
`f2screened` and `sweetUsable`, the elite–donor gap at the start and end (`gapStart`, `gapEnd`), and `visualSort`.

`out/matched_controls/per_run_summary.csv` (`R/04_analyse_matched_controls.R`): `gainEnd`, `genicEnd` and
`admitted` per replicate, scenario and arm, including the ELITE2M arms.

## Analysis tables

- `out/application/founder_summary.csv` (`R/02_analyse_application.R`): one row per replicate, `rep`, `seed`
  and the founder-level columns of the per-cycle records.
- `out/application/founder_scale.csv` (`R/03_founder_scale.R`): per founder object, `scale` (t ha⁻¹ per
  internal unit), `sdA_tha` (additive sd after burn-in, t ha⁻¹), `h2_cycle1`, `varE` (internal scale),
  `locusChr`, `nRedraw`, `baseMean`, `donorMean`, `gap` (t ha⁻¹), `AlphaSimR`, `R`.
- `out/application/paired_contrasts.csv` (`R/02_analyse_application.R`): per scenario and arm against CLOSED,
  for each statistic `<f>` of the per-run summary, `d_<f>` (mean paired difference), `lo_<f>`, `hi_<f>` (95% t
  interval), `bonLo_<f>`, `bonHi_<f>` (Bonferroni over 16 contrasts) and `sd_<f>`.
- `out/application/power.csv` (`R/02_analyse_application.R`): endpoint contrast `diff`, `sd`, `mde80`
  (t ha⁻¹) and `significant` (nominal 5%).
- `out/application/estimator_comparison.csv` (`R/02_analyse_application.R`): per scenario, arm and
  `estimator`, `meanDiff`, `sdDiff` and `cv` = sdDiff / |meanDiff|.
- `out/application/locus_cost.csv` (`R/02_analyse_application.R`): per-run means of `f2screened`,
  `sweetUsable`, `plotsPerCyc` and `admitted`; `plotsPerLine` (total plots / total admitted lines) with its
  bootstrap interval `plotsPerLineLo`, `plotsPerLineHi` (4,000 resamples, seed 20260925); `zeroAdmitReps`
  of `nReps`.
- `out/application/supply.csv` (`R/02_analyse_application.R`): means of `cand`, `clearBar` and `clearTrue`
  per cycle and of `accDonor`, `biasDonor` and `admitTrueGap` (internal genetic scale).
- `out/application/prediction_bias_tha.csv` (`R/03_founder_scale.R`): Table 5 in t ha⁻¹, `bias` (from
  `biasDonor`) and `gap` (from `admitTrueGap`) with `biasN`, `biasLo`, `biasHi`, `gapN`, `gapLo`, `gapHi`
  (95% t interval over replicates), and `clearBar`, `clearTrue` per cycle.
- `out/matched_controls/matched_contrasts.csv` (`R/04_analyse_matched_controls.R`): per scenario, arm and
  `contrast`, gain at cycle 20: `n`, `diff`, `sd`, `lo`, `hi` (95% t), `bonLo`, `bonHi` (Bonferroni over
  six), and the mean lines admitted by the two arms compared, `admittedA`, `admittedB`.
- `out/matched_controls/estimator_distinct.csv` (`R/04_analyse_matched_controls.R`): Table 4 over the seven
  distinct contrasts: `nContrasts`, `medianRelSd`, `minRelSd`, `maxRelSd`, `leastVariable`, `mostVariable`,
  `ratioToEndpoint`, `bonferroniT`.

## Founder objects: `out/founders/`

`founder_<seed>.rds` (seeds 4001 to 4010, about 0.5 MB each), written by `sim/run_application.R`
(`save_founders()` in `sim/R/founders.R`), is the founder object of one replicate: its seed, the
software versions, the founder part of the configuration, the AlphaSimR `SimParam`, the cycle-1 elite
nursery with its phenotypes, the kept donor founders, the endosperm locus, the QTL and marker names,
`nRedraw`, the cull shares, the t ha⁻¹ scale and the diagnostic markers behind `donorPct`. Each file has
a SHA-256 sidecar `founder_<seed>.rds.sha256`, and `SHA256SUMS` lists all ten. A saved object is
reloaded only after its seed, configuration and checksum have been checked. Under AlphaSimR 2.1.0 the
seed does not regenerate the founders (Experiment 1), so these files, not the seeds, define the
replicates.

## Experiments 1 and 2

- `out/exp1/within_session.csv` (`sim/run_exp1.R`): `layer`, `rep` (call), `signature` (mean dosage |
  dosage sum of the first individual | number of sites for a genome; final gain | final genic variance
  for a programme run, ten decimals), `note` (for the programme layers, whether the founder object was
  in memory or recovered with `saveRDS`/`readRDS`).
- `out/exp1/across_sessions.csv`: `session`, `pid` (process id of the R session), `layer`, `signature`,
  `AlphaSimR`, `R`.
- `out/exp1/table2.csv`: Table 2, per layer `calls`, `distinct` (distinct signatures) and
  `reproducible`.
- `out/exp1/founder_20260925.rds` and its `.sha256`: the founder object of the Experiment 1 programme
  layers (founder seed 20260925; Table 1, Experiment 1 column). The `runMacs` and `quickHaplo` genomes of
  layers 1 and 2 are not saved.
- `out/exp2/shard[0-3].csv` (`sim/run_exp2.R`): `design` (paired or unpaired), `seedBase`, `rep`,
  `arm`, `gain` (cumulative gain at cycle 10, t ha⁻¹), `genicPct` (genic variance at cycle 10, % of
  its post-burn-in value), `accLate` (mean accuracy over cycles 6 to 10) and `seed` (founder and
  programme seed of the run: seedBase + rep for a paired replicate; seedBase + 100000 k + rep for arm
  k = 1, 2, 3, 4 = CLOSED, DIRECT, BRIDGE, MAS of an unpaired replicate). The founder populations of
  Experiment 2 are not saved.
- `out/exp2/pooled.csv` (`R/exp2_pool.py`): the four shards joined, with `uid` (= seedBase_rep) in
  place of `seed`.
- `out/exp2/variance_reduction_ci.csv` (`R/exp2_summary.py`): per outcome and arm, `n`, the sd of
  the unpaired and paired contrasts (`sdUnpaired`, `sdPaired`), the paired correlation (`corPaired`),
  the variance-reduction factor (`varianceReduction`) with its F interval (`ciLo`, `ciHi`),
  `nUnpairedForSamePrecision` and a `verdict`.

## Statistical inference: `out/inference/` (`R/06_inference.py`)

- `exp2_robust.csv`: per outcome and arm, the variance-reduction factor with F, bootstrap and
  Bonferroni-over-nine limits, Brown–Forsythe p, Shapiro–Wilk p and excess kurtosis of both
  contrasts, leave-one-out range, correlation with its Fisher interval, the paired-only estimate
  of equation (17) with its bootstrap interval, and the exact form of equation (17) from the
  correlation and the ratio of standard deviations.
- `rw_contrasts.csv`: per distinct contrast, the drift-scale variance ratios of the late slope
  and the late mean to the endpoint with bootstrap intervals, Pitman–Morgan p against the
  constant-increment prediction and against 1, mean late slope per cycle, mean endpoint / 19,
  lag-one correlation of the increments and Var(D₂₀)/Var(D₁₀).
- `rw_by_cycle.csv`, `rw_by_cycle_contrast.csv`: variance of the contrast and of its increments
  by cycle, with the random-walk and random-drift fits.
- `contrasts_inference.csv`: every contrast of the application by family (`closed7`, `donor6`,
  `recycling6`, and the descriptive `unmatched6` and `elite2m_vs_elite2_6`): mean, sd, t interval,
  t, raw, Bonferroni and sign-flip max-|t| p, both critical values, nominal and corrected
  equivalence bounds and minimum detectable effects.
- `published_reanalysis.csv`: Section 3.5. Per combination of programme size (`crosses`),
  `heritability`, `parameter` (`genetic_mean`, `genetic_variance`), pipeline (`program`) and
  genotype-by-environment interaction (`ge`) in Peixoto et al. (2024), with the ratios read in the
  order of their Table 3 (`ratio` in the input): the two ratios compared (`first`, `second`: on
  the genetic mean the best and the second best, on genetic variance the lowest ratio other than
  0:1 and 0:1) and their means (`mean_first`, `mean_second`); the difference `diff`, as a
  percentage of `mean_first` (`diff_pct`) and in founder additive standard deviations, sqrt(30)
  (`diff_founder_sd`); `z` with the spreads read as standard deviations over 50 independent
  replicates, its range within the rounding of the published values (`z_rounding_lo`,
  `z_rounding_hi`), the Welch p-value (`welch_p`) and z with a correlation of 0.7 between the two
  scenarios (`z_rho07`); the smallest difference 50 independent replicates detect with 80% power
  (`mde50`); the replicates per scenario that detect a difference of 1% of `mean_first` with 80%
  power (`n_1pct`); and z if the spreads were standard errors (`z_if_se`).

## Published input: `data/peixoto2024_tables2_3.csv`

Tables 2 (50 crosses) and 3 (200 crosses) of Peixoto et al. (2024), G3 14: jkae128,
doi:10.1093/g3journal/jkae128, CC BY 4.0, read by `scripts/extract_published_tables.py` from the
coordinates of the words on page 6 of the published PDF and checked against the page image. One
row per cell: `table`, `crosses`, `heritability`, `parameter`, `program` (Conv, GSTC, GSF1), `ge`
(absent, present), `ratio_printed` (the column label as printed: Table 2 prints 1:3, 1:1, 3:1,
0:1 and Table 3 prints 3:1, 1:1, 1:3, 0:1), `ratio` (both tables in the order of Table 3; Section 2.7), `mean`
over 50 replicates and `spread`, the number in parentheses, read as a standard deviation
(Section 2.7). Values are unchanged from the article; the layout and the
`ratio` column are added here. Source: Peixoto MA, Coelho IF, Leach KA, Lübberstedt T, Bhering LL,
Resende MFR (2024), copyright the authors, licensed under CC BY 4.0
(https://creativecommons.org/licenses/by/4.0/).

## `out/numbers.json`

A dictionary `key → {"value": number, "text": string}`. `text` is the number as the paper prints it.
