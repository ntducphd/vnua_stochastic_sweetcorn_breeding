# sim/: breeding-programme simulation engine (clean-room rewrite)

This folder holds a new, independent R implementation of the stochastic sweet corn breeding
programme described in Section 2 of the manuscript (`internal/R02_stochastic_sweetcorn_workflow.md`):
the founder object, the cycle operator of Algorithm 1, the five programmes and the count-matched
controls, the four scenarios, and the drivers for the application and Experiments 1 and 2. It was
written from the paper, `DATA_DICTIONARY.md` and the AlphaSimR/rrBLUP help pages only. It writes the
per-cycle records and the Experiment 1/2 tables in the layouts the existing analysis scripts read.

Every implementation choice where the paper is silent or ambiguous is listed, with the sentence it
interprets, in **[SPEC_NOTES.md](SPEC_NOTES.md)**. Read that file before aligning the Methods text.

## Requirements

R 4.5.1 with AlphaSimR 2.1.0 and rrBLUP 4.6.3 (the versions in Section 2.8). No other packages are
needed; checksums use `tools::sha256sum()` (R ≥ 4.5). Every script sets `OMP_NUM_THREADS=1` before
loading AlphaSimR and passes `nThreads = 1` to AlphaSimR.

## How to run

All commands are run from the repository root.

```sh
# tests (about 3 minutes; writes sim/tests/test_results.txt)
Rscript sim/tests/run_tests.R

# one replicate of the application: all four scenarios, the five arms and the ELITE2M controls
Rscript sim/run_application.R --seed 4001 --out out_sim/application --founders out_sim/founders

# all ten replicates (each call skips a replicate whose CSV already exists, so the loop can be
# restarted after an interruption; replicates can also run as separate processes)
for s in $(seq 4001 4010); do
  Rscript sim/run_application.R --seed $s --out out_sim/application --founders out_sim/founders \
      --pipeline-layout
done

# Experiment 1 (layers 1-4 in one session, then layers 1-2 in three separate sessions)
Rscript sim/run_exp1.R all --out out_sim/exp1

# Experiment 2: four shards of ten replicates (base seeds 500000, 510000, 520000, 530000)
for k in 0 1 2 3; do Rscript sim/run_exp2.R --shard $k --out out_sim/exp2; done
```

Useful options of `run_application.R`:

| Option | Effect |
|---|---|
| `--config application\|exp1\|exp2\|test` | parameter set (default `application`, Table 1) |
| `--scenarios a,b` | subset of `baseline,resource_matched,visual_sort,visual_matched` |
| `--arms a,b` | subset of `CLOSED,ELITE2,DIRECT,BRIDGE,MAS` |
| `--elite2m-scenarios a,b\|all\|none` | where ELITE2M arms run (default baseline and resource_matched) |
| `--pipeline-layout` | also write `raw_w<seed>.csv` and `raw_m<seed>.csv` (see below) |
| `--cycles N` | shorten the run (timing only) |
| `--force` | rerun a replicate whose CSV exists |

A saved founder file **defines** its replicate: if `--founders <dir>/founder_<seed>.rds` exists it
is reloaded (its seed, founder configuration and SHA-256 checksum are checked), otherwise it is built
under `set.seed(seed)` and saved. Under AlphaSimR 2.1.0 `runMacs` does not reproduce from a seed, so
keep the founder files with the outputs.

### Outputs read by the analyses

With `--pipeline-layout` the application driver writes, per seed, `raw_w<seed>.csv` (the five
programmes in the four scenarios) and `raw_m<seed>.csv` (the count-matched controls). The deposit keeps
them in `out/application/` and `out/matched_controls/`, where `R/01_pool_application.py` and
`R/04_analyse_matched_controls.R` read them. Experiment 2 shards (`shard<k>.csv`, 80 rows each) are read
by `R/exp2_pool.py`, and `R/03_founder_scale.R` reads the founder objects (fields listed in SPEC_NOTES A10).

## Files

| File | Contents |
|---|---|
| `R/load.R` | `load_sim(simDir)`: one thread, packages, then the modules in order |
| `R/config.R` | parameter sets `config_application()`, `config_exp1()`, `config_exp2()`, `config_test()`; `validate_config()`, `modify_config()`, resource-matching derivation, the scenario table |
| `R/utils.R` | seeding with an explicit generator kind, generator-state capture, logging, ranking, variance and correlation helpers, CLI parser |
| `R/founders.R` | the founder object: coalescent draw, endosperm locus, cull, QTL/marker redraws, residual variance, t/ha scale, burn-in; save/reload with SHA-256 checksums |
| `R/genetics.R` | locus dosage and marketability (equation 3), scale conversion, genic variance (equation 11), donor-allele content |
| `R/prediction.R` | training archive (latest records kept) and ridge regression with `rrBLUP::mixed.solve` (equation 4) |
| `R/arms.R` | the arms as strategy objects (`offer`, `admit`, nursery size): CLOSED, ELITE2, DIRECT, BRIDGE, MAS, `ELITE2M_<arm>`; admission rule (equation 6) |
| `R/cycle.R` | `run_cycle()`: Algorithm 1, steps 1–6, with an optional trace of generator states |
| `R/programme.R` | `run_programme()` (one arm, one scenario) and `run_replicate()` (all arms and scenarios plus ELITE2M replay) |
| `R/outcomes.R` | record schema (`RECORD_COLUMNS`), step-3 measurements, supply and bias (equations 8–10), plot accounting |
| `R/io.R` | atomic CSV writing, the pipeline layout, append helpers, timing log |
| `run_application.R` | driver: one replicate per call, resumable, timing log |
| `run_exp1.R` | driver: Experiment 1 (`within`, `across`, `sessions`, `all`) |
| `run_exp2.R` | driver: one Experiment 2 shard, paired and unpaired designs, resumable per block |
| `tests/run_tests.R` | test suite (13 tests); `tests/rerun_in_new_session.R` is its helper |
| `tests/test_results.txt` | output of the last test run |
| `SPEC_NOTES.md` | implementation choices, with the paper text each interprets |

## Outputs

| File | Rows | Columns |
|---|---|---|
| `raw_<seed>.csv` (application) | one per scenario × arm × cycle: (4 × 5 + 2 × 3) × 20 = 520 | the 42 per-cycle columns of DATA_DICTIONARY.md, in its order |
| `raw_w<seed>.csv`, `raw_m<seed>.csv` | main arms (400) and ELITE2M arms (120) | same |
| `timing_<seed>.tsv` | founder build or load, each run, total | time, seed, scenario, arm, seconds |
| `founders/founder_<seed>.rds`, `.rds.sha256`, `SHA256SUMS` | one per replicate | founder object (SPEC_NOTES A13) |
| `exp1/within_session.csv` | 2 × 5 genome calls + 4 × 6 programme calls | layer, rep, signature, note |
| `exp1/across_sessions.csv` | 2 per session | session, pid, layer, signature, AlphaSimR, R |
| `exp1/table2.csv` | one per layer | layer, calls, distinct, reproducible |
| `exp2/shard<k>.csv` | 10 replicates × 2 designs × 4 arms = 80 | design, seedBase, rep, arm, gain, genicPct, accLate, seed |

## Design in brief

* **State.** A programme state is the elite nursery (an AlphaSimR `Pop` with phenotypes drawn at
  creation), the progeny per cross it was grown with, and the training archive.
* **Cycle.** `run_cycle()` is one function for every arm. Its random draws come only from R's
  generator, in the order of Algorithm 1: training sample (step 1); donor and mate draws, crosses
  and selfs (step 4, inflow arms only); random mating (step 6). Steps 2, 3 and 5 draw nothing.
* **Arms.** An arm supplies three things: what it offers (`offer`), what it admits (`admit`) and the
  size of the nursery it grows. CLOSED offers nothing; DIRECT, BRIDGE and MAS offer donor-derived
  sweet homozygotes and admit by the threshold and cap; ELITE2 and ELITE2M offer the elite nursery
  and admit the next-ranked candidates on a schedule (all twos for ELITE2, the inflow arm's
  admissions for ELITE2M). Because ELITE2 is the elite control with a constant schedule, and the
  mating pool is built the same way in every arm, ELITE2M reproduces CLOSED (all zeros) and ELITE2
  (all twos) exactly.
* **Pairing.** Every arm of a replicate starts from the same founder object and programme seed, so
  arms share the first training sample and parents, and part at step 4 of cycle 1 (inflow arms) or
  at step 6 (elite controls). The test suite checks this on the generator state itself.
* **Scale.** Genetic values stay on AlphaSimR's internal scale and are converted to t/ha only when a
  record is written; `biasDonor` and `admitTrueGap` stay on the internal scale, as the dictionary says.

## Compute

Single-threaded on this laptop (R 4.5.1, reference BLAS):

| Unit | Time |
|---|---|
| founder object, application size (360 individuals × 10 chromosomes in MaCS, 8 burn-in cycles) | about 4 min (236 s) |
| one 20-cycle run of one arm, application size | 75–92 s, of which about 85% is `mixed.solve` at 1,000 records (≈ 4.5 s per cycle from cycle 4) |
| baseline scenario of one replicate (5 arms + 3 ELITE2M, 20 cycles), measured for seed 4001 | 698 s for the runs + 236 s founder build = **15.6 min** |
| full replicate (4 scenarios × 5 arms + ELITE2M in 2 scenarios, founder build included) | **about 42 min** (estimated: 236 + 4 × 453 + 2 × 245 s) |
| ten replicates of the application | about 7 h sequentially; replicates are independent processes |

Measured on 16 October 2026 with `run_application.R --seed 4001 --scenarios baseline` (per-arm
seconds: CLOSED 79, ELITE2 91, DIRECT 92, BRIDGE 98, MAS 94, ELITE2M_DIRECT 75, ELITE2M_BRIDGE 89,
ELITE2M_MAS 81). In that run the elite pool was homozygous sweet in every cycle, DIRECT and MAS
recovered 0.2531 and 0.2517 of their F2 as sweet homozygotes and BRIDGE 0.6245, DIRECT offered 228
and MAS 906 candidates per cycle, and the founder gap was 5.75 t/ha (base 21.17, donor 15.42).
Experiment 1 at full size would take about 8 min (founders about 2.6 min) and one Experiment 2
shard about 75 min, most of it in `runMacs` (50 founder populations per shard).
