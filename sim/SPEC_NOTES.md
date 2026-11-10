# Implementation notes

Choices the engine makes where the method needs a concrete rule. Code references are to files in
`sim/R/`; Table 1 holds the parameter values.

## A. Founders (`founders.R`)

**A1. Founders simulated and kept.** Both numbers in Table 1 are per pool: `runMacs` simulates
2 × 180 individuals in the application and keeps 80 cull survivors in each pool (Experiment 1:
2 × 120 and 60; Experiment 2: 2 × 150 and 40).

**A2. Pool layout.** Of the individuals returned by `runMacs(split = 100)`, 1..n/2 form one pool and
n/2+1..n the other (mean absolute frequency difference 0.19 between the halves, 0.03 between odd and
even individuals).

**A3. Endosperm locus.** The segregating site with the largest absolute frequency difference between
the pools (the first in genome order on ties). The elite pool is the half in which allele 1 is more
frequent at that site, so allele 1 is the sweet allele of equation (1). The two pools are exchangeable
in the coalescent history.

**A4. Cull.** Each pool keeps its first survivors in MaCS order. If either pool has too few, the
coalescent draw is repeated (`nCoalescentDraws`, at most 10).

**A5. QTL and markers.** Each draw creates a fresh `SimParam` from the founder haplotypes;
`addTraitA(nQtlPerChr, mean = 0, var = 1)` draws the QTL and their effects, then
`addSnpChip(nSnpPerChr)` draws the markers from the remaining sites, so the two sets are disjoint. A
draw in which either set contains the endosperm locus is rejected; `nRedraw` counts the rejected
draws.

**A6. Trait scale.** `addTraitA(mean = 0, var = 1, gamma = FALSE)` scales the trait in the combined
founder population (all 2n simulated individuals, before the cull), with variance on divisor N.
AlphaSimR writes the genetic value as intercept + Σ α_ℓ (X_iℓ − 1), which is equation (2) with
a_0 = intercept − Σ α_ℓ.

**A7. Residual variance.** σ²_e = V_G × (1 − h²)/h² = 7/3 for h² = 0.30 (V_G = 1 by construction),
stored as the `SimParam` default, so every new population receives phenotypes.

**A8. Burn-in.** The kept elite founders are phenotyped; each burn-in cycle mates the n_p individuals
of highest phenotype at random into n_c crosses of n_o progeny (`randCross`; 40 × 15 in the
application, 20 × 10 in Experiment 2). The last burn-in generation is the cycle-1 elite nursery X_1.
The donor pool is not burned in.

**A9. t ha⁻¹ scale.** Mean and additive variance of the kept elite founders (after the cull, before
the burn-in), variance on divisor N. The trait is additive, so additive and genetic variance coincide.

**A10. Founder-level record columns.** `baseMean` = mean of X_1; `donorMean` = mean of the kept donor
founders; `realisedGap` = `baseMean` − `donorMean`; `maxFreqDiff` = absolute frequency difference at
the locus; `cullDonor` = share of simulated donor-pool individuals not homozygous non-sweet. The
founder object also stores `cullElite`, the sweet-allele frequencies of both pools, the QTL and
marker names, the elite founders' genetic values, `scale$sdA0`, `scale$varE`, `base$sdA` (additive
SD of X_1 in t ha⁻¹) and `base$h2Cycle1` = V_G(X_1)/(V_G(X_1) + σ²_e).

**A11. Sugar trait.** Not simulated, because only the marketability indicator is read. Its
parameters are kept in `cfg$sugar`, and `locusA` records the additive effect (5.82).

**A12. Seeds.** The founder object of replicate s is built after `set.seed(s)`. The program seed is
s + `seeds$programmeOffset` (offset 0) and is set with an explicit generator kind
(Mersenne-Twister, Inversion, Rejection) immediately before cycle 1, so a program run does not depend
on whether the founders were built or read from disk. Experiment 1 uses founder seed 20260925 and
program seed 7.

**A13. Saved founder object.** `founder_<seed>.rds` holds the seed, the software versions, the founder
configuration, the `SimParam`, X_1 with its phenotypes, the donor pool and the quantities of A10
(about 0.5 MB). Each file has a SHA-256 sidecar, and the directory a `SHA256SUMS` list. On reload, the
seed, founder configuration and checksum are checked. Every run works on a deep copy of the
`SimParam` (an R6 object), so runs share no state.

## B. The cycle (`cycle.R`, `prediction.R`)

**B1. Training archive.** Empty at the start of cycle 1; burn-in phenotypes are not used. Each cycle
adds 250 records and the archive keeps the latest 1,000, so cycle t trains on min(250t, 1,000)
records.

**B2. Training sample.** Drawn with `sample.int` without replacement from X_t, with the phenotypes
drawn when X_t was created.

**B3. Ridge regression.** `rrBLUP::mixed.solve(y, Z = M − 1)` with its defaults (REML, K = I,
intercept as the only fixed effect). A marker monomorphic in the archive receives an effect of zero.

**B4. Ties.** Rankings on GEBV keep nursery order for ties (stable sort); the same rule ranks
parents, bridge parents and admission candidates.

**B5. Mating pool.** The parents in GEBV-rank order followed by the admitted lines in admission
order, mated by `randCross` (random distinct pairs, no selfing). The order is the same in every arm,
so ELITE2M with an all-zero schedule reproduces CLOSED exactly.

**B6. Recombination.** The MaCS maize map (2 M per chromosome) with AlphaSimR's default crossover
model (gamma model, v = 2.6, p = 0, close to the Kosambi map).

**B7. Phenotypes.** Every population created (F1, F2, BRIDGE's second step, MAS seedlings, the new
nursery) receives phenotypes at creation, which consume the generator; donor-derived phenotypes are
never used.

## C. Arms (`arms.R`)

**C1. Donor crosses.** Six donors are drawn first, then six mates from the selected parents, each
with replacement; each donor × mate cross gives one F1 (`makeCross2`, donor as female), selfed into
150 F2.

**C2. BRIDGE.** One F1 per bridge parent, selfed into 45 lines. With fewer than 20 sweet first-step
F2, all sweet F2 are used and the second step has fewer families; with none, nothing is offered.

**C3. MAS.** 4 × 150 = 600 seedlings per F1 (132 under plot matching); every sweet homozygote is
grown and offered. The cap applies to admissions, not to the number offered.

**C4. Admission.** Equation (6) with "≥": a candidate exactly at b_t is eligible; largest GEBV first.

**C5. Elite controls.** An elite control admits the candidates ranked n_p + 1, ..., n_p + k_t on
GEBV: k_t = κ = 2 in every cycle for ELITE2, the inflow arm's `nAdmitted` for ELITE2M. For both,
`nCand` is the nursery size and `nClearBar`/`nClearTrue` count the whole nursery; `accDonor` and
`biasDonor` are NA. CLOSED records zero counts and NA for the three pool statistics.

**C6. ELITE2M.** `run_replicate()` runs each inflow arm first and passes its `nAdmitted` vector to
its ELITE2M arm in the same call, then checks the replay. ELITE2M arms run in the `baseline` and
`resource_matched` scenarios (`--elite2m-scenarios all` runs all four).

**C7. CLOSED and ELITE2.** Run in every scenario, so that "visual sorting changes nothing genetic"
and "plot matching leaves ELITE2 unchanged" are checked.

## D. Scenarios and plot accounting (`outcomes.R`, `config.R`)

**D1. Plot-matched counts.** From the 200-plot budget: DIRECT ⌊200/6⌋ = 33 F2 per cross; BRIDGE
⌊100/6⌋ = 16 per cross in step 1 and 20 families × 5 progeny in step 2; MAS 4 × 33 = 132 seedlings
per cross. The elite nursery has 40 crosses × 10 progeny.

**D2. Nursery size.** X_1 (600 plants) is common to every arm. A plot-matched inflow arm and its
ELITE2M grow 400 plants from the nursery produced in cycle 1 onward.

**D3. Plots.** `plotsGrown` = plants of the elite nursery X_t + plants of that cycle's introgression
nurseries; F1 plants are not counted. Under field discard every F2 and second-step line is a plot;
under visual sorting only sweet homozygotes are; MAS grows only its sweet seedlings.

**D4. Visual sorting.** A flag read only by the plot count; no draw and no decision depends on it.
`visual_matched` is plot matching with visual sorting.

## E. Records (`outcomes.R`)

**E1. Timing.** Population columns (`meanG` to `pctSweet`, `gapCycle`) are measured at step 3 on X_t;
admission columns (`nAdmitted` to `admitTrueGap`, `nF2screened`, `nSweetUsable`) describe steps 4–5 of
the same cycle; `plotsGrown` counts the plots of the same cycle; `year` = 3t.

**E2. `nProgeny`.** Progeny per cross of X_t: 15 at cycle 1 in every arm, and 10 from cycle 2 in the
plot-matched inflow arms and their ELITE2M arms.

**E3. Variances and accuracy.** `totalVarG` = variance (divisor N) of the genetic values of X_t, in
(t ha⁻¹)²; `genicVar` = equation (11), in (t ha⁻¹)²; `genicPct` = 100 × V_t/V(X_1); `accuracy` =
correlation of GEBV and genetic value over X_t; `accDonor` and `biasDonor` (equation 10, internal
scale) are computed over the offered donor-derived pool, NA with fewer than two members or no
variance; `admitTrueGap` = mean genetic value of the admitted lines − b*_t, internal scale.

**E4. `donorPct`.** A marker is diagnostic when X_1 is fixed for one allele and the donor pool
carries the other; `donorPct` = 100 × Σ (frequency of the diagnostic alleles in X_t) / Σ (their
frequency in the donor pool). It is 0 in CLOSED, ELITE2 and ELITE2M and 100 in the donor pool.

**E5. Run parameters.** `locusA` = 5.82 (A11), `h2` = 0.30, `nAdmit` = κ = 2, `nDonorCross` = 6;
`rep` = seed − 4000; `RVersion`, `AlphaSimR` and `ompThreads` record the software and thread count.

**E6. Output files.** `run_application.R` writes `raw_<seed>.csv` with every arm and scenario
(columns in `RECORD_COLUMNS` order); with `--pipeline-layout` it also writes `raw_w<seed>.csv`
(CLOSED, ELITE2, DIRECT, BRIDGE, MAS) and `raw_m<seed>.csv` (ELITE2M), kept in `out/application/` and
`out/matched_controls/`.

## F. Experiments (`run_exp1.R`, `run_exp2.R`)

**F1. Experiment 1, generator layers.** `runMacs(nInd = 40, nChr = 3, segSites = 150,
species = "MAIZE", split = 100, nThreads = 1)` and `quickHaplo(nInd = 40, nChr = 3, segSites = 150)`,
each after `set.seed(20260925)`: five calls in one session and one call in each of three sessions.
Signature: mean dosage | dosage sum of individual 1 | number of sites.

**F2. Experiment 1, program layer.** One founder object built under seed 20260925 with the
Experiment 1 column of Table 1, saved and read back; CLOSED, DIRECT, BRIDGE and MAS run five times on
the in-memory object and once on the reloaded one, program seed 7, four cycles. Signature: final gain
| final genic variance, ten decimals.

**F3. Experiment 2, seeds.** Paired replicate r uses seed base + r for its founders and every arm's
program. Unpaired replicate r, arm k uses base + 100000k + r, with k = 1, 2, 3, 4 for CLOSED, DIRECT,
BRIDGE, MAS, for both its founders and its program.

**F4. Experiment 2, outcomes.** `gain` and `genicPct` at cycle 10 and `accLate`, the mean accuracy
over cycles 6–10, in the baseline scenario. Founders are not saved, so a rerun reproduces the outputs
only in distribution.

## G. Software

`OMP_NUM_THREADS = 1` is set before AlphaSimR is loaded, and every `SimParam` and `runMacs` call
uses `nThreads = 1`. With R's reference BLAS the results are as deposited; an optimised or threaded
BLAS can change the last digits of `mixed.solve`.
