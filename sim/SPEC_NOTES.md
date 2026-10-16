# Implementation notes: where the paper is silent or ambiguous

This engine (`sim/`) was written as a clean-room implementation from the manuscript
(`internal/R10_stochastic_sweetcorn_workflow.md`, Sections 2.1–2.4, Table 1 and the Results),
`DATA_DICTIONARY.md`, the column names read by the downstream scripts, and the AlphaSimR 2.1.0 and
rrBLUP 4.6.3 help pages. Each entry below quotes the text it interprets, states what the engine does,
and, where it matters, what the authors may want to change in the Methods. Items marked **[check]**
are places where the engine's behaviour and a number in the paper may not agree.

Code references are to files in `sim/R/`.

---

## A. Founders (`founders.R`)

**A1. Founders simulated versus kept per pool.**
Paper: Table 1, "Founders simulated and kept per pool | 180 and 80 | 120 and 60 | 150 and 40";
Section 2.3, "run for four cycles after a three-cycle burn-in from 60 elite and 60 donor founders".
Engine: both numbers are per pool. runMacs simulates 2 × 180 = 360 individuals (application), and
80 survivors of the cull are kept in each pool. Experiment 1: 2 × 120 simulated, 60 kept per pool;
Experiment 2: 2 × 150 simulated, 40 kept per pool.

**A2. Which simulated individuals form which pool.**
Paper: Table 1, "MaCS maize history in AlphaSimR, pools split 100 generations ago".
The runMacs help page does not say how the two populations of `split` are laid out. Checked
empirically: individuals 1..n/2 form one pool and n/2+1..n the other (mean |frequency difference|
0.19 between halves against 0.03 between odd and even individuals). The engine uses that layout.

**A3. The endosperm locus and the direction of the sweet allele.**
Paper: "It is the segregating site with the largest difference in allele frequency between the two
founder pools"; equation (1), "X_iℓ counts copies of the allele coded 1 ... (at the endosperm
locus, the sweet allele)".
Engine: the largest *absolute* difference over all segregating sites (first site in genome order on
ties). The **elite pool is the half in which allele 1 is the more frequent** at that site, so that
allele 1 is the sweet allele as equation (1) requires. The two pools are exchangeable in the MaCS
history, so assigning the elite role this way does not change the distribution of anything.

**A4. Cull failure.**
Paper: "The founders kept in each pool are the first that survive the cull" and "no run failed";
silent on what happens if a pool has fewer survivors than must be kept.
Engine: keeps the first survivors in MaCS order. If either pool has too few, the coalescent draw is
repeated (counter `nCoalescentDraws` in the founder object; at most 10). At the application sizes
this did not occur in the replicate built for timing (elite cull 27%, donor cull 0.6%).

**A5. Drawing QTL and markers, and counting redraws.**
Paper: "the QTL and marker sets are redrawn until neither contains it"; Table 1, "QTL and markers per
chromosome (disjoint) | 100 and 100".
Engine: for each draw, a fresh `SimParam` is created from the founder haplotypes;
`addTraitA(nQtlPerChr, mean = 0, var = 1)` draws the QTL and their effects, then
`addSnpChip(nSnpPerChr)` draws the markers from the sites not taken as QTL (AlphaSimR's default
eligibility rule), so the sets are disjoint. If either set contains the locus the whole draw is
rejected. `nRedraw` is the number of rejected draws (0 when the first draw is valid).
AlphaSimR's `restrSegSites(excludeQtl, excludeSnp)` would avoid redraws altogether but would not
be the procedure the paper describes, so it is not used.
**[check]** Under this procedure a draw is valid with probability (250 − 200)/250 = 0.2 at the
application sizes, so `nRedraw` is geometric with mean 4 and median 3. Section 3.4.1 reports
"0 to 33 redraws per replicate (median 2.5)"; 33 or more has probability 0.8^33 ≈ 6 × 10⁻⁴ per
replicate under this scheme. If the original code drew the sets differently (for example
chromosome by chromosome), the Methods should say how.

**A6. Trait scale, intercept and effect distribution.**
Paper: equation (2), "α_ℓ ~ N(0, σ²_α) ... a_0 and σ²_α are set so that the genetic values of the
founder population have mean 0 and variance 1".
Engine: AlphaSimR `addTraitA(mean = 0, var = 1, gamma = FALSE)`, scaled in the *combined* founder
population (all 2 × n simulated individuals, both pools, before the cull), with variance on divisor
N (AlphaSimR's convention). AlphaSimR writes the genetic value as intercept + Σ α_ℓ (X_iℓ − 1);
this is equation (2) with a_0 = intercept − Σ α_ℓ.

**A7. Residual variance.**
Paper: "The residual variance is set from a heritability of 0.30 in the combined founder population
of both pools, before the cull at the endosperm locus, and then held fixed."
Engine: σ²_e = V_G(founders) × (1 − h²)/h² = 7/3 (V_G = 1 by construction), stored as the
`SimParam` default, so AlphaSimR draws a phenotype for every individual of every new population.

**A8. Burn-in scheme.**
Paper: Table 1, "Burn-in by phenotypic selection (cycles) | 8 | 3 | 5"; no further detail.
Engine: the kept elite founders are phenotyped; in each burn-in cycle the n_p individuals of highest
phenotype are mated at random into n_c crosses of n_o progeny (`randCross`, Table 1 programme sizes:
40 × 15 in the application, 20 × 10 in Experiment 2); the phenotypes of each new generation are
drawn as it is created. With 80 kept founders, the first cycle selects 40 of 80. The last burn-in
generation is the cycle-1 elite nursery X_1, with the phenotypes from which the first training
sample is drawn. The donor pool is not burned in ("unselected donor founders", "the donor pool is
frozen").

**A9. Which "elite founders" define the t/ha scale.**
Paper: "ḡ_0 and σ²_A0 the mean and additive variance of the elite founders".
Engine: the **kept** elite founders (after the cull, before the burn-in), variance on divisor N.
The trait is purely additive, so additive and genetic variance coincide. Using all simulated
elite-pool individuals instead would shift the scale slightly.

**A10. Founder-level quantities written to every record.**
`baseMean` = mean of X_1 (t/ha); `donorMean` = mean of the kept donor founders; `realisedGap` =
`baseMean` − `donorMean`; `maxFreqDiff` = the absolute frequency difference at the locus;
`cullDonor` = share of the simulated donor-pool individuals that were not homozygous non-sweet.
The founder object also stores `cullElite`, the frequencies of the sweet allele in both pools, the
QTL and marker names, the elite founders' genetic values, and, for a founder-scale adapter,
`scale$sdA0` (internal units per t/ha), `scale$varE`, `base$sdA` (additive SD of X_1 in t/ha) and
`base$h2Cycle1` = V_G(X_1)/(V_G(X_1) + σ²_e).

**A11. The sugar trait.**
Paper: "The locus also carries a sugar trait ... (additive effect 5.82 and dominance deviation −5.82,
residual standard deviation 0.46 ...). Only m_i is read, so these parameters affect no result."
Engine: not simulated. As a second AlphaSimR trait it would add residual draws to every population
and so shift every later draw while changing nothing that is read. The parameters are kept in the
configuration (`cfg$sugar`) and `locusA` (see E5) records the additive effect.

**A12. Seeds of founders and programme.**
Paper: "Its independent unit is the replicate, one founder object and one programme seed shared by
every arm and scenario. Ten replicates were run (seeds 4001 to 4010)"; Table 1, Experiment 1
"20260925; programme seed 7".
Engine: the founder object of replicate s is built after `set.seed(s)`. The programme seed is
s + `seeds$programmeOffset` (offset 0 by default, so "one seed per replicate") and is set with an
explicit generator kind (Mersenne-Twister, Inversion, Rejection) immediately before cycle 1, after
the founder object exists. A programme run therefore does not depend on whether the founders were
just built or read from disk (tested). With offset 0, the founder build and the programme start
from the same seed value; a non-zero offset separates the two streams if the authors prefer.
Experiment 1 uses founder seed 20260925 and programme seed 7.

**A13. The saved founder object.**
Paper: "The application driver cached the founder population of each replicate to disk on first use
and reloaded it for every arm and scenario. Each object is about 1.2 MB and records the seed".
Engine: `founder_<seed>.rds` holds the seed, the software versions, the founder part of the
configuration, the `SimParam`, X_1 (with phenotypes), the donor pool and the quantities in A10.
Each file gets a SHA-256 sidecar and the directory a `SHA256SUMS` list. On reload the seed, the
founder configuration and the checksum are checked. Every run works on a deep copy of the
`SimParam` (an R6 reference object), so runs never share state. An application-size file is about
0.5 MB (519 kB for seed 4001), smaller than the 1.2 MB quoted in Section 2.8 for the earlier objects.

## B. The cycle (`cycle.R`, `prediction.R`)

**B1. The training archive before cycle 1.**
Paper: "The archive holds the most recent 1,000 records. Each cycle adds 250 candidates drawn at
random from the elite nursery".
Engine: the archive is empty when cycle 1 starts; burn-in phenotypes are not used for training.
Cycle t therefore trains on min(250 t, 1,000) records (250, 500, 750, then 1,000 from cycle 4).

**B2. Training sample.**
Drawn with `sample.int` without replacement from X_t; the phenotypes are those drawn when X_t was
created. When the archive overflows, the oldest records are dropped.

**B3. Ridge regression.**
Paper: equation (4), rrBLUP, "markers coded X − 1".
Engine: `rrBLUP::mixed.solve(y, Z = M − 1)` with its defaults (REML, K = I, intercept as the only
fixed effect, default bounds on the ridge parameter). A marker monomorphic in the archive is
absorbed by the intercept and receives an effect of zero (to 10⁻¹⁴), consistent with Section 4.6.

**B4. Ranking ties.**
The parents are the first n_p of the GEBV ranking. Ties keep nursery order (stable radix sort); the
same rule ranks bridge parents and admission candidates.

**B5. Order of the mating pool.**
Paper: equation (7), "random mating within the union of the selected parents and the admitted lines".
Engine: the parents in GEBV-rank order followed by the admitted lines in admission order, mated by
AlphaSimR `randCross` (random distinct pairs, no selfing). The order fixes which individuals the
random pairs map onto; it is the same in every arm, which is what makes ELITE2M with an all-zero
schedule reproduce CLOSED exactly.

**B6. Recombination.**
Paper: "Mendelian segregation with recombination along a uniform genetic map of two morgans per
chromosome".
Engine: the MaCS maize map (2 M per chromosome, positions proportional to physical position) with
AlphaSimR's default crossover model (gamma model with interference parameter v = 2.6 and p = 0,
approximating the Kosambi map). The paper does not mention crossover interference; if it was meant
to be absent (Haldane), `SP$v` should be set to 1 and the Methods should say so.

**B7. Phenotypes of non-training populations.**
Paper: Algorithm 1, "[R] cross them and self each F1 (meiosis and phenotypes)".
Engine: every population created (F1, F2, BRIDGE's second-step F1 and progeny, MAS seedlings, the
new nursery) receives AlphaSimR phenotypes at creation, which consume the generator; donor-derived
phenotypes are never used ("Donor-derived candidates are never phenotyped before admission").

## C. Arms (`arms.R`)

**C1. Order of the draws at step 4.**
Paper: "draw six donor founders and six mates from π_t, each with replacement".
Engine: the six donors are drawn first, then the six mates; each donor × mate cross gives one F1
(`makeCross2`, donor as female), and each F1 is selfed into 150 F2. Which parent is female has no
genetic consequence (one genetic map).

**C2. BRIDGE when fewer than 20 F2 are sweet.**
Paper: "BRIDGE then crosses its 20 sweet F2 of highest ĝ to mates drawn with replacement from π_t,
selfs them and keeps the sweet homozygotes"; Table 1 "BRIDGE second step, families × progeny | 20 × 45".
Engine: one F1 per bridge parent, selfed into 45 lines. If fewer than 20 first-step F2 are sweet
(possible under resource matching: P = 0.14 per cycle with 96 F2), all sweet F2 are used and the
second step has fewer families. If none is sweet, nothing is offered.

**C3. MAS.**
Paper: "screens four times as many seedlings by marker and keeps the sweet homozygotes".
Engine: 4 × 150 = 600 seedlings per F1 (132 under resource matching); every sweet homozygote is grown
and offered. No cap applies to the number offered, only to the number admitted.

**C4. Admission threshold.**
Equation (6) is applied with "≥" (a candidate exactly at b_t is eligible), largest GEBV first.

**C5. Elite controls.**
Paper: "admit the next-ranked elite candidates below the parent cut, two in ELITE2 and the scheduled
number in ELITE2M"; "For ELITE2, whose offered lines cannot clear b_t, S_t and S*_t are counted over
the whole elite nursery"; Table 5 lists 600 candidates for ELITE2 and no accuracy or bias.
Engine: an elite control admits the candidates ranked n_p + 1, ..., n_p + k_t on GEBV. ELITE2 is the
elite control with k_t = κ = 2 in every cycle; ELITE2M uses its inflow arm's `nAdmitted` vector.
For both, `nCand` is the nursery size and `nClearBar`/`nClearTrue` count the whole nursery;
`accDonor` and `biasDonor` are NA (no donor-derived pool); `admitTrueGap` is computed for the lines
they admit. CLOSED records zeros for the counts and NA for the three pool statistics.

**C6. ELITE2M schedules and scenarios.**
Paper: "The ELITE2M arms replay the admissions recorded in the main run".
Engine: `run_replicate()` runs the inflow arm first and passes its `nAdmitted` vector to its ELITE2M
arm in the same call (same founders, programme seed and scenario), then checks the replay. By
default ELITE2M arms are run in `baseline` and `resource_matched` only, the two scenarios the
analysis uses; under visual sorting they would repeat those exactly (tested).
`--elite2m-scenarios all` runs them in all four.

**C7. CLOSED and ELITE2 in every scenario.**
They are run again in every scenario, not copied, so that the identities "visual sorting changes
nothing genetic" and "resource matching leaves ELITE2 unchanged" are checked rather than assumed.

## D. Scenarios and plot accounting (`outcomes.R`, `config.R`)

**D1. Resource-matched counts.**
Paper: "The 200 plots give 33 F2 per donor cross in DIRECT, 16 in the first step of BRIDGE, whose
two steps share them equally, and 132 screened seedlings per cross in MAS"; elite nursery 40 × 10.
Engine: derived from the budget and checked against Table 1 during validation:
DIRECT ⌊200/6⌋ = 33; BRIDGE step 1 ⌊100/6⌋ = 16 per cross, step 2 = 20 families × ⌊100/20⌋ = **5
progeny**; MAS 4 × 33 = 132. The step-2 size (20 × 5) is the engine's reading of "share them equally".

**D2. Smaller nursery from cycle 2.**
Paper: "in the first cycle their introgression plots are added to the common 600-plot nursery".
Engine: X_1 (600 plants) is common to every arm; a resource-matched inflow arm and its ELITE2M grow
40 × 10 = 400 plants from the nursery produced in cycle 1 onwards.
**[check]** Resulting per-cycle plots from cycle 2: DIRECT 400 + 198 = 598; MAS 400 + its sweet
seedlings (≈ 198); BRIDGE 400 + 96 + 100 = 596 (fewer when it has under 20 sweet F2). Over 20
cycles DIRECT and MAS grow about 1.3% more than the control, but BRIDGE about 1.0% more (at most
796 + 19 × 596 = 12,120 against 12,000). Section 2.2 says "about 598 plots" and "1.2 to 1.3% more".
Either BRIDGE's second step was sized differently in the original runs or the sentence should read
"1.0 to 1.3%".

**D3. What a plot is.**
Paper: "baseline DIRECT 1,500 plots per cycle, of which 900 introgression; BRIDGE 2,400 of which
1,800; under visual sorting 827 and 1,726; MAS 1,502".
Engine: `plotsGrown` = plants of the elite nursery X_t + plants of that cycle's introgression
nurseries. F1 plants (6 per cycle, and BRIDGE's ≤ 20 second-step F1) are not counted, which is what
reproduces 1,500 and 2,400. Under field discard every F2 and every second-step line is a plot;
under visual sorting only sweet homozygotes are; MAS grows only its sweet seedlings in both. The
expected values are 600 + 900/4 = 825 (DIRECT, visual), 600 + 900/4 + 900 = 1,725 (BRIDGE, visual)
and 600 + 3,600/4 = 1,500 (MAS); the paper's 827, 1,726 and 1,502 are realised means.

**D4. Visual sorting.**
Paper: "Under visual sorting they are removed as seed"; "with the number of kernels screened fixed
the assumption changes only which plants are grown".
Engine: a flag read only by the plot count (`stage_plots()`); no draw and no decision depends on it.
`visual_matched` is resource matching plus visual sorting.

## E. Records (`outcomes.R`)

**E1. When columns are measured.**
Population columns (`meanG` ... `pctSweet`, `gapCycle`) are measured at step 3 on X_t. The admission
columns (`nAdmitted` ... `admitTrueGap`, `nF2screened`, `nSweetUsable`) describe steps 4–5 of the
same cycle, and `plotsGrown` the plots of the same cycle. `year` = 3 t.

**E2. `nProgeny`.**
DATA_DICTIONARY: "progeny per cross in the elite nursery (15; 10 in the inflow arms and their
ELITE2M controls under resource matching)".
Engine: the progeny per cross of the nursery X_t measured in that row. It is 15 at cycle 1 in every
arm and scenario (X_1 is common) and 10 from cycle 2 in the resource-matched inflow arms and their
ELITE2M arms, consistent with `plotsGrown`.

**E3. Variances and accuracy.**
`totalVarG` = population variance (divisor N) of the genetic values of X_t, in (t/ha)².
`genicVar` = equation (11) over the QTL, in (t/ha)²; `genicPct` = 100 × V_t/V(X_1).
`accuracy` = Pearson correlation of GEBV and genetic value over all individuals of X_t.
`accDonor` (correlation) and `biasDonor` (equation 10, internal scale) are computed over the offered
donor-derived pool A_t; NA when it has fewer than two members or no variance.
`admitTrueGap` = mean genetic value of the admitted lines − b*_t, internal scale; NA when none admitted.

**E4. `donorPct` (not defined in the paper).**
DATA_DICTIONARY: "donor-allele content of the nursery at the markers, %".
Engine: a marker is diagnostic when X_1 is fixed for one allele and the donor pool carries the
other. `donorPct` = 100 × Σ (frequency of the diagnostic alleles in X_t) / Σ (their frequency in the
donor pool). It is exactly 0 in CLOSED, ELITE2 and ELITE2M (no mutation, no inflow), 100 for the
donor pool itself and about 50 for unselected donor × elite F2. Other definitions (IBD tracking,
all markers) would give different numbers; the Methods should define the one used.

**E5. `locusA`.**
DATA_DICTIONARY: "locusA, h2, nAdmit, nDonorCross | parameters of the run".
Engine: the sugar trait's additive effect at the endosperm locus, 5.82 (see A11). `h2` = 0.30,
`nAdmit` = κ = 2, `nDonorCross` = 6.

**E6. Identification columns.** `rep` = seed − 4000 (application); `RVersion` = "4.5.1";
`AlphaSimR` = "2.1.0"; `ompThreads` = the value of `OMP_NUM_THREADS` ("1").

**E7. Output files.**
`run_application.R` writes one CSV per replicate, `raw_<seed>.csv`, with every arm and scenario
(RECORD_COLUMNS order, the DATA_DICTIONARY columns). With `--pipeline-layout` it also writes
`raw_w<seed>.csv` (CLOSED, ELITE2, DIRECT, BRIDGE, MAS) and `raw_m<seed>.csv` (ELITE2M), which the
deposit keeps in `out/application/` and `out/matched_controls/` for `R/01_pool_application.py` and
`R/04_analyse_matched_controls.R`.

## F. Experiments (`run_exp1.R`, `run_exp2.R`)

**F1. Experiment 1 generator layers.**
Paper: "runMacs ... called for 40 individuals on three chromosomes of 150 segregating sites";
"quickHaplo ... with the same settings".
Engine: `runMacs(nInd = 40, nChr = 3, segSites = 150, species = "MAIZE", split = 100, nThreads = 1)`
and `quickHaplo(nInd = 40, nChr = 3, segSites = 150)` (genetic length 1 M, its default), each after
`set.seed(20260925)`, five calls each in one session and one call in each of three Rscript
sessions. Species and split are not stated for this layer; the application's history is used.
Genome signature: "mean dosage|dosage sum of individual 1|number of sites" (mean dosage to 10
decimals).

**F2. Experiment 1 programme layer.**
Engine: one founder object built under seed 20260925 with the Experiment 1 column of Table 1 (120
simulated and 60 kept per pool, 3 burn-in cycles), written to the output directory and read back in
the same session for layer 4; CLOSED, DIRECT, BRIDGE and MAS run five times on the in-memory object
and once on the reloaded one, programme seed 7, four cycles. Signature: "final gain|final genicVar"
(t/ha and (t/ha)², 10 decimals). The paper's name `runProgram` refers to the earlier code; the
function here is `run_programme()`.

**F3. Experiment 2 seeds.**
Paper: "The unpaired design gives every arm of every replicate its own founder population and seed";
four shards, base seeds 500000–530000.
Engine: paired replicate r of a shard uses seed base + r for its founders and for every arm's
programme. Unpaired replicate r, arm k uses base + 100000 k + r with **k = 1, 2, 3, 4 for CLOSED,
DIRECT, BRIDGE, MAS** (1-based, so that an unpaired seed never equals a paired one), for both its
founders and its programme. The shard file has an extra column `seed`.

**F4. Experiment 2 outcomes.**
`gain` and `genicPct` at cycle 10; `accLate` = mean accuracy over cycles 6–10 (Table 3 note).
Only the baseline scenario is run. Founders are not saved (Section 2.4), so a rerun reproduces the
outputs only in distribution.

## G. Software

`OMP_NUM_THREADS = 1` is set before AlphaSimR is loaded, every `SimParam` gets `nThreads = 1` and
`runMacs` is called with `nThreads = 1`. R's reference BLAS is used; with an optimised or threaded
BLAS, `mixed.solve` may give results that differ in the last digits.
