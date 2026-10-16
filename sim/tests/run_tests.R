#!/usr/bin/env Rscript
## =============================================================================
## run_tests.R -- tests of the simulation engine
##
## Uses two small configurations so that the whole suite runs in a few minutes:
##   tiny        config_test(): 2 chromosomes x 50 sites, 60-plant nursery, 4 cycles
##   accounting  the application's programme, introgression and resource-matching layers
##               (600-plot nursery, 6 x 150 F2, 20 x 45 bridge, 4x MAS screen) on a small
##               genome for 2 cycles, to check plot counts at their published values
## Usage: Rscript sim/tests/run_tests.R [--work <dir>]
## Writes sim/tests/test_results.txt and exits with status 1 if any test fails.
## =============================================================================

Sys.setenv(OMP_NUM_THREADS = "1")
TEST_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")])
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath("sim/tests")
})
SIM_DIR <- dirname(TEST_DIR)
source(file.path(SIM_DIR, "R", "load.R"))
load_sim(SIM_DIR)

opt <- parse_cli(commandArgs(trailingOnly = TRUE))
WORK <- if (is.null(opt$work)) file.path(tempdir(), "sim_tests") else opt$work
unlink(WORK, recursive = TRUE)
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
RSCRIPT <- file.path(R.home("bin"), "Rscript")

## ---- a minimal harness -----------------------------------------------------------
RESULTS <- list()
DETAILS <- character(0)
check <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)
note <- function(...) DETAILS <<- c(DETAILS, sprintf(...))
run_test <- function(name, code) {
  DETAILS <<- character(0)
  t0 <- Sys.time()
  err <- tryCatch({ force(code); NULL }, error = function(e) conditionMessage(e))
  s <- seconds_since(t0)
  status <- if (is.null(err)) "PASS" else "FAIL"
  RESULTS[[length(RESULTS) + 1]] <<- list(name = name, status = status, seconds = s,
                                          details = c(DETAILS, if (!is.null(err)) paste("error:", err)))
  cat(sprintf("[%s] %s (%.1f s)\n", status, name, s))
  for (d in RESULTS[[length(RESULTS)]]$details) cat("       ", d, "\n")
}
same_genetics <- function(a, b, cols = GENETIC_COLUMNS) {
  a <- a[order(a$cycle), cols]; b <- b[order(b$cycle), cols]
  rownames(a) <- rownames(b) <- NULL
  identical(a, b)
}
pick <- function(d, scenario, strategy) {
  z <- d[d$scenario == scenario & d$strategy == strategy, ]
  z[order(z$cycle), ]
}

## ---- configurations -----------------------------------------------------------------
cfgT <- config_test()
cfgA <- modify_config(config_application(), list(
  name = "accounting",
  genome = list(nChr = 2L, segSites = 60L, nQtlPerChr = 15L, nSnpPerChr = 15L),
  founders = list(nSimulatedPerPool = 80L, nKeptPerPool = 40L, burnInCycles = 1L),
  programme = list(nCycles = 2L)))
SEED_T <- 11L
SEED_A <- 21L
FDIR <- file.path(WORK, "founders")

## =============================================================================
run_test("configuration: Table 1 sets validate, bad sets are rejected", {
  app <- config_application(); e1 <- config_exp1(); e2 <- config_exp2()
  check(app$programme$nCrosses * app$programme$nProgeny == 600, "application nursery is not 600 plots")
  check(identical(derive_resource_matching(200L, 6L, 20L, 4L),
                  list(directF2PerCross = 33L, bridgeF2PerCross = 16L, bridgeProgeny = 5L,
                       masSeedlingsPerCross = 132L)), "resource-matched counts differ from Table 1")
  check(e1$programme$nCycles == 4 && e1$founders$burnInCycles == 3 && e1$seeds$programmeSeed == 7,
        "Experiment 1 settings")
  check(e2$genome$nChr == 4 && e2$programme$nCrosses * e2$programme$nProgeny == 200 &&
          e2$programme$nTrain == 150, "Experiment 2 settings")
  bad <- list(list(genome = list(nQtlPerChr = 150L, nSnpPerChr = 100L)),
              list(founders = list(nKeptPerPool = 500L)),
              list(resourceMatching = list(directF2PerCross = 34L)),
              list(programme = list(nTrain = 700L)),
              list(founders = list(h2 = 1.2)))
  for (b in bad) check(inherits(try(modify_config(app, b), silent = TRUE), "try-error"),
                       paste("accepted an invalid change:", paste(names(unlist(b)), collapse = ",")))
  note("application, exp1, exp2, test and accounting configurations valid; 5 invalid variants rejected")
})

run_test("training archive keeps the latest records; ridge model predicts", {
  a <- new_archive(100)
  for (k in 1:3) a <- archive_add(a, matrix(k, 40, 5) + matrix(rep(1:40, 5), 40), (k - 1) * 40 + 1:40)
  check(archive_size(a) == 100, "archive size is not the capacity")
  check(identical(a$y, as.numeric(21:120)), "archive did not keep the latest 100 records in order")
  seed_rng(1)
  M <- matrix(sample(0:2, 300 * 50, TRUE), 300); M[, 50] <- 0L       # one monomorphic marker
  u <- rnorm(50); u[50] <- 0
  y <- as.numeric((M - 1) %*% u + rnorm(300))
  arc <- archive_add(new_archive(1000), M, y)
  fit <- fit_marker_model(arc)
  r <- cor(predict_gebv(fit, M), (M - 1) %*% u)
  check(r > 0.8, sprintf("prediction accuracy %.2f on a simple simulated set", r))
  check(abs(fit$u[50]) < 1e-8, "a marker monomorphic in the archive received a non-zero effect")
  note("accuracy %.3f; effect of monomorphic marker %.1e", r, fit$u[50])
})

run_test("admission rule, equation (6): threshold, cap, largest first", {
  g <- c(5, 1, 7, 3, 7, 2)
  check(identical(admit_threshold_cap(g, 3, 2), c(3L, 5L)), "cap 2")
  check(identical(admit_threshold_cap(g, 3, 10), c(3L, 5L, 1L, 4L)), "no binding cap")
  check(identical(admit_threshold_cap(g, 10, 2), integer(0)), "nothing clears")
  check(identical(admit_threshold_cap(g, 3, 0), integer(0)), "cap 0")
})

## =============================================================================
founderT <- NULL
run_test("founders: locus, cull, disjoint sets, scale, residual variance, save/reload", {
  t0 <- Sys.time()
  founderT <<- load_or_build_founders(cfgT, SEED_T, FDIR, verbose = FALSE)
  note("built in %.0f s: locus %s, %d redraws, cull elite %.1f%% donor %.1f%%, max freq diff %.3f",
       seconds_since(t0), founderT$locus$name, founderT$nRedraw, 100 * founderT$cull$cullElite,
       100 * founderT$cull$cullDonor, founderT$locus$maxFreqDiff)
  f <- founderT; sp <- f$SP
  check(identical(attr(f, "source"), "built"), "first call did not build")
  check(all(locus_dosage(f$nursery, f, sp) == 2L), "post-burn-in elite nursery not all sweet")
  check(all(locus_dosage(f$donors, f, sp) == 0L), "donor founders not all homozygous non-sweet")
  check(nInd(f$donors) == cfgT$founders$nKeptPerPool, "wrong number of donor founders kept")
  check(!(f$locus$name %in% c(f$qtl, f$snp)), "locus inside the QTL or marker set")
  check(!length(intersect(f$qtl, f$snp)), "QTL and markers overlap")
  check(length(f$qtl) == cfgT$genome$nChr * cfgT$genome$nQtlPerChr &&
          length(f$snp) == cfgT$genome$nChr * cfgT$genome$nSnpPerChr, "set sizes")
  check(identical(f$qtl, loci_names(sp$traits[[1]])), "QTL names do not match the trait")
  g <- to_tha(f$eliteFounderGv, f$scale)
  check(isTRUE(all.equal(mean(g), 14)) && isTRUE(all.equal(pop_var(g), 1)),
        "kept elite founders are not at mean 14 and additive variance 1 t/ha")
  check(isTRUE(all.equal(f$scale$varGFounders, 1)), "combined founders not scaled to variance 1")
  check(isTRUE(all.equal(f$scale$varE, f$scale$varGFounders * 0.7 / 0.3)), "varE not from h2 = 0.30")
  check(isTRUE(all.equal(sp$varE, f$scale$varE)), "SimParam residual variance differs")
  check(isTRUE(verify_founder_checksums(FDIR)), "SHA256SUMS does not verify")
  again <- load_or_build_founders(cfgT, SEED_T, FDIR, verbose = FALSE)
  check(identical(attr(again, "source"), "file"), "second call did not reload the saved file")
  check(inherits(try(read_founders(founder_path(FDIR, SEED_T), cfgA, SEED_T), silent = TRUE), "try-error"),
        "a founder file was accepted under a different founder configuration")
  check(identical(f$software$AlphaSimR, as.character(packageVersion("AlphaSimR"))) &&
          f$seed == SEED_T, "seed or versions not recorded")
  note("varE %.4f (internal), base %.2f t/ha, donor %.2f t/ha, file %.0f kB",
       f$scale$varE, f$base$meanG, f$base$donorMean, file.size(founder_path(FDIR, SEED_T)) / 1024)
})

## =============================================================================
run_test("Mendelian segregation at the endosperm locus (selfed donor x elite F1)", {
  ## ten donor x elite F1 families, each selfed into 1,000 F2
  f <- founderT; sp <- f$SP$clone(deep = TRUE)
  seed_rng(99)
  plan <- cbind(sample.int(nInd(f$donors), 10, replace = TRUE), sample.int(nInd(f$nursery), 10))
  f1 <- makeCross2(f$donors, f$nursery, plan, nProgeny = 1, simParam = sp)
  f2 <- self(f1, nProgeny = 1000, simParam = sp)
  d <- locus_dosage(f2, f, sp)
  obs <- tabulate(d + 1L, 3); n <- sum(obs)
  p <- chisq.test(obs, p = c(0.25, 0.5, 0.25))$p.value
  check(all(locus_dosage(f1, f, sp) == 1L), "an F1 is not heterozygous")
  check(p > 0.001, sprintf("1:2:1 rejected (p = %.2g)", p))
  check(abs(obs[3] / n - 0.25) < 4 * sqrt(0.25 * 0.75 / n), "sweet fraction far from 1/4")
  note("10 F1 families x 1,000 F2: dosage 0/1/2 = %d/%d/%d, sweet %.4f, chi-square p = %.3f",
       obs[1], obs[2], obs[3], obs[3] / n, p)
})

## =============================================================================
tiny <- NULL
run_test("tiny replicate: all arms and scenarios run; elite pool always sweet", {
  t0 <- Sys.time()
  tiny <<- run_replicate(cfgT, SEED_T, founder = founderT, verbose = FALSE)
  check(nrow(tiny) == (4 * 5 + 2 * 3) * cfgT$programme$nCycles, "unexpected number of records")
  check(identical(names(tiny), RECORD_COLUMNS), "columns differ from the data dictionary")
  check(all(tiny$pctSweet == 100), "an elite nursery was not homozygous sweet")
  check(all(tiny$gain[tiny$cycle == 1] == 0), "gain at cycle 1 is not zero")
  check(all(tiny$donorPct[tiny$strategy %in% c("CLOSED", "ELITE2") | startsWith(tiny$strategy, "ELITE2M")] == 0),
        "donor-allele content non-zero in an arm without donors")
  note("%d records in %.1f s", nrow(tiny), seconds_since(t0))
})

run_test("ELITE2M identities: all-zero schedule = CLOSED, all-twos schedule = ELITE2", {
  pseed <- programme_seed(cfgT, SEED_T)
  nC <- cfgT$programme$nCycles
  closed <- run_programme(founderT, cfgT, "CLOSED", "baseline", pseed)
  elite2 <- run_programme(founderT, cfgT, "ELITE2", "baseline", pseed)
  cols <- setdiff(GENETIC_COLUMNS, c("nAdmitted", "nCand", "nClearBar", "nClearTrue", "admitTrueGap"))
  for (a in INFLOW_ARMS) {
    z <- run_programme(founderT, cfgT, paste0(MATCHED_PREFIX, a), "baseline", pseed, schedule = rep(0L, nC))
    check(same_genetics(z, closed, c(cols, "plotsGrown")), paste("zero schedule differs from CLOSED for", a))
    w <- run_programme(founderT, cfgT, paste0(MATCHED_PREFIX, a), "baseline", pseed, schedule = rep(2L, nC))
    check(same_genetics(w, elite2, c(GENETIC_COLUMNS, "plotsGrown")), paste("twos schedule differs from ELITE2 for", a))
  }
  for (sc in c("baseline", "resource_matched")) for (a in INFLOW_ARMS)
    check(identical(pick(tiny, sc, a)$nAdmitted, pick(tiny, sc, paste0(MATCHED_PREFIX, a))$nAdmitted),
          paste("ELITE2M did not replay", a, "in", sc))
  note("exact equality (identical()) of every genetic column over %d cycles; replay checked in 6 arm-scenarios", nC)
})

run_test("paired start: arms share cycle-1 training sample and parents; streams part at step 4 (inflow) or 6 (elite)", {
  pseed <- programme_seed(cfgT, SEED_T)
  tr <- lapply(setNames(MAIN_ARMS, MAIN_ARMS), function(a)
    run_programme(founderT, cfgT, a, "baseline", pseed, trace = TRUE))
  first <- lapply(tr, function(x) attr(x, "trace")[[1]])
  for (a in MAIN_ARMS[-1]) {
    check(identical(first[[a]]$trainIdx, first$CLOSED$trainIdx), paste(a, "training sample differs"))
    check(identical(first[[a]]$parents, first$CLOSED$parents), paste(a, "parents differ"))
    check(identical(first[[a]]$rngBeforeOffer, first$CLOSED$rngBeforeOffer), paste(a, "generator state differs before step 4"))
    pop1 <- c("meanG", "gain", "genicVar", "genicPct", "totalVarG", "accuracy", "donorPct", "pctSweet")
    check(identical(tr[[a]][1, pop1], tr$CLOSED[1, pop1]), paste(a, "cycle-1 outcomes differ"))
  }
  for (a in INFLOW_ARMS)
    check(!identical(first[[a]]$rngBeforeMating, first$CLOSED$rngBeforeMating), paste(a, "drew nothing at step 4"))
  check(identical(first$ELITE2$rngBeforeMating, first$CLOSED$rngBeforeMating), "ELITE2 drew at step 5")
  np <- cfgT$programme$nParents
  check(identical(first$ELITE2$admittedIds, founderT$nursery@id[first$CLOSED$rank[np + 1:2]]),
        "ELITE2 did not admit the two candidates ranked just below the parent cut")
  check(tr$ELITE2$meanG[2] != tr$CLOSED$meanG[2], "ELITE2 did not diverge from CLOSED at step 6 of cycle 1")
  note("training sample, parents and generator state identical at step 4 of cycle 1 in all five arms")
})

run_test("visual sorting is genetically identical to field discard; resource matching leaves CLOSED and ELITE2 unchanged", {
  for (a in MAIN_ARMS) {
    check(same_genetics(pick(tiny, "baseline", a), pick(tiny, "visual_sort", a)), paste(a, "baseline vs visual_sort"))
    check(same_genetics(pick(tiny, "resource_matched", a), pick(tiny, "visual_matched", a)),
          paste(a, "resource_matched vs visual_matched"))
  }
  ## ELITE2M under visual sorting (run on request with --elite2m-scenarios all)
  vs <- run_replicate(cfgT, SEED_T, founder = founderT, scenarios = c("baseline", "visual_sort"),
                      arms = c("CLOSED", "DIRECT"), elite2mScenarios = c("baseline", "visual_sort"), verbose = FALSE)
  check(same_genetics(pick(vs, "baseline", "ELITE2M_DIRECT"), pick(vs, "visual_sort", "ELITE2M_DIRECT"),
                      c(GENETIC_COLUMNS, "plotsGrown")), "ELITE2M_DIRECT differs under visual sorting")
  for (a in c("CLOSED", "ELITE2"))
    check(same_genetics(pick(tiny, "baseline", a), pick(tiny, "resource_matched", a), c(GENETIC_COLUMNS, "plotsGrown")),
          paste(a, "changed under resource matching"))
  note("all %d genetic columns identical, for every arm", length(GENETIC_COLUMNS))
})

run_test("a programme run from a saved founder object is identical in a new R session", {
  ## founders were built in this session; run from the in-memory object and from the file
  fromMemory <- run_replicate(cfgT, SEED_T, founder = founderT, verbose = FALSE)
  fromFile <- run_replicate(cfgT, SEED_T, founder = read_founders(founder_path(FDIR, SEED_T), cfgT, SEED_T),
                            verbose = FALSE)
  a <- file.path(WORK, "session_parent.csv"); write_records(fromFile, a)
  m <- file.path(WORK, "session_memory.csv"); write_records(fromMemory, m)
  cfgPath <- file.path(WORK, "cfg_tiny.rds"); saveRDS(cfgT, cfgPath)
  b <- file.path(WORK, "session_child.csv")
  status <- system2(RSCRIPT, c(shQuote(file.path(TEST_DIR, "rerun_in_new_session.R")), shQuote(SIM_DIR),
                               shQuote(cfgPath), shQuote(FDIR), SEED_T, shQuote(b)), stdout = TRUE, stderr = TRUE)
  check(file.exists(b), paste("child session failed:", paste(status, collapse = " | ")))
  h <- c(memory = sha256_file(m), parent = sha256_file(a), child = sha256_file(b))
  check(h[["memory"]] == h[["parent"]], "run from the in-memory founders differs from the run from the file")
  check(h[["parent"]] == h[["child"]], "records differ between R sessions")
  check(identical(fromFile, tiny[, names(fromFile)]) || isTRUE(all.equal(fromFile, tiny, check.attributes = FALSE)),
        "rerun in the same session differs")
  note("%d records; SHA-256 of the three CSVs %s", nrow(fromFile), substr(h[["child"]], 1, 16))
})

## =============================================================================
acct <- NULL
run_test("plot accounting at the application's programme sizes (all scenarios, 2 cycles)", {
  t0 <- Sys.time()
  fA <- load_or_build_founders(cfgA, SEED_A, FDIR, verbose = FALSE)
  acct <<- run_replicate(cfgA, SEED_A, founder = fA, verbose = FALSE)
  note("accounting founders + replicate in %.0f s", seconds_since(t0))
  check(isTRUE(verify_founder_checksums(FDIR)), "SHA256SUMS does not verify after a second founder file")
  P <- function(sc, a) pick(acct, sc, a)
  ## baseline (field discard)
  check(all(P("baseline", "DIRECT")$plotsGrown == 1500) && all(P("baseline", "DIRECT")$nF2screened == 900),
        "baseline DIRECT is not 1,500 plots with 900 introgression plots")
  check(all(P("baseline", "BRIDGE")$plotsGrown == 2400) && all(P("baseline", "BRIDGE")$nF2screened == 1800),
        "baseline BRIDGE is not 2,400 plots with 1,800 introgression plots")
  m <- P("baseline", "MAS")
  check(all(m$plotsGrown == 600 + m$nSweetUsable) && all(m$nF2screened == 3600), "baseline MAS accounting")
  for (sc in names(SCENARIOS)) for (a in c("CLOSED", "ELITE2"))
    check(all(P(sc, a)$plotsGrown == 600), paste(a, "does not grow 600 plots in", sc))
  ## visual sorting
  vd <- P("visual_sort", "DIRECT"); vb <- P("visual_sort", "BRIDGE"); vm <- P("visual_sort", "MAS")
  check(all(vd$plotsGrown == 600 + vd$nSweetUsable), "visual DIRECT is not 600 + sweet F2")
  check(all(vb$plotsGrown == 600 + vb$nSweetUsable) && all(vb$nSweetUsable >= 900),
        "visual BRIDGE is not 600 + sweet F2 + 900 second-step lines")
  check(identical(vm$plotsGrown, m$plotsGrown), "MAS plots differ between field discard and visual sorting")
  ## expected values against the published means (827, 1,726, 1,502)
  sdMean <- function(n, k) sqrt(n * 0.25 * 0.75 / k)
  k <- cfgA$programme$nCycles
  check(abs(mean(vd$plotsGrown) - 825) < 4 * sdMean(900, k), "visual DIRECT plots far from 600 + 900/4")
  check(abs(mean(vb$plotsGrown) - 1725) < 4 * sdMean(900, k), "visual BRIDGE plots far from 600 + 900/4 + 900")
  check(abs(mean(m$plotsGrown) - 1500) < 4 * sdMean(3600, k), "MAS plots far from 600 + 3600/4")
  ## resource matching: 200 introgression plots; nursery 600 in cycle 1, 400 afterwards
  rd <- P("resource_matched", "DIRECT"); rb <- P("resource_matched", "BRIDGE"); rmas <- P("resource_matched", "MAS")
  check(identical(rd$plotsGrown, c(798L, 598L)) && all(rd$nF2screened == 198), "matched DIRECT is not 798 then 598")
  check(all(rb$plotsGrown == c(600, 400) + rb$nF2screened) && all(rb$nF2screened <= 196) &&
          all((rb$nF2screened - 96) %% 5 == 0), "matched BRIDGE is not nursery + 96 F2 + 5 x families")
  check(all(rmas$plotsGrown == c(600, 400) + rmas$nSweetUsable) && all(rmas$nF2screened == 792), "matched MAS")
  for (a in INFLOW_ARMS) {
    check(identical(P("resource_matched", a)$nProgeny, c(15L, 10L)), paste("matched", a, "nursery progeny"))
    check(identical(P("resource_matched", paste0(MATCHED_PREFIX, a))$plotsGrown, c(600L, 400L)),
          paste("matched ELITE2M", a, "does not grow its arm's smaller nursery"))
  }
  vmd <- P("visual_matched", "DIRECT")
  check(all(vmd$plotsGrown == c(600, 400) + vmd$nSweetUsable), "visual_matched DIRECT")
  tab <- aggregate(plotsGrown ~ strategy + scenario, acct[!startsWith(acct$strategy, MATCHED_PREFIX), ], mean)
  for (sc in names(SCENARIOS)) note("%-16s %s", sc, paste(sprintf("%s %.0f", tab$strategy[tab$scenario == sc],
                                                                    tab$plotsGrown[tab$scenario == sc]), collapse = ", "))
})

run_test("recovery of sweet homozygotes: 1/4 for DIRECT and MAS F2, 0.625 for BRIDGE", {
  base <- acct[acct$scenario %in% c("baseline", "resource_matched"), ]
  for (a in c("DIRECT", "MAS")) {
    z <- base[base$strategy == a, ]
    bt <- binom.test(sum(z$nSweetUsable), sum(z$nF2screened), 0.25)
    check(bt$p.value > 0.001, sprintf("%s recovery %.4f rejects 1/4", a, bt$estimate))
    note("%s: %d of %d sweet (%.4f), exact binomial p = %.3f", a, sum(z$nSweetUsable), sum(z$nF2screened),
         bt$estimate, bt$p.value)
  }
  b <- pick(acct, "baseline", "BRIDGE")
  ratio <- sum(b$nSweetUsable) / sum(b$nF2screened)
  step1 <- binom.test(sum(b$nSweetUsable - 900), 900 * nrow(b), 0.25)     # second step: 900 lines, all sweet
  check(step1$p.value > 0.001, "BRIDGE first-step recovery rejects 1/4")
  check(abs(ratio - 0.625) < 4 * sqrt(900 * nrow(b) * 0.1875) / (1800 * nrow(b)), "BRIDGE recovery far from 0.625")
  note("BRIDGE: %d of %d sweet (%.4f; expected (900 x 0.25 + 900)/1800 = 0.625)", sum(b$nSweetUsable),
       sum(b$nF2screened), ratio)
  check(all(acct$pctSweet == 100), "an elite nursery of the accounting runs was not all sweet")
  for (a in MAIN_ARMS) check(same_genetics(pick(acct, "baseline", a), pick(acct, "visual_sort", a)),
                     paste(a, "differs between field discard and visual sorting at application sizes"))
})

## =============================================================================
run_test("drivers: run_application.R (resumable, pipeline layout), run_exp2.R, run_exp1.R", {
  out <- file.path(WORK, "app")
  args <- c(shQuote(file.path(SIM_DIR, "run_application.R")), "--config", "test", "--seed", "5",
            "--out", shQuote(out), "--founders", shQuote(file.path(WORK, "app_founders")), "--pipeline-layout")
  log1 <- system2(RSCRIPT, args, stdout = TRUE, stderr = TRUE)
  f <- replicate_path(out, 5)
  check(file.exists(f), paste("run_application.R wrote no CSV:", paste(utils::tail(log1, 3), collapse = " | ")))
  d <- read_records(f)
  check(identical(names(d), RECORD_COLUMNS) && nrow(d) == 26 * cfgT$programme$nCycles, "replicate CSV shape")
  w <- read.csv(file.path(out, "raw_w5.csv")); m <- read.csv(file.path(out, "raw_m5.csv"))
  check(nrow(w) == 20 * cfgT$programme$nCycles && nrow(m) == 6 * cfgT$programme$nCycles &&
          identical(names(w), names(m)), "pipeline layout")
  check(file.exists(file.path(out, "timing_5.tsv")), "no timing log")
  check(file.exists(file.path(WORK, "app_founders", "founder_5.rds")) &&
          isTRUE(verify_founder_checksums(file.path(WORK, "app_founders"))), "founder file or checksum missing")
  mt <- file.mtime(f)
  log2 <- system2(RSCRIPT, args, stdout = TRUE, stderr = TRUE)
  check(any(grepl("skipped", log2)) && identical(file.mtime(f), mt), "an existing replicate was not skipped")
  note("run_application.R: %d records, rerun skipped", nrow(d))

  e2 <- file.path(WORK, "exp2")
  a2 <- c(shQuote(file.path(SIM_DIR, "run_exp2.R")), "--config", "test", "--shard", "1", "--reps", "1",
          "--out", shQuote(e2))
  log3 <- system2(RSCRIPT, a2, stdout = TRUE, stderr = TRUE)
  s <- file.path(e2, "shard1.csv")
  check(file.exists(s), paste("run_exp2.R wrote no shard:", paste(utils::tail(log3, 3), collapse = " | ")))
  sh <- read.csv(s)
  check(nrow(sh) == 8 && all(c("design", "seedBase", "rep", "arm", "gain", "genicPct", "accLate") %in% names(sh)),
        "shard layout")
  check(all(sh$seedBase == 510000) && isTRUE(all.equal(as.numeric(sort(unique(sh$seed[sh$design == "unpaired"]))),
                                                         510001 + 100000 * 1:4)), "shard seeds")
  check(length(unique(sh$seed[sh$design == "paired"])) == 1, "paired arms do not share one seed")
  log4 <- system2(RSCRIPT, a2, stdout = TRUE, stderr = TRUE)
  check(nrow(read.csv(s)) == 8, "a completed Experiment 2 block was rerun")
  note("run_exp2.R: 8 rows (2 designs x 4 arms), rerun skipped")

  e1 <- file.path(WORK, "exp1")
  log5 <- system2(RSCRIPT, c(shQuote(file.path(SIM_DIR, "run_exp1.R")), "all", "--config", "test",
                             "--out", shQuote(e1)), stdout = TRUE, stderr = TRUE)
  ws <- file.path(e1, "within_session.csv")
  check(file.exists(ws), paste("run_exp1.R failed:", paste(utils::tail(log5, 3), collapse = " | ")))
  x <- read.csv(ws, stringsAsFactors = FALSE)
  t2 <- read.csv(file.path(e1, "table2.csv"), stringsAsFactors = FALSE)
  prog <- t2[startsWith(t2$layer, "programme"), ]
  check(nrow(x) == 2 * 2 + 2 * 3 && all(prog$reproducible), "programme layer not reproducible")
  check(all(t2$reproducible[grepl("quickHaplo", t2$layer)]), "quickHaplo not reproducible")
  ac <- read.csv(file.path(e1, "across_sessions.csv"))
  check(nrow(ac) == 4 && length(unique(ac$pid)) == 2, "across-session layout")
  note("run_exp1.R: %s", paste(sprintf("%s: %d/%d distinct", t2$layer, t2$distinct, t2$calls), collapse = "; "))
})

## ---- report ----------------------------------------------------------------------
nFail <- sum(vapply(RESULTS, function(r) r$status == "FAIL", logical(1)))
lines <- c(sprintf("sim/tests/run_tests.R  %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
           sprintf("%s | AlphaSimR %s | rrBLUP %s | %s | OMP_NUM_THREADS=%s", R.version.string,
                   packageVersion("AlphaSimR"), packageVersion("rrBLUP"), R.version$platform,
                   Sys.getenv("OMP_NUM_THREADS")),
           "",
           unlist(lapply(RESULTS, function(r) c(sprintf("[%s] %s (%.1f s)", r$status, r$name, r$seconds),
                                                if (length(r$details)) paste("       ", r$details)))),
           "",
           sprintf("%d tests, %d passed, %d failed, %.1f min", length(RESULTS), length(RESULTS) - nFail, nFail,
                   sum(vapply(RESULTS, `[[`, numeric(1), "seconds")) / 60))
writeLines(lines, file.path(TEST_DIR, "test_results.txt"))
cat("\n", utils::tail(lines, 1), "\n", sep = "")
quit(status = if (nFail) 1 else 0)
