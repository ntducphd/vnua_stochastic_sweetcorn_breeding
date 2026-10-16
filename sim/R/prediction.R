## =============================================================================
## prediction.R -- training archive and ridge-regression genomic prediction
##
## The archive is the second half of the Markov state (Section 2.1.1): marker genotypes
## and phenotypes of the training candidates, newest last, at most `capacity` records.
## The model is equation (4): GEBV = mu + sum_l u_l (X_l - 1), fitted by REML ridge
## regression with rrBLUP::mixed.solve (marker effects i.i.d., intercept as the only
## fixed effect). Nothing here draws a random number.
## =============================================================================

new_archive <- function(capacity) {
  list(M = NULL, y = numeric(0), capacity = as.integer(capacity))
}

## Append records and keep the latest `capacity` of them.
archive_add <- function(archive, M, y) {
  stopifnot(nrow(M) == length(y))
  M <- rbind(archive$M, M)
  y <- c(archive$y, y)
  n <- length(y)
  if (n > archive$capacity) {
    keep <- seq.int(n - archive$capacity + 1L, n)
    M <- M[keep, , drop = FALSE]
    y <- y[keep]
  }
  archive$M <- M
  archive$y <- y
  archive
}

archive_size <- function(archive) length(archive$y)

## Fit equation (4) to the archive. Markers are coded X - 1 in {-1, 0, 1}.
fit_marker_model <- function(archive) {
  if (archive_size(archive) < 2) stop("the training archive holds fewer than two records")
  Z <- archive$M - 1
  storage.mode(Z) <- "double"
  fit <- rrBLUP::mixed.solve(y = archive$y, Z = Z)
  list(mu = as.numeric(fit$beta), u = as.numeric(fit$u), Vu = fit$Vu, Ve = fit$Ve,
       nRecords = archive_size(archive))
}

## GEBV of genotypes M (individuals x markers, coded 0/1/2).
predict_gebv <- function(model, M) {
  if (!nrow(M)) return(numeric(0))
  as.numeric(model$mu + (M - 1) %*% model$u)
}
