################################################################################
## session_info.R -- record the computing environment in environment/.
##
## Writes environment/sessionInfo.txt (R, packages, BLAS and LAPACK, platform, threads) and
## environment/R_packages.csv (name and version of every package the scripts load). The
## Python versions of the deposited run are in environment/python_versions.txt.
## Usage: Rscript scripts/session_info.R
################################################################################

Sys.setenv(OMP_NUM_THREADS = "1")
suppressMessages({ library(AlphaSimR); library(rrBLUP); library(jsonlite) })
dir.create("environment", showWarnings = FALSE)

sink("environment/sessionInfo.txt")
cat("Recorded", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), "\n\n")
print(sessionInfo())
cat("\nextSoftVersion():\n"); print(extSoftVersion())
cat("\nLa_library(): ", La_library(), "\n", sep = "")
cat("La_version(): ", La_version(), "\n", sep = "")
cat("BLAS (sessionInfo()$BLAS): ", sessionInfo()$BLAS, "\n", sep = "")
cat("OMP_NUM_THREADS = ", Sys.getenv("OMP_NUM_THREADS"), "\n", sep = "")
sink()

pk <- c("AlphaSimR", "rrBLUP", "jsonlite", "Rcpp", "RcppArmadillo", "R6")
v  <- vapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_),
             character(1))
built <- vapply(pk, function(p) tryCatch(packageDescription(p)$Built, error = function(e) NA_character_),
                character(1))
write.csv(data.frame(package = c("R", pk), version = c(paste(R.version$major, R.version$minor, sep = "."), v),
                     built = c(R.version.string, built)),
          "environment/R_packages.csv", row.names = FALSE)
cat("wrote environment/sessionInfo.txt and environment/R_packages.csv\n")
