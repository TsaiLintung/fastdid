# Run the tinytest suite with one test file per worker, from the package root:
#
#   Rscript development/run_tests.R [workers]
#
# Exits 1 when a test fails. dev/deploy-fastdid.sh in the mono repo and its CI
# call it. The files run in any order, so a test file must not need another.

args <- commandArgs(trailingOnly = TRUE)
cores <- parallel::detectCores()
n_files <- length(list.files("inst/tinytest", "^test.*\\.[rR]$"))
workers <- if (length(args)) as.integer(args[1]) else min(6L, n_files, max(1L, cores - 1L))
# the workers share the cores, so each worker takes its part of the threads
threads <- max(1L, cores %/% workers)

cl <- parallel::makeCluster(workers)
parallel::clusterExport(cl, "threads")
invisible(parallel::clusterEvalQ(cl, {
  suppressMessages(devtools::load_all(quiet = TRUE))
  data.table::setDTthreads(threads)
  NULL
}))

suppressMessages(devtools::load_all(quiet = TRUE))
r <- tinytest::run_test_dir("inst/tinytest", cluster = cl, verbose = 0)
parallel::stopCluster(cl)
print(summary(r))
if (any(!as.logical(r))) quit(status = 1)
