argv <- commandArgs(trailingOnly = TRUE)
workers <- as.numeric(argv[1])

if (workers == 1) {
  targets::tar_make()
} else if (workers == 100) {
  targets::tar_make_clustermq(
   workers = parallel::detectCores())
} else if (is.numeric(workers)) {
# for HPC
  targets::tar_make_clustermq(
   workers = workers)
}
