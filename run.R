arg = commandArgs(T)

targets::tar_make_clustermq(
  workers = as.numeric(arg[1])
 # workers = parallel::detectCores()#,
# reporter = "silent"
)

# targets::tar_make(
# reporter = "silent"
# )
