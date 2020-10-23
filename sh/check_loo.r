library(tidyverse)
library(rstan)
library(bayesplot)
options(mc.cores = parallel::detectCores())

argv <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(argv[1])
print(n)

files <- list.files("./rda")

n_samp <- NULL
res_loo <- NULL
file_name <- str_c("./rda/", files[n])
file_name
load(file_name)
tmp <- loo(fit)$estimates[3, 1]
res_loo <- rbind(res_loo,
    c(file_name = files[n], 
      season = dry,
      hab = hab,
      trait = trait_data,
      n_samp = list_dat_d$N, 
      looic = tmp)) %>%
  as_tibble



# R values
output <- "loo.csv"
#out <- file(paste(output), "w") # write

#writeLines(paste0("file: ", res_loo$file_name),
#           out,
#           sep = "\n")
#writeLines(paste0(  "season: ", res_loo$season),
#           out,
#           sep = "\n")
#writeLines(paste0(  "hab: ", res_loo$hab),
#           out,
#           sep = "\n")
#writeLines(paste0(  "trait: ", res_loo$trait),
#           out,
#           sep = "\n")
#writeLines(paste0(  "n_samp: ", res_loo$n_samp),
#           out,
#           sep = "\n")
#writeLines(paste0(  "looic: ", res_loo$looic),
#           out,
#           sep = "\n")
cat(paste(res_loo$file_name, 
          res_loo$season,
          res_loo$hab,
          res_loo$trait,
          res_loo$n_samp,
          tmp, sep = ","), file = output, append = TRUE, sep = "\n")
#writeLines(paste0("file: ", file_name),
#           out,
#           sep = "\n")
#writeLines(paste0("     n: ", n),
#           out,
#           sep = "\n")
#close(out)
