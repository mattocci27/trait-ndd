library(targets)
library(tarchetypes)
library(tidyverse)
library(stantargets)
library(cmdstanr)
library(furrr)

source("R/data_clean.R")
#source("R/stan.R")

plan(multicore)
options(clustermq.scheduler = "multicore")

tar_option_set(packages = c(
  "tidyverse",
  "patchwork",
  "parallel",
  "janitor",
  "extrafont",
  "loo",
  "jsonlite",
  "doParallel",
  "foreach",
  "httpgd",
  "multcompView"
))

# check if it's inside a container
# if (file.exists("/.dockerenv") | file.exists("/.singularity.d/startscript")) {
#   Sys.setenv(CMDSTAN = "/opt/cmdstan/cmdstan-2.29.0")
#   set_cmdstan_path("/opt/cmdstan/cmdstan-2.29.0")
# }

cmdstan_version()

list(
  # data cleaning ----------------------------------
  tar_target(
    seedling_csv,
    "data-raw/seedlingmatrix.csv",
    format = "file"
  ),
  tar_target(
    habitat_csv,
    "data-raw/habitat150.csv",
    format = "file"
  ),
  tar_target(
    trait_csv,
    "data-raw/trait_sp_code.csv",
    format = "file"
  ),
  tar_target(
    data_list,
    gen_seedling(seedling_csv, trait_csv, habitat_csv, n_ab = 50)
  ),


  # stan -------------------------------------------------
  tar_target(
    dry_full,
    gen_stan_dat(data_list, season = "dry", habitat = "all",
      neighbor = "full", trait_set = "full"),
  ),
  tar_target(
    wet_full,
    gen_stan_dat(data_list, season = "rainy", habitat = "all",
      neighbor = "full", trait_set = "full"),
  ),
  tar_target(
    dry_wd,
    gen_stan_dat(data_list, season = "dry", habitat = "all",
      neighbor = "full", trait_set = "wd"),
  ),
  tar_target(
    wet_wd,
    gen_stan_dat(data_list, season = "rainy", habitat = "all",
      neighbor = "full", trait_set = "wd"),
  ),
  tar_target(
    dry_pca,
    gen_stan_dat(data_list, season = "dry", habitat = "all",
      neighbor = "full", trait_set = "pca"),
  ),
  tar_target(
    wet_pca,
    gen_stan_dat(data_list, season = "rainy", habitat = "all",
      neighbor = "full", trait_set = "pca"),
  )

  # tar_target(
  #   pa_stan_dat,
  #   generate_pa_stan(read_csv(pa_csv)),
  # ),
  # tar_stan_mcmc(
  #   fit_1,
  #   "stan/GL_LMA.stan",
  #   data = gl_stan_dat,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 2000,
  #   iter_sampling = 2000,
  #   adapt_delta = 0.9,
  #   max_treedepth = 15,
  #   seed = 123),
  # tar_stan_mcmc(
  #   fit_2,
  #   "stan/PA_LMA.stan",
  #   data = pa_stan_dat,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 2000,
  #   iter_sampling = 2000,
  #   adapt_delta = 0.9,
  #   max_treedepth = 15,
  #   seed = 123)

 #)#,
  # tar_render(
  #   report,
  #   "report.Rmd"
  # )
)
