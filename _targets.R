library(targets)
library(tarchetypes)
library(tidyverse)
library(stantargets)
library(cmdstanr)
library(furrr)
library(languageserver)

source("R/data_clean.R")
source("R/stan.R")
source("R/figs.R")

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

tar_option_set(
  garbage_collection = TRUE,
  memory = "transient"
)

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
    dry_full_int,
    gen_stan_dat(data_list, season = "dry", habitat = "all",
      inter = TRUE, trait_set = "full"),
  ),
  tar_target(
    wet_full_int,
    gen_stan_dat(data_list, season = "rainy", habitat = "all",
      inter = TRUE, trait_set = "full"),
  ),
  tar_target(
    dry_wd_int,
    gen_stan_dat(data_list, season = "dry", habitat = "all",
      inter = TRUE, trait_set = "wd"),
  ),
  tar_target(
    wet_wd_int,
    gen_stan_dat(data_list, season = "rainy", habitat = "all",
      inter = TRUE, trait_set = "wd"),
  ),
  tar_target(
    dry_pca_int,
    gen_stan_dat(data_list, season = "dry", habitat = "all",
      inter = TRUE, trait_set = "pca"),
  ),
  tar_target(
    wet_pca_int,
    gen_stan_dat(data_list, season = "rainy", habitat = "all",
      inter = TRUE, trait_set = "pca"),
  ),
  tar_target(
    dry_full_noint,
    gen_stan_dat(data_list, season = "dry", habitat = "all",
      inter = FALSE, trait_set = "full"),
  ),
  tar_target(
    wet_full_noint,
    gen_stan_dat(data_list, season = "rainy", habitat = "all",
      inter = FALSE, trait_set = "full"),
  ),
  tar_target(
    dry_wd_noint,
    gen_stan_dat(data_list, season = "dry", habitat = "all",
      inter = FALSE, trait_set = "wd"),
  ),
  tar_target(
    wet_wd_noint,
    gen_stan_dat(data_list, season = "rainy", habitat = "all",
      inter = FALSE, trait_set = "wd"),
  ),
  tar_target(
    dry_pca_noint,
    gen_stan_dat(data_list, season = "dry", habitat = "all",
      inter = FALSE, trait_set = "pca"),
  ),
  tar_target(
    wet_pca_noint,
    gen_stan_dat(data_list, season = "rainy", habitat = "all",
      inter = FALSE, trait_set = "pca"),
  ),

  tar_stan_mcmc(
    fit_1_dry_full_int,
    "stan/model_ind.stan",
    data = dry_full_int,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 2000,
    iter_sampling = 1000,
    draws = FALSE,
    diagnostics = FALSE,
    summary = FALSE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_2_wet_full_int,
    "stan/model_ind.stan",
    data = wet_full_int,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 2000,
    iter_sampling = 1000,
    draws = FALSE,
    diagnostics = FALSE,
    summary = FALSE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_3_dry_wd_int,
    "stan/model_ind.stan",
    data = dry_wd_int,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 2000,
    iter_sampling = 1000,
    draws = FALSE,
    diagnostics = TRUE,
    summary = TRUE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_4_wet_wd_int,
    "stan/model_ind.stan",
    data = wet_wd_int,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 2000,
    iter_sampling = 1000,
    draws = FALSE,
    diagnostics = TRUE,
    summary = TRUE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    test_dry,
    "stan/model_no_like.stan",
    data = dry_wd_int,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 1000,
    iter_sampling = 1000,
    draws = FALSE,
    diagnostics = FALSE,
    summary = FALSE,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    test_wet,
    "stan/model_no_like.stan",
    data = wet_wd_int,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 1000,
    iter_sampling = 1000,
    draws = FALSE,
    diagnostics = FALSE,
    summary = FALSE,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123)#,
  # tar_stan_mcmc(
  #   fit_5_dry_pca_int,
  #   "stan/model_ind.stan",
  #   data = dry_pca_int,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 2000,
  #   iter_sampling = 1000,
  #   draws = FALSE,
  #   diagnostics = FALSE,
  #   summary = FALSE,
  #   adapt_delta = 0.95,
  #   max_treedepth = 15,
  #   seed = 123),
  # tar_stan_mcmc(
  #   fit_6_wet_pca_int,
  #   "stan/model_ind.stan",
  #   data = wet_pca_int,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 2000,
  #   iter_sampling = 1000,
  #   draws = FALSE,
  #   diagnostics = FALSE,
  #   summary = FALSE,
  #   adapt_delta = 0.95,
  #   max_treedepth = 15,
  #   seed = 123),
  # tar_stan_mcmc(
  #   fit_7_dry_full_noint,
  #   "stan/model_ind.stan",
  #   data = dry_full_noint,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 2000,
  #   iter_sampling = 1000,
  #   draws = FALSE,
  #   diagnostics = FALSE,
  #   summary = FALSE,
  #   adapt_delta = 0.95,
  #   max_treedepth = 15,
  #   seed = 123),
  # tar_stan_mcmc(
  #   fit_8_wet_full_noint,
  #   "stan/model_ind.stan",
  #   data = wet_full_noint,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 2000,
  #   iter_sampling = 1000,
  #   draws = FALSE,
  #   diagnostics = FALSE,
  #   summary = FALSE,
  #   adapt_delta = 0.95,
  #   max_treedepth = 15,
  #   seed = 123),
  # tar_stan_mcmc(
  #   fit_9_dry_wd_noint,
  #   "stan/model_ind.stan",
  #   data = dry_wd_noint,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 2000,
  #   iter_sampling = 1000,
  #   draws = FALSE,
  #   diagnostics = FALSE,
  #   summary = FALSE,
  #   adapt_delta = 0.95,
  #   max_treedepth = 15,
  #   seed = 123),
  # tar_stan_mcmc(
  #   fit_10_wet_wd_noint,
  #   "stan/model_ind.stan",
  #   data = wet_wd_noint,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 2000,
  #   iter_sampling = 1000,
  #   draws = FALSE,
  #   diagnostics = FALSE,
  #   summary = FALSE,
  #   adapt_delta = 0.95,
  #   max_treedepth = 15,
  #   seed = 123),
  # tar_stan_mcmc(
  #   fit_11_dry_pca_noint,
  #   "stan/model_ind.stan",
  #   data = dry_pca_noint,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 2000,
  #   iter_sampling = 1000,
  #   draws = FALSE,
  #   diagnostics = FALSE,
  #   summary = FALSE,
  #   adapt_delta = 0.95,
  #   max_treedepth = 15,
  #   seed = 123),
  # tar_stan_mcmc(
  #   fit_12_wet_pca_noint,
  #   "stan/model_ind.stan",
  #   data = wet_pca_noint,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 2000,
  #   iter_sampling = 1000,
  #   draws = FALSE,
  #   diagnostics = FALSE,
  #   summary = FALSE,
  #   adapt_delta = 0.95,
  #   max_treedepth = 15,
  #   seed = 123)#,

  # tar_target(
  #   dry_full_coef_data,
  #   create_stan_tab(fit_1_dry_full_int_draws_model_ind)
  # ),
  # tar_target(
  #   wet_full_coef_data,
  #   create_stan_tab(fit_2_wet_full_int_draws_model_ind)
  # ),
  # tar_target(
  #   dry_wd_coef_data,
  #   create_stan_tab(fit_3_dry_wd_int_draws_model_ind)
  # ),
  # tar_target(
  #   wet_wd_coef_data,
  #   create_stan_tab(fit_4_wet_wd_int_draws_model_ind)
  # ),
  # tar_target(
  #   dry_pca_coef_data,
  #   create_stan_tab(fit_5_dry_pca_int_draws_model_ind)
  # ),
  # tar_target(
  #   wet_pca_coef_data,
  #   create_stan_tab(fit_6_wet_pca_int_draws_model_ind)
  # ),

  # tar_render(
  #   report,
  #   "report.Rmd"
  # )
)
