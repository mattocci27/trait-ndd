library(targets)
library(tarchetypes)
library(tidyverse)
library(stantargets)
library(cmdstanr)
library(furrr)
library(languageserver)

source("R/data_clean.R")
source("R/data_check.R")
source("R/stan.R")
source("R/figs.R")

plan(multicore)
options(clustermq.scheduler = "multicore")

tar_option_set(packages = c(
  "tidyverse",
  "GGally",
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
    seedling_clean_csv,
    seedling_clean(seedling_csv),
    format = "file"
  ),
  tar_target(
    trait_clean_csv,
    trait_clean(trait_csv),
    format = "file"
  ),

  # data check ------------------------------------------

  tar_target(
    demo_cov, {
    read_csv(seedling_clean_csv) |>
      dplyr::select(cons, cona, hets, heta) |>
      ggpairs()
    }
  ),

  tar_target(
    trait_log, {
    read_csv(trait_clean_csv) |>
      mutate(log_la = log(la)) |>
      mutate(log_sla = log(sla)) |>
      mutate(log_lt = log(lt)) |>
      dplyr::select(-la, -sla, -lt)
    }
  ),

  tar_target(
    trait_pair_plot, {
      p <- trait_pair(trait_log, full = TRUE)
      ggsave(
        "figs/trait_pair.png",
        p,
        dpi = 300,
        width = 7.25,
        height = 7.25
      )
      ggsave(
        "figs/trait_pair.pdf",
        p,
        device = cairo_pdf,
        width = 7.25,
        height = 7.25
      )
      paste0("figs/trait_pair", c(".png", ".pdf"))
    },
    format = "file"
  ),

  # stan -------------------------------------------------




  tar_target(
    data_list,
    gen_seedling(seedling_csv, trait_csv, habitat_csv, n_ab = 50)
  ),


  tar_target(
    dry_cn_int,
    gen_stan_dat(data_list, season = "dry", habitat = "all",
      inter = TRUE, trait_set = "cn"),
  ),
  tar_target(
    wet_cn_int,
    gen_stan_dat(data_list, season = "rainy", habitat = "all",
      inter = TRUE, trait_set = "cn"),
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
    dry_cn_noint,
    gen_stan_dat(data_list, season = "dry", habitat = "all",
      inter = FALSE, trait_set = "cn"),
  ),
  tar_target(
    wet_cn_noint,
    gen_stan_dat(data_list, season = "rainy", habitat = "all",
      inter = FALSE, trait_set = "cn"),
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

  # tar_stan_mcmc(
  #   fit_1_dry_cn_int,
  #   "stan/model_ind.stan",
  #   data = dry_cn_int,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 1000,
  #   iter_sampling = 1000,
  #   draws = FALSE,
  #   diagnostics = FALSE,
  #   summary = FALSE,
  #   adapt_delta = 0.9,
  #   max_treedepth = 15,
  #   seed = 123),
  # tar_stan_mcmc(
  #   fit_2_wet_cn_int,
  #   "stan/model_ind.stan",
  #   data = wet_cn_int,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 1000,
  #   iter_sampling = 1000,
  #   draws = FALSE,
  #   diagnostics = FALSE,
  #   summary = FALSE,
  #   adapt_delta = 0.9,
  #   max_treedepth = 15,
  #   seed = 123),
  # tar_stan_mcmc(
  #   fit_3_dry_wd_int,
  #   "stan/model_ind.stan",
  #   data = dry_wd_int,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 1000,
  #   iter_sampling = 1000,
  #   draws = FALSE,
  #   diagnostics = TRUE,
  #   summary = TRUE,
  #   adapt_delta = 0.9,
  #   max_treedepth = 15,
  #   seed = 123),
  # tar_stan_mcmc(
  #   fit_4_wet_wd_int,
  #   "stan/model_ind.stan",
  #   data = wet_wd_int,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 1000,
  #   iter_sampling = 1000,
  #   draws = FALSE,
  #   diagnostics = TRUE,
  #   summary = TRUE,
  #   adapt_delta = 0.9,
  #   max_treedepth = 15,
  #   seed = 123),

  tar_stan_mcmc(
    fit_5_dry_cn_noint,
    "stan/model_ind.stan",
    data = dry_cn_noint,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 10,
    iter_sampling = 10,
    draws = FALSE,
    diagnostics = FALSE,
    summary = FALSE,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_6_wet_cn_noint,
    "stan/model_ind.stan",
    data = wet_cn_noint,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 10,
    iter_sampling = 10,
    draws = FALSE,
    diagnostics = FALSE,
    summary = FALSE,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_7_dry_wd_noint,
    "stan/model_ind.stan",
    data = dry_wd_noint,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 10,
    iter_sampling = 10,
    draws = FALSE,
    diagnostics = FALSE,
    summary = FALSE,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_8_wet_wd_noint,
    "stan/model_ind.stan",
    data = wet_wd_noint,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 10,
    iter_sampling = 10,
    draws = FALSE,
    diagnostics = FALSE,
    summary = FALSE,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123),

  tar_stan_mcmc(
    fit_9_dry_pca_noint,
    "stan/model_ind.stan",
    data = dry_pca_noint,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 10,
    iter_sampling = 10,
    draws = FALSE,
    diagnostics = FALSE,
    summary = FALSE,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_10_wet_pca_noint,
    "stan/model_ind.stan",
    data = wet_pca_noint,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 10,
    iter_sampling = 10,
    draws = FALSE,
    diagnostics = FALSE,
    summary = FALSE,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_11_dry_pca_noint,
    "stan/model_ind.stan",
    data = dry_pca_noint,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 10,
    iter_sampling = 10,
    draws = FALSE,
    diagnostics = FALSE,
    summary = FALSE,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_12_wet_pca_noint,
    "stan/model_ind.stan",
    data = wet_pca_noint,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 10,
    iter_sampling = 10,
    draws = FALSE,
    diagnostics = FALSE,
    summary = FALSE,
    adapt_delta = 0.9,
    max_treedepth = 15,
    seed = 123),

  tar_render(
    data_check_html,
    "docs/data_check.Rmd",
    output_format = "html_document"
    #knit_root_dir = here::here()
  )


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
