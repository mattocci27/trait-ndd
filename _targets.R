library(targets)
library(tarchetypes)
library(tidyverse)
library(stantargets)
library(cmdstanr)
library(furrr)

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
  "loo",
  "jsonlite",
  "httpgd",
  "multcompView",
  "RColorBrewer",
  "ggridges"
))

tar_option_set(
  garbage_collection = TRUE,
  memory = "transient"
)

# check if it's inside a container
if (file.exists("/.dockerenv") | file.exists("/.singularity.d/startscript")) {
  Sys.setenv(CMDSTAN = "/opt/cmdstan/cmdstan-2.29.2")
  set_cmdstan_path("/opt/cmdstan/cmdstan-2.29.2")
}

cmdstan_version()

data_names <- expand_grid(a1 = c("phy", "het"), a2 = c("season", "rain"), a3 = c("ab", "ba", "both")) |>
  mutate(data = str_c(a1, a2, a3, sep = "_")) |>
  pull(data)


mcmc_names <- str_c("fit_mcmc_logistic_", data_names)

loo_map <- tar_map(
    values = list(mcmc = rlang::syms(mcmc_names)),
    tar_target(
      loo,
      my_loo(mcmc)
    )
)

# data cleaning ----------------------------------
data_ <- list(
  tar_target(
    data_rda,
    "data-raw/dataCNDD.Rdata",
    format = "file"
  ),
  tar_target(
    trait_csv,
    clean_data("data-raw/dataCNDD.Rdata", write_trait = TRUE),
    format = "file"
  ),
  tar_target(
    seedling_csv,
    clean_data("data-raw/dataCNDD.Rdata", write_trait = FALSE),
    format = "file"
  ),

  NULL
)

main_ <- list(
  # stan
  tar_target(
    scale_cc,
    calc_scale_cc(seedling_csv)
  ),

  tar_map(
    list(abund = c("ab", "ba", "both")),
    tar_target(phy_season,
      generate_stan_data(seedling_csv, trait_csv, scale_cc,
       model = "phy_season", abund = abund)),
    tar_target(phy_rain,
      generate_stan_data(seedling_csv, trait_csv, scale_cc,
       model = "phy_rain", abund = abund)),
    tar_target(het_season,
      generate_stan_data(seedling_csv, trait_csv, scale_cc,
       model = "het_season", abund = abund)),
    tar_target(het_rain,
      generate_stan_data(seedling_csv, trait_csv, scale_cc,
       model = "het_rain", abund = abund))
  ),

  tar_map(
    values = list(stan_data = rlang::syms(data_names)),
    tar_stan_mcmc(
      fit,
      "stan/logistic.stan",
      data = stan_data,
      refresh = 0,
      chains = 4,
      parallel_chains = getOption("mc.cores", 4),
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.9,
      max_treedepth = 15,
      seed = 123,
      return_draws = TRUE,
      return_diagnostics = TRUE,
      return_summary = TRUE,
      summaries = list(
        mean = ~mean(.x),
        sd = ~sd(.x),
        mad = ~mad(.x),
        ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
        posterior::default_convergence_measures()
      )
    )
  ),

  # loo_map,
  # tar_combine(
  #   loo_list,
  #   loo_map,
  #   command = list(!!!.x)
  # ),

  # tar_quarto(
  #   bayes_check_html,
  #   "docs/bayes_check.qmd"
  # ),

  NULL
 )


list(data_, main_)

