library(targets)
library(tarchetypes)
library(tidyverse)
library(stantargets)
library(cmdstanr)
library(furrr)
library(languageserver)
# library(quarto)

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

list(
  # data cleaning ----------------------------------
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

  # stan
  tar_target(
    scale_cc,
    calc_scale_cc(seedling_csv)
  ),

  tar_target(
    test,
    generate_stan_data(seedling_csv, trait_csv, scale_cc, model = "het_season", abund = "ba")
  ),

  tar_map(
    list(model = c("phy_season", "het_season", "phy_rain", "het_rain")),
    tar_target(ab,
      generate_stan_data(seedling_csv, trait_csv, scale_cc,
       model = model, abund = "ab")),
    tar_target(ba,
      generate_stan_data(seedling_csv, trait_csv, scale_cc,
       model = model, abund = "ba")),
    tar_target(both,
      generate_stan_data(seedling_csv, trait_csv, scale_cc,
       model = model, abund = "both"))
  ),


  NULL
 )

