library(targets)
library(tarchetypes)
library(tidyverse)
library(stantargets)
library(cmdstanr)
library(furrr)
library(languageserver)
library(quarto)

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

  # scaling across seasons
  tar_target(
    dry_each_oneint,
    gen_stan_dat(data_list, season = "dry",
      inter = TRUE, trait_set = "each", one_inter = TRUE),
  ),
  tar_target(
    dry_each_int,
    gen_stan_dat(data_list, season = "dry",
      inter = TRUE, trait_set = "each"),
  ),
  tar_target(
    wet_each_int,
    gen_stan_dat(data_list, season = "rainy",
      inter = TRUE, trait_set = "each"),
  ),
  tar_target(
    dry_each_noint,
    gen_stan_dat(data_list, season = "dry",
      inter = FALSE, trait_set = "each"),
  ),
  tar_target(
    wet_each_noint,
    gen_stan_dat(data_list, season = "rainy",
      inter = FALSE, trait_set = "each"),
  ),
  tar_target(
    dry_pca_int,
    gen_stan_dat(data_list, season = "dry",
      inter = TRUE, trait_set = "pca"),
  ),
  tar_target(
    wet_pca_int,
    gen_stan_dat(data_list, season = "rainy",
      inter = TRUE, trait_set = "pca"),
  ),
  tar_target(
    dry_pca_noint,
    gen_stan_dat(data_list, season = "dry",
      inter = FALSE, trait_set = "pca"),
  ),
  tar_target(
    wet_pca_noint,
    gen_stan_dat(data_list, season = "rainy",
      inter = FALSE, trait_set = "pca"),
  ),

  # scaling within seasons
  tar_target(
    dry_each_int_s,
    gen_stan_dat(data_list, season = "dry",
      inter = TRUE, trait_set = "each",
      scaling_within_seasons = TRUE),
  ),
  tar_target(
    wet_each_int_s,
    gen_stan_dat(data_list, season = "rainy",
      inter = TRUE, trait_set = "each",
      scaling_within_seasons = TRUE),
  ),
  tar_stan_mcmc(
    fit_0_dry_each_oneint,
    "stan/model_ind.stan",
    data = dry_each_oneint,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 1000,
    iter_sampling = 1000,
    draws = TRUE,
    diagnostics = TRUE,
    summary = TRUE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_1_dry_each_int,
    "stan/model_ind.stan",
    data = dry_each_int,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 1000,
    iter_sampling = 1000,
    draws = TRUE,
    diagnostics = TRUE,
    summary = TRUE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_2_wet_each_int,
    "stan/model_ind.stan",
    data = wet_each_int,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 1000,
    iter_sampling = 1000,
    draws = TRUE,
    diagnostics = TRUE,
    summary = TRUE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_3_dry_each_noint,
    "stan/model_ind.stan",
    data = dry_each_noint,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 1000,
    iter_sampling = 1000,
    draws = TRUE,
    diagnostics = TRUE,
    summary = TRUE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_4_wet_each_noint,
    "stan/model_ind.stan",
    data = wet_each_noint,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 1000,
    iter_sampling = 1000,
    draws = TRUE,
    diagnostics = TRUE,
    summary = TRUE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_5_dry_pca_int,
    "stan/model_ind.stan",
    data = dry_pca_int,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 1000,
    iter_sampling = 1000,
    draws = TRUE,
    diagnostics = TRUE,
    summary = TRUE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_6_wet_pca_int,
    "stan/model_ind.stan",
    data = wet_pca_int,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 1000,
    iter_sampling = 1000,
    draws = TRUE,
    diagnostics = TRUE,
    summary = TRUE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_7_dry_pca_noint,
    "stan/model_ind.stan",
    data = dry_pca_noint,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 1000,
    iter_sampling = 1000,
    draws = TRUE,
    diagnostics = TRUE,
    summary = TRUE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_8_wet_pca_noint,
    "stan/model_ind.stan",
    data = wet_pca_noint,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 1000,
    iter_sampling = 1000,
    draws = TRUE,
    diagnostics = TRUE,
    summary = TRUE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_9_dry_each_int_s,
    "stan/model_ind.stan",
    data = dry_each_int_s,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 1000,
    iter_sampling = 1000,
    draws = TRUE,
    diagnostics = TRUE,
    summary = TRUE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  tar_stan_mcmc(
    fit_10_wet_each_int_s,
    "stan/model_ind.stan",
    data = wet_each_int_s,
    refresh = 0,
    chains = 4,
    parallel_chains = getOption("mc.cores", 4),
    iter_warmup = 1000,
    iter_sampling = 1000,
    draws = TRUE,
    diagnostics = TRUE,
    summary = TRUE,
    adapt_delta = 0.95,
    max_treedepth = 15,
    seed = 123),
  # # Better not to use `mclapply`. It requries too much RAM
  tar_target(
    dry_loo,
    lapply(
      list(
          fit_0_dry_each_oneint_mcmc_model_ind = fit_0_dry_each_oneint_mcmc_model_ind,
          fit_1_dry_each_int_mcmc_model_ind = fit_1_dry_each_int_mcmc_model_ind,
          fit_3_dry_each_noint_mcmc_model_ind = fit_3_dry_each_noint_mcmc_model_ind,
          fit_5_dry_pca_int_mcmc_model_ind = fit_5_dry_pca_int_mcmc_model_ind,
          fit_7_dry_pca_noint_mcmc_model_ind = fit_7_dry_pca_noint_mcmc_model_ind
        ),
    \(x)x$loo(cores = parallel::detectCores())
    )
  ),
  tar_target(
    wet_loo,
    lapply(
      list(
          fit_2_wet_each_int_mcmc_model_ind = fit_2_wet_each_int_mcmc_model_ind,
          fit_4_wet_each_noint_mcmc_model_ind = fit_4_wet_each_noint_mcmc_model_ind,
          fit_6_wet_pca_int_mcmc_model_ind = fit_6_wet_pca_int_mcmc_model_ind,
          fit_8_wet_pca_noint_mcmc_model_ind = fit_8_wet_pca_noint_mcmc_model_ind
        ),
    \(x)x$loo(cores = parallel::detectCores())
    )
  ),
  # tar_render(
  #   data_check_html,
  #   "docs/data_check.Rmd",
  #   output_format = "html_document"
  #   #knit_root_dir = here::here()
  # ),
  tar_render(
    bayes_check_html,
    "docs/bayes_check.Rmd",
    output_format = "html_document"
  ),
  tar_render(
    bayes_check_all_html,
    "docs/bayes_check_all.Rmd",
    output_format = "html_document"
  ),
  # tar_quarto(
  #   si_pdf,
  #   "ms/SI.qmd"
  #   # output_format = "pdf"
  # ),
  # tar_target(
  #   si_pdf, {
  #     quarto_render(
  #       "ms/SI.qmd",
  #       output_format = "pdf"
  #     )
  #     paste("ms/SI.pdf")
  #   },
  #   format = "file"
  # ),
  # tar_render(
  #   vis_idea_html,
  #   "docs/vis_idea.Rmd",
  #   output_format = "html_document"
  #   #knit_root_dir = here::here()
  # ),

  tar_target(
    fit0_tab,
    create_stan_tab(fit_0_dry_each_oneint_draws_model_ind)
  ),
  tar_target(
    fit1_tab,
    create_stan_tab(fit_1_dry_each_int_draws_model_ind)
  ),
  tar_target(
    fit3_tab,
    create_stan_tab(fit_3_dry_each_noint_draws_model_ind)
  ),
  tar_target(
    fit5_tab,
    create_stan_tab(fit_5_dry_pca_int_draws_model_ind)
  ),
  tar_target(
    fit7_tab,
    create_stan_tab(fit_7_dry_pca_noint_draws_model_ind)
  ),
  tar_target(
    fit2_tab,
    create_stan_tab(fit_2_wet_each_int_draws_model_ind)
  ),
  tar_target(
    fit4_tab,
    create_stan_tab(fit_4_wet_each_noint_draws_model_ind)
  ),
  tar_target(
    fit6_tab,
    create_stan_tab(fit_6_wet_pca_int_draws_model_ind)
  ),
  tar_target(
    fit8_tab,
    create_stan_tab(fit_8_wet_pca_noint_draws_model_ind)
  ),
  tar_target(
    fit9_tab,
    create_stan_tab(fit_9_dry_each_int_s_draws_model_ind)
  ),
  tar_target(
    fit10_tab,
    create_stan_tab(fit_10_wet_each_int_s_draws_model_ind)
  ),
  tar_target(
    fit1_gamma,
    create_gamma_tab(fit1_tab, dry_each_int)
  ),
  tar_target(
    fit2_gamma,
    create_gamma_tab(fit2_tab, wet_each_int)
  ),
  tar_target(
    fit3_gamma,
    create_gamma_tab(fit3_tab, dry_each_noint)
  ),
  tar_target(
    fit4_gamma,
    create_gamma_tab(fit4_tab, wet_each_noint)
  ),
  tar_target(
    fit9_gamma,
    create_gamma_tab(fit9_tab, dry_each_int_s)
  ),
  tar_target(
    fit10_gamma,
    create_gamma_tab(fit10_tab, wet_each_int_s)
  ),
  tar_target(
    dry_gamma_csv, {
      write_csv(fit9_gamma, "data/dry_gamma.csv")
      paste("data/dry_gamma.csv")
    },
    format = "file"
  ),
  tar_target(
    wet_gamma_csv, {
      write_csv(fit10_gamma, "data/wet_gamma.csv")
      paste("data/wet_gamma.csv")
    },
    format = "file"
  ),

  tar_target(
    test_ridge_plot, {
    p <- coef_ridge(fit1_tab, dry_each_int, fit_1_dry_each_int_draws_model_ind)
      ggsave(
        "figs/test_ridge.png",
        p,
        dpi = 300,
        width = 3.5,
        height = 3.5)
      ggsave(
        "figs/test_ridge.pdf",
        p,
        width = 3.5,
        height = 3.5)
      paste0("figs/test_ridge", c(".png", ".pdf"))
    },
    format = "file"
  ),

  tar_target(
    coef_trait_int_plot, {
      p <- coef_pointrange(fit1_gamma, fit2_gamma)
      ggsave(
        "figs/coef_trait_int.png",
        p,
        dpi = 300,
        width = 6,
        height = 3)
      ggsave(
        "figs/coef_trait_int.pdf",
        p,
        width = 6,
        height = 3)
      paste0("figs/coef_trait_int", c(".png", ".pdf"))
    },
    format = "file"
  ),
  tar_target(
    coef_trait_int_s_plot, {
      p <- coef_pointrange(fit9_gamma, fit10_gamma)
      ggsave(
        "figs/coef_trait_int_s.png",
        p,
        dpi = 300,
        width = 6,
        height = 3)
      ggsave(
        "figs/coef_trait_int_s.pdf",
        p,
        width = 6,
        height = 3)
      paste0("figs/coef_trait_int_s", c(".png", ".pdf"))
    },
    format = "file"
  ),
  tar_target(
    coef_trait_noint_plot, {
      p <- coef_pointrange(fit3_gamma, fit4_gamma)
      ggsave(
        "figs/coef_trait_noint.png",
        p,
        dpi = 300,
        width = 6,
        height = 3)
      ggsave(
        "figs/coef_trait_noint.pdf",
        p,
        width = 6,
        height = 3)
      paste0("figs/coef_trait_noint", c(".png", ".pdf"))
    },
    format = "file"
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

)
