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

values <- expand_grid(
  season = c("dry", "wet"),
  het = c("phy", "het"),
  rain = c("norain", "intrain", "rain"),
  sp_pred = c("nlog", "n", "ab", "ba", "ab1ba", "ab2ba")
  )
values_inter <- expand_grid(
  season = c("dry", "wet"),
  het = "het",
  rain = "intrain2",
  sp_pred = "nlog")

data_names <- values |>
  mutate(data_names = str_c(season, het, rain, sp_pred, sep = "_")) |>
  pull(data_names)

mcmc_names <- values |>
  mutate(mcmc_names = str_c("fit_mcmc_logistic_stan_data_", data_names)) |>
  pull(mcmc_names)

mcmc_names2 <- str_replace_all(mcmc_names, "_stan", "_simple_stan")

mcmc_names3 <- c(
  "fit_mcmc_logistic_simple_stan_data_dry_het_intrain2_nlog",
  "fit_mcmc_logistic_simple_stan_data_wet_het_intrain2_nlog",
  "fit2_mcmc_logistic_simple_stan_data_wet_het_norain_nlog",
  "fit2_mcmc_logistic_simple_stan_data_wet_phy_norain_nlog",
  "fit2_mcmc_logistic_simple_stan_data_wet_phy_rain_nlog")

loo_map <- tar_map(
    values = list(mcmc = rlang::syms(c(mcmc_names, mcmc_names2, mcmc_names3))),
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
  tar_target(
    scale_wet,
    calc_scale_cc(seedling_csv, wet = TRUE),
  ),
  tar_target(
    scale_dry,
    calc_scale_cc(seedling_csv, wet = FALSE),
  ),
  tar_target(
    wet_cc_data,
    generate_cc_data(seedling_csv, wet = TRUE),
  ),
  tar_target(
    dry_cc_data,
    generate_cc_data(seedling_csv, wet = FALSE),
  ),

  tar_target(
    cc_line_plot, {
      p <- cc_line(wet = wet_cc_data, dry = dry_cc_data)
      my_ggsave(
        "figs/cc_line",
        p,
        dpi = 300,
        width = 173,
        height = 173,
        units = "mm"
      )
    },
    format = "file"
  ),

  tar_map(
    values = values,
    tar_target(stan_data,
      generate_stan_data(
        seedling_csv, trait_csv,
        scale_cc = list(wet = scale_wet, dry = scale_dry),
        season, het, rain, sp_pred))
  ),
  tar_map(
    values = values_inter,
    tar_target(stan_data,
      generate_stan_data(
        seedling_csv, trait_csv,
        scale_cc = list(wet = scale_wet, dry = scale_dry),
        season, het, rain, sp_pred))
  ),
  # compile stan model so that targets can track
  # eventually, I need to run `_targets.R` twice.
  tar_target(
    logistic_stan,
    compile_model("stan/logistic.stan"),
    format = "file",
    deployment = "main"
  ),

  tar_map(
    values = list(stan_data = rlang::syms(str_c("stan_data_", data_names))),
    tar_stan_mcmc(
      fit,
      c("stan/logistic.stan", "stan/logistic_simple.stan"),
      data = stan_data,
      refresh = 0,
      chains = 4,
      parallel_chains = getOption("mc.cores", 1),
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.9,
      max_treedepth = 15,
      seed = 123,
      return_draws = FALSE,
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

  tar_map(
    values = list(stan_data = rlang::syms(str_c("stan_data_", c("dry", "wet"), "_het_intrain2_nlog"))),
    tar_stan_mcmc(
      fit,
      "stan/logistic_simple.stan",
      data = stan_data,
      refresh = 0,
      chains = 4,
      parallel_chains = getOption("mc.cores", 4),
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.9,
      max_treedepth = 15,
      seed = 123,
      return_draws = FALSE,
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

  tar_map(
    values = list(stan_data = rlang::syms(c(str_c("stan_data_wet_", c("het", "phy"), "_norain_nlog"), "stan_data_wet_phy_rain_nlog"))),
    tar_stan_mcmc(
      fit2,
      "stan/logistic_simple.stan",
      data = stan_data,
      refresh = 0,
      chains = 4,
      parallel_chains = getOption("mc.cores", 4),
      iter_warmup = 1000,
      iter_sampling = 2000,
      adapt_delta = 0.95,
      max_treedepth = 15,
      seed = 123,
      return_draws = FALSE,
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

  # tar_stan_mcmc(
  #   check_ess,
  #   "stan/logistic.stan",
  #   data = stan_data_dry_het_intrain_ab,
  #   refresh = 0,
  #   chains = 4,
  #   parallel_chains = getOption("mc.cores", 4),
  #   iter_warmup = 2000,
  #   iter_sampling = 2000,
  #   adapt_delta = 0.9,
  #   max_treedepth = 15,
  #   seed = 123,
  #   return_draws = TRUE,
  #   return_diagnostics = TRUE,
  #   return_summary = TRUE,
  #   summaries = list(
  #     mean = ~mean(.x),
  #     sd = ~sd(.x),
  #     mad = ~mad(.x),
  #     ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
  #     posterior::default_convergence_measures()
  #   )
  # ),

  loo_map,
  tar_combine(
    loo_list,
    loo_map,
    command = list(!!!.x)
  ),

  # tar_quarto(
  #   bayes_check_html,
  #   "docs/bayes_check.qmd"
  # ),

  tar_target(
    loo_tbl,
    generate_loo_tbl(loo_list)
  ),
  tar_target(
    dry_trait,
    load_mcmc_summary(loo_tbl, season = "dry", trait = "n")
  ),
  tar_target(
    wet_trait,
    load_mcmc_summary(loo_tbl, season = "wet", trait = "n")
  ),

  tar_target(
    dry_trait_suv_contour_plot, {
      p <- dry_trait_suv_contour(dry_trait, alpha = 0.05)
      my_ggsave(
        "figs/dry_trait_suv_contour",
        p,
        dpi = 300,
        width = 173,
        height = 260,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    wet_trait_suv_contour_plot, {
      p <- wet_trait_suv_contour(wet_trait, alpha = 0.05)
      my_ggsave(
        "figs/wet_trait_suv_contour",
        p,
        dpi = 300,
        width = 173,
        height = 180,
        units = "mm"
      )
    },
    format = "file"
  ),


  tar_target(
    beta_wet_rain_n,
    generate_beta_list(
      wet_trait$draws,
      wet_trait$data,
      x_lab = "N",
      y_lab = "Rain~effect",
      ind_pred = 7,
      sp_pred = 8)
  ),
  tar_target(
    beta_wet_rain_tlp,
    generate_beta_list(
      wet_trait$draws,
      wet_trait$data,
      x_lab = expression(pi[tlp]),
      y_lab = "Rain~effect",
      ind_pred = 7,
      sp_pred = 9)
  ),

  tar_target(
    beta_wet_consrain_sla,
    generate_beta_list(
      wet_trait$draws,
      wet_trait$data,
      x_lab = "SLA",
      y_lab  = "ConS%*%Rainfall~effect",
      ind_pred = 8,
      sp_pred = 5)
  ),
  tar_target(
    beta_wet_consrain_n,
    generate_beta_list(
      wet_trait$draws,
      wet_trait$data,
      x_lab = "N",
      y_lab  = "ConS%*%Rainfall~effect",
      ind_pred = 8,
      sp_pred = 8)
  ),
  tar_target(
    beta_wet_consrain_tlp,
    generate_beta_list(
      wet_trait$draws,
      wet_trait$data,
      x_lab = expression(pi[tlp]),
      y_lab  = "ConS%*%Rainfall~effect",
      ind_pred = 8,
      sp_pred = 9)
  ),

  tar_target(
    beta_dry_rain_ldmc,
    generate_beta_list(
      dry_trait$draws,
      dry_trait$data,
      x_lab = "LDMC",
      y_lab = "Rain~effect",
      ind_pred = 7,
      sp_pred = 2)
  ),
  tar_target(
    beta_dry_consrain_ldmc,
    generate_beta_list(
      dry_trait$draws,
      dry_trait$data,
      x_lab = "LDMC",
      y_lab  = "ConS%*%Rainfall~effect",
      ind_pred = 8,
      sp_pred = 2)
  ),
  tar_target(
    beta_dry_rain_sdmc,
    generate_beta_list(
      dry_trait$draws,
      dry_trait$data,
      x_lab = "SDMC",
      y_lab = "Rain~effect",
      ind_pred = 7,
      sp_pred = 3)
  ),
  tar_target(
    beta_dry_cons_sdmc,
    generate_beta_list(
      dry_trait$draws,
      dry_trait$data,
      x_lab = "SDMC",
      y_lab  = "ConS~effect",
      ind_pred = 3,
      sp_pred = 3)
  ),
  tar_target(
    beta_dry_rain_lt,
    generate_beta_list(
      dry_trait$draws,
      dry_trait$data,
      x_lab = "LT",
      y_lab = "Rain~effect",
      ind_pred = 7,
      sp_pred = 6)
  ),
  tar_target(
    beta_dry_consrain_lt,
    generate_beta_list(
      dry_trait$draws,
      dry_trait$data,
      x_lab = "LT",
      y_lab  = "ConS%*%Rainfall~effect",
      ind_pred = 8,
      sp_pred = 6)
  ),
  tar_target(
    beta_dry_rain_c13,
    generate_beta_list(
      dry_trait$draws,
      dry_trait$data,
      x_lab = expression(delta*C[13]),
      y_lab = "Rain~effect",
      ind_pred = 7,
      sp_pred = 7)
  ),

  tar_target(
    beta_par_wet_traits, {
      p <- beta_plot(beta_wet_rain_n) +
        beta_plot(beta_wet_rain_tlp) +
        plot_spacer() +
        beta_plot(beta_wet_consrain_sla) +
        beta_plot(beta_wet_consrain_n) +
        beta_plot(beta_wet_consrain_tlp) +
        plot_layout(ncol = 3, nrow = 2) +
        plot_annotation(tag_levels = "a") &
        theme(
          text = element_text(size = 8),
          plot.tag = element_text(face = "bold"))
      my_ggsave(
        "figs/beta_par_wet_traits",
        p,
        dpi = 300,
        width = 110,
        height = 75,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    beta_raw_wet_traits, {
      p <- beta_plot(beta_wet_rain_n, partial = FALSE) +
        beta_plot(beta_wet_rain_tlp, partial = FALSE) +
        plot_spacer() +
        beta_plot(beta_wet_consrain_sla, partial = FALSE) +
        beta_plot(beta_wet_consrain_n, partial = FALSE) +
        beta_plot(beta_wet_consrain_tlp, partial = FALSE) +
        plot_layout(ncol = 3, nrow = 2) +
        plot_annotation(tag_levels = "a") &
        theme(
          text = element_text(size = 8),
          plot.tag = element_text(face = "bold"))
      my_ggsave(
        "figs/beta_raw_wet_traits",
        p,
        dpi = 300,
        width = 110,
        height = 75,
        units = "mm"
      )
    },
    format = "file"
  ),

  tar_target(
    beta_par_dry_traits, {
      p <- beta_plot(beta_dry_cons_sdmc) +
        beta_plot(beta_dry_rain_sdmc) +
        beta_plot(beta_dry_rain_ldmc) +
        beta_plot(beta_dry_rain_lt) +
        beta_plot(beta_dry_rain_c13) +
        beta_plot(beta_dry_consrain_ldmc) +
        beta_plot(beta_dry_consrain_lt) +
        plot_spacer() +
        plot_layout(ncol = 4, nrow = 2) +
        plot_annotation(tag_levels = "a") &
        theme(
          text = element_text(size = 8),
          plot.tag = element_text(face = "bold"))
      my_ggsave(
        "figs/beta_par_dry_traits",
        p,
        dpi = 300,
        width = 173,
        height = 85,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    beta_raw_dry_traits, {
      p <- beta_plot(beta_dry_cons_sdmc, partial = FALSE) +
        beta_plot(beta_dry_rain_sdmc, partial = FALSE) +
        beta_plot(beta_dry_rain_ldmc, partial = FALSE) +
        beta_plot(beta_dry_rain_lt, partial = FALSE) +
        beta_plot(beta_dry_rain_c13, partial = FALSE) +
        beta_plot(beta_dry_consrain_ldmc, partial = FALSE) +
        beta_plot(beta_dry_consrain_lt, partial = FALSE) +
        plot_spacer() +
        plot_layout(ncol = 4, nrow = 2) +
        plot_annotation(tag_levels = "a") &
        theme(
          text = element_text(size = 8),
          plot.tag = element_text(face = "bold"))
      my_ggsave(
        "figs/beta_raw_dry_traits",
        p,
        dpi = 300,
        width = 173,
        height = 85,
        units = "mm"
      )
    },
    format = "file"
  ),

# best models
tar_map(
  values = list(
    x = rlang::syms(c(
      "fit_summary_logistic_simple_stan_data_dry_het_intrain_ab",
      "fit_summary_logistic_simple_stan_data_dry_het_intrain2_nlog",
      "fit_summary_logistic_simple_stan_data_wet_het_intrain_ab",
      "fit_summary_logistic_simple_stan_data_wet_het_intrain2_nlog")),
    stan_data = rlang::syms(c(
      "stan_data_dry_het_intrain_ab",
      "stan_data_dry_het_intrain2_nlog",
      "stan_data_wet_het_intrain_ab",
      "stan_data_wet_het_intrain2_nlog")),
    path =
      str_c(
      "data/",
       c("dry_abund", "dry_traits", "wet_abund", "wet_traits"),
      "_gamma.csv")),
  tar_target(
    gamma_out_csv, {
      create_gamma_tab(x, stan_data)  |>
        my_write_csv(path)
    },
    format = "file"
  )
),


  # tar_target(
  #   wet_gamma_csv, {
  #     write_csv(fit10_gamma, "data/wet_gamma.csv")
  #     paste("data/wet_gamma.csv")
  #   },
  #   format = "file"
  # ),

  # tar_quarto(
  #   bayes_check_html,
  #   "docs/bayes_check.qmd",
  # ),

  NULL
 )


list(data_, main_)
