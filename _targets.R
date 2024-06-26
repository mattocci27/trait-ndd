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
  "ggridges",
  "bayesplot",
  "quarto",
  "lme4",
  "here",
  "kableExtra",
  "tictoc",
  "DT",
  "ggrepel",
  "factoextra",
  "FactoMineR",
  "ggpointdensity",
  "ggpmisc",
  "broom"
))

tar_option_set(
  garbage_collection = TRUE,
  memory = "transient"
)

# check if it's inside a container
if (file.exists("/.dockerenv") | file.exists("/.singularity.d/startscript")) {
  Sys.setenv(CMDSTAN = "/opt/cmdstan/cmdstan-2.33.1")
  set_cmdstan_path("/opt/cmdstan/cmdstan-2.33.1")
}

cmdstan_version()

values <- expand_grid(
  season = c("dry", "wet"),
  rain = c("norain", "rain", "intrain", "intrain2", "intrain3", "intrain4"),
  sp_pred = c("nlog", "ab", "ba", "pc12")
  )

# values <- expand_grid(
#   season = c("dry", "wet"),
#   rain = "norain",
#   sp_pred = "ab")

data_names <- values |>
  mutate(data_names = str_c(season, rain, sp_pred, sep = "_")) |>
  pull(data_names)

mcmc_names <- values |>
  mutate(mcmc_names = str_c("fit_mcmc_suv_ind_", data_names)) |>
  pull(mcmc_names)

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
    traits_csv,
    clean_data("data-raw/dataCNDD.Rdata", write_trait = TRUE),
    format = "file"
  ),
  tar_target(
    traits_df,
    read_csv(traits_csv) |> janitor::clean_names() |> dplyr::select(-cname)
  ),
  tar_target(
    seedling_df,
    read_csv(seedling_csv) |> janitor::clean_names()
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
    all_cc_data,
    generate_cc_all_data(seedling_csv),
  ),

  tar_target(
    cc_line_plot, {
      p <- cc_line(wet = wet_cc_data, dry = dry_cc_data, all = all_cc_data)
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
  tar_target(
    pca_panel_plot, {
      p <- pca_panel(traits_df)
      my_ggsave(
        "figs/pca_panel",
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
    ggpairs_plot, {
      p <- my_ggpairs(traits_csv)
      my_ggsave(
        "figs/pairs",
        p,
        dpi = 300,
        width = 180,
        height = 180,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    phy_het_plot, {
      p <- phy_het_points(seedling_df)
      my_ggsave(
        "figs/phy_het",
        p,
        dpi = 300,
        width = 173,
        height = 173,
        units = "mm"
      )
    },
    format = "file"
  ),

  tar_target(
    test_stan_data,
    generate_stan_data(
      seedling_df, traits_df,
      scale_cc = list(wet = scale_wet, dry = scale_dry),
      season = "dry", rain = "intrain2", sp_pred = "pc12")
  ),
  tar_stan_mcmc(
    test_fit,
    c("stan/suv_ind.stan", "stan/suv.stan", "stan/suv_simple.stan"),
    data = test_stan_data,
    refresh = 0,
    chains = 1,
    parallel_chains = getOption("mc.cores", 3),
    iter_warmup = 1,
    iter_sampling = 1,
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
  ),
  # tar_target(
  #   test_loo_suv,
  #   my_loo(test_fit_mcmc_suv)
  # ),
  # tar_target(
  #   test_loo_suv_ind,
  #   my_loo(test_fit_mcmc_suv_ind)
  # ),
  # tar_target(
  #   test_loo_suv_simple,
  #   my_loo(test_fit_mcmc_suv_simple)
  # ),
  tar_map(
    values = tibble(phy = c("phy", "het")),
    tar_target(
      fit_glm,
      pre_glm(seedling_df, phy = phy)
    )
  ),
  tar_map(
    values = expand_grid(wet = c("wet", "dry"), phy = c("phy", "het")),
    tar_target(
      fit_glmer,
      pre_glmm(seedling_df, wet = wet, phy = phy)
    )
  ),
  tar_map(
    values = expand_grid(wet = c("wet", "dry"), phy = c("phy", "het")),
    tar_target(
      fit2_glmer,
      pre_glmm2(seedling_df, wet = wet, phy = phy)
    )
  ),
  tar_map(
    values = values |>
      filter((season == "dry" & sp_pred %in% c("nlog", "pc12")) | (season == "wet" & sp_pred == "nlog")),
    tar_target(stan_data,
      generate_stan_data(
        seedling_df, traits_df,
        scale_cc = list(wet = scale_wet, dry = scale_dry),
        season, rain, sp_pred)),
    tar_stan_mcmc(
      fit,
      "stan/suv_ind.stan",
      data = stan_data,
      refresh = 0,
      chains = 4,
      parallel_chains = getOption("mc.cores", 4),
      iter_warmup = 1000,
      iter_sampling = 2000,
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
    values = values |>
      filter(season == "wet" & sp_pred == "pc12"),
    tar_target(stan_data,
      generate_stan_data(
        seedling_df, traits_df,
        scale_cc = list(wet = scale_wet, dry = scale_dry),
        season, rain, sp_pred)),
    tar_stan_mcmc(
      fit,
      "stan/suv_ind.stan",
      data = stan_data,
      refresh = 0,
      chains = 4,
      parallel_chains = getOption("mc.cores", 4),
      iter_warmup = 1000,
      iter_sampling = 2000,
      adapt_delta = 0.99,
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
    values = values |> filter(sp_pred %in% c("ab", "ba")) |> filter(season == "dry"),
    tar_target(stan_data,
      generate_stan_data(
        seedling_df, traits_df,
        scale_cc = list(wet = scale_wet, dry = scale_dry),
        season, rain, sp_pred)),
    tar_stan_mcmc(
      fit,
      "stan/suv_ind.stan",
      data = stan_data,
      refresh = 0,
      chains = 4,
      parallel_chains = getOption("mc.cores", 4),
      iter_warmup = 1000,
      iter_sampling = 2000,
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
    values = values |> filter(sp_pred %in% c("ab", "ba")) |> filter(season == "wet"),
    tar_target(stan_data,
      generate_stan_data(
        seedling_df, traits_df,
        scale_cc = list(wet = scale_wet, dry = scale_dry),
        season, rain, sp_pred)),
    tar_stan_mcmc(
      fit,
      "stan/suv_ind.stan",
      data = stan_data,
      refresh = 0,
      chains = 4,
      parallel_chains = getOption("mc.cores", 4),
      iter_warmup = 1000,
      iter_sampling = 2000,
      adapt_delta = 0.99,
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
  loo_map,
  tar_combine(
    loo_list,
    loo_map,
    command = list(!!!.x)
  ),
  tar_target(
    loo_tbl,
    generate_loo_tbl(loo_list)
  ),
  NULL
)

fig_list <- list(
  # best models
  tar_target(
    dry_trait,
    list(
      data = stan_data_dry_intrain2_nlog,
      summary = fit_summary_suv_ind_dry_intrain2_nlog,
      draws =  posterior::as_draws_df(fit_mcmc_suv_ind_dry_intrain2_nlog)
    )
  ),
  tar_target(
    wet_trait,
    list(
      data = stan_data_wet_intrain2_nlog,
      summary = fit_summary_suv_ind_wet_intrain2_nlog,
      draws =  posterior::as_draws_df(fit_mcmc_suv_ind_wet_intrain2_nlog)
    )
  ),
  tar_target(
    dry_abund,
    list(
      data = stan_data_dry_intrain_ab,
      summary = fit_summary_suv_ind_dry_intrain_ab,
      draws =  posterior::as_draws_df(fit_mcmc_suv_ind_dry_intrain_ab)
    )
  ),
  tar_target(
    wet_abund,
    list(
      data = stan_data_wet_intrain_ab,
      summary = fit_summary_suv_ind_wet_intrain_ab,
      draws =  posterior::as_draws_df(fit_mcmc_suv_ind_wet_intrain_ab)
    )
  ),
  tar_target(
    dry_trait_suv_contour_plot, {
      p <- dry_trait_suv_contour(dry_trait, alpha = 0.05)
      my_ggsave(
        "figs/dry_trait_suv_contour",
        p,
        dpi = 600,
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
        dpi = 600,
        width = 173,
        height = 130,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    dry_data_coef,
    generate_coef_data(dry_trait$draws, dry_trait$data, season = "Dry")
  ),
  tar_target(
    wet_data_coef,
    generate_coef_data(wet_trait$draws, wet_trait$data, season = "Rainy")
  ),
  tar_target(
    coef_trait_plot, {
      p <- coef_pointrange(dry_data_coef, wet_data_coef, comb = FALSE)
      my_ggsave(
        "figs/coef_trait",
        p,
        dpi = 300,
        width = 173,
        height = 86,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    sdmc_res_data,
    generate_beta_list(dry_trait$draws, dry_trait$data, x_lab = "SDMC", y_lab = "ConS", ind_pred = 3, sp_pred = 3)
  ),
  tar_target(
    sdmc_partial_plot, {
      p <- beta_plot(sdmc_res_data)
      my_ggsave(
        "figs/sdmc_partial",
        p,
        dpi = 300,
        width = 82,
        height = 82,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    abund_res_data,
    generate_beta_list(dry_abund$draws, dry_abund$data, x_lab = "ln Abundance", y_lab = "ConS", ind_pred = 3, sp_pred = 2)
  ),
  tar_target(
    abund_partial_plot, {
      p <- beta_plot(abund_res_data)
      my_ggsave(
        "figs/abund_partial",
        p,
        dpi = 300,
        width = 82,
        height = 82,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    sdmc_abund_parital_plot, {
      p1 <- beta_plot(sdmc_res_data)
      p2 <- beta_plot(abund_res_data)
      p <- p1 + p2 +
        plot_annotation(tag_levels = "a") +
        theme(
          plot.tag = element_text(face = "bold")
        )
      my_ggsave(
        "figs/sdmc_abund_partial",
        p,
        dpi = 300,
        width = 110,
        height = 55,
        units = "mm"
      )
    }
  ),
  NULL
)

diagnostics_mapped <- tar_map(
    values = list(value = rlang::syms(str_c("fit_diagnostics_suv_ind_", data_names)), data_names = data_names),
    tar_target(
      diagnostics,
      value |> mutate(model = data_names)
    )
  )
tar_combined_diagnostics_data <- tar_combine(
  combined_diagnostics,
  diagnostics_mapped[["diagnostics"]],
  command = dplyr::bind_rows(!!!.x)
)

summary_mapped <- tar_map(
    values = list(value = rlang::syms(str_c("fit_summary_suv_ind_", data_names)), data_names = data_names),
    tar_target(
      summary,
      value |> mutate(model = data_names)
    )
  )
tar_combined_summary_data <- tar_combine(
  combined_summary,
  summary_mapped[["summary"]],
  command = dplyr::bind_rows(!!!.x)
)

# Define a list of model names
models <- c(
  "dry_intrain_pc12",
  "dry_intrain2_nlog",
  "wet_intrain_pc12",
  "wet_intrain2_nlog",
  "dry_intrain_ab",
  "wet_intrain_ab"
)

# Construct 'x' and 'stan_data' lists based on the model names
x <- rlang::syms(str_c("fit_summary_suv_ind_", models))
stan_data <- rlang::syms(str_c("stan_data_", models))

# Define a helper function to map model names to paths
get_path <- function(model_name) {
  str_c("data/", model_name, "_gamma.csv")
}

# Construct the path list based on the model names
path <- purrr::map_chr(models, get_path)

# Use 'tar_map' function with the derived variables
best_csv_mapped <- tar_map(
  values = list(x = x, stan_data = stan_data, path = path),
  tar_target(
    gamma_out_csv, {
      create_gamma_tab(x, stan_data)  |> my_write_csv(path)
    },
    format = "file"
  )
)

util_list <- list(
  diagnostics_mapped,
  tar_combined_diagnostics_data,
  summary_mapped,
  tar_combined_summary_data,
  tar_target(
    diagnostic_tables_csv,
    write_diagnostics_tables(combined_summary, combined_diagnostics, loo_tbl, "data/diagnostics_tables.csv"),
    format = "file"
  ),
  best_csv_mapped,
  NULL
)


list(data_, main_) |>
  append(fig_list) |>
  append(util_list)
