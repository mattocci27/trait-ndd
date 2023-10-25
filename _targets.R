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
  "bayesplot"
))

# tar_option_set(
#   garbage_collection = TRUE,
#   memory = "transient"
# )

# check if it's inside a container
if (file.exists("/.dockerenv") | file.exists("/.singularity.d/startscript")) {
  Sys.setenv(CMDSTAN = "/opt/cmdstan/cmdstan-2.29.2")
  set_cmdstan_path("/opt/cmdstan/cmdstan-2.29.2")
}

cmdstan_version()

values <- expand_grid(
  season = c("dry", "wet"),
  rain = c("norain", "rain", "intrain", "intrain2", "intrain3", "intrain4"),
  sp_pred = c("nlog", "ab", "ba", "ab1ba", "pc12")
  )

data_names <- values |>
  mutate(data_names = str_c(season, rain, sp_pred, sep = "_")) |>
  pull(data_names)

mcmc_names <- values |>
  mutate(mcmc_names = str_c("fit_mcmc_logistic_simple_stan_data_", data_names)) |>
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
        seedling_df, traits_df,
        scale_cc = list(wet = scale_wet, dry = scale_dry),
        season, rain, sp_pred))
  ),
  NULL)

hoge <- list(
  # tar_map(
  #   values = list(stan_data = rlang::syms(str_c("stan_data_", data_names))),
  #   tar_stan_mcmc(
  #     fit,
  #     "stan/logistic_simple.stan",
  #     data = stan_data,
  #     refresh = 0,
  #     chains = 4,
  #     parallel_chains = getOption("mc.cores", 4),
  #     iter_warmup = 1000,
  #     iter_sampling = 2000,
  #     adapt_delta = 0.95,
  #     max_treedepth = 15,
  #     seed = 123,
  #     return_draws = FALSE,
  #     return_diagnostics = TRUE,
  #     return_summary = TRUE,
  #     summaries = list(
  #       mean = ~mean(.x),
  #       sd = ~sd(.x),
  #       mad = ~mad(.x),
  #       ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
  #       posterior::default_convergence_measures()
  #     )
  #   )
  # ),

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
  # best models
  tar_target(
    dry_trait,
    load_mcmc_summary(loo_tbl, season = "dry", trait = "n")
  ),
  tar_target(
    wet_trait,
    load_mcmc_summary(loo_tbl, season = "wet", trait = "n")
  ),
  tar_target(
    dry_abund,
    load_mcmc_summary(loo_tbl, season = "dry", trait = "ab")
  ),
  tar_target(
    wet_abund,
    load_mcmc_summary(loo_tbl, season = "wet", trait = "ab")
  ),

  # best models
  tar_target(
    dry_het_intrain2_trait,
    generate_mcmc_summary(
      fit_summary_logistic_simple_stan_data_dry_het_intrain2_nlog,
      fit_mcmc_logistic_simple_stan_data_dry_het_intrain2_nlog,
      stan_data_dry_het_intrain2_nlog)
  ),
  tar_target(
    wet_phy_norain_trait,
    generate_mcmc_summary(
      fit_summary_logistic_simple_stan_data_wet_phy_norain_nlog,
      fit_mcmc_logistic_simple_stan_data_wet_phy_norain_nlog,
      stan_data_wet_phy_norain_nlog)
  ),
  tar_target(
    dry_het_intrain_abund,
    generate_mcmc_summary(
      fit_summary_logistic_simple_stan_data_dry_het_intrain_ab,
      fit_mcmc_logistic_simple_stan_data_dry_het_intrain_ab,
      stan_data_dry_het_intrain_ab)
  ),
  tar_target(
    wet_phy_rain_abund,
    generate_mcmc_summary(
      fit_summary_logistic_simple_stan_data_wet_phy_rain_ab,
      fit_mcmc_logistic_simple_stan_data_wet_phy_rain_ab,
      stan_data_wet_phy_rain_ab)
  ),

  # not best
  tar_target(
    wet_phy_intrain2_trait,
    generate_mcmc_summary(
      fit_summary_logistic_simple_stan_data_wet_phy_intrain2_nlog,
      fit_mcmc_logistic_simple_stan_data_wet_phy_intrain2_nlog,
      stan_data_wet_phy_intrain2_nlog)
  ),

  tar_target(
    dry_het_intrain2_trait_suv_contour_plot, {
      p <- dry_trait_suv_contour(dry_het_intrain2_trait, alpha = 0.05)
      my_ggsave(
        "figs/dry_het_intrain2_trait_suv_contour",
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
    wet_phy_intrain2_trait_suv_contour_plot, {
      p <- wet_trait_suv_contour(wet_phy_intrain2_trait, alpha = 0.05)
      my_ggsave(
        "figs/wet_phy_intrain2_trait_suv_contour",
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
    wet_phy_intrain2_trait_suv_contour_plot_cons, {
      p <- wet_trait_suv_contour(wet_phy_intrain2_trait, alpha = 0.05, keep_cons = TRUE)
      my_ggsave(
        "figs/wet_phy_intrain2_trait_suv_contour_cons",
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
    coef_trait_plot, {
      p <- coef_pointrange(dry_het_intrain2_trait, wet_phy_norain_trait, comb = FALSE)
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
    ggpairs_plot, {
      p <- my_ggpairs(traits_csv)
      my_ggsave(
        "figs/pairs",
        p,
        dpi = 300,
        width = 210,
        height = 210,
        units = "mm"
      )
    },
    format = "file"
  ),
  tar_target(
    abund, {
      p1 <- generate_suv_pred2(dry_het_intrain_abund$summary,
        dry_het_intrain_abund$data, alpha = 0.05, 2) |>
        subplot_fun(low = FALSE) +
        ggtitle("Abundant species") +
        xlab("Conspecific adult density")

      p2 <- generate_suv_pred2(dry_het_intrain_abund$summary,
        dry_het_intrain_abund$data, alpha = 0.05, 2) |>
        subplot_fun(low = TRUE) +
        ggtitle("Rare species") +
        xlab("Conspecific adult density")

      p <- p1 + p2

      my_ggsave(
        "figs/abund",
        p,
        dpi = 300,
        width = 173,
        height = 65,
        units = "mm"
      )
    },
    format = "file"
  ),


# # best models
# tar_map(
#   values = list(
#     x = rlang::syms(c(
#       "fit_summary_logistic_simple_stan_data_dry_het_intrain2_nlog",
#       "fit_summary_logistic_simple_stan_data_wet_phy_norain_nlog",
#       "fit_summary_logistic_simple_stan_data_wet_phy_intrain2_nlog",
#       "fit_summary_logistic_simple_stan_data_dry_het_intrain_ab",
#       "fit_summary_logistic_simple_stan_data_wet_phy_rain_ab")),
#     stan_data = rlang::syms(c(
#       "stan_data_dry_het_intrain2_nlog",
#       "stan_data_wet_phy_norain_nlog",
#       "stan_data_wet_phy_intrain2_nlog",
#       "stan_data_dry_het_intrain_ab",
#       "stan_data_wet_phy_rain_ab")),
#     path =
#       str_c(
#       "data/",
#        c("dry_het_intrain2_traits",
#         "wet_phy_norain_traits",
#         "wet_phy_intrain2_traits",
#         "dry_het_intrain_abund",
#         "wet_phy_rain_abund"),
#       "_gamma.csv")),
#   tar_target(
#     gamma_out_csv, {
#       create_gamma_tab(x, stan_data)  |>
#         my_write_csv(path)
#     },
#     format = "file"
#   )
# ),

  # tar_quarto(
  #   bayes_check_html,
  #   "docs/bayes_check.qmd",
  # ),
  # tar_quarto(
  #   si_pdf,
  #   "ms/SI.qmd"
  # ),
  # tar_quarto(
  #   main_docx,
  #   "ms/main.qmd"
  # ),

  NULL
 )

diagnostics_mapped <- tar_map(
    values = list(value = rlang::syms(str_c("fit_diagnostics_logistic_simple_stan_data_", data_names)), data_names = data_names),
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
    values = list(value = rlang::syms(str_c("fit_summary_logistic_simple_stan_data_", data_names)), data_names = data_names),
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
  "dry_het_intrain4_pc12",
  "dry_het_intrain2_nlog",
  "wet_het_rain_pc12",
  "wet_phy_norain_nlog",
  "dry_phy_intrain3_ab",
  "wet_phy_intrain3_ab"
)

# Construct 'x' and 'stan_data' lists based on the model names
x <- rlang::syms(str_c("fit_summary_logistic_simple_stan_data_", models))
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
  best_csv_mapped
)


list(data_, main_) #|>
  # append(util_list)
