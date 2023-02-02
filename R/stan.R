#' @title Check divergence from draws
div_check <- function(diags) {
  n1 <- diags |>
    filter(divergent__ == 1) |>
    nrow()
  n2 <- diags |>
    nrow()
  print(paste(
    n1, "of", n2,
    "iterations ended with a divergence", n1 / n2 * 100, "%"
  ))
}

#' @title Compile a Stan model and return a path to the compiled model output.
#' @description We return the paths to the Stan model specification
#'   and the compiled model file so `targets` can treat them as
#'   dynamic files and compile the model if either file changes.
#' @return Path to the compiled Stan model, which is just an RDS file.
#'   To run the model, you can read this file into a new R session with
#'   `readRDS()` and feed it to the `object` argument of `sampling()`.
#' @param model_file Path to a Stan model file.
#'   This is a text file with the model spceification.
#' @references https://github.com/wlandau/targets-stan
#' @examples
#' library(cmdstanr)
#' compile_model("stan/model.stan")
compile_model <- function(model_file) {
  quiet(cmdstan_model(model_file))
  model_file
}

#' @title Suppress output and messages for code.
#' @description Used in the pipeline.
#' @return The result of running the code.
#' @param code Code to run quietly.
#' @references https://github.com/wlandau/targets-stan
#' @examples
#' library(cmdstanr)
#' library(tidyverse)
#' compile_model("stan/model.stan")
#' quiet(fit_model("stan/model.stan", simulate_data_discrete()))
#' out
quiet <- function(code) {
  sink(nullfile())
  on.exit(sink())
  suppressMessages(code)
}

my_loo <- function(x) x$loo(cores = parallel::detectCores())

#' @title Create summary table for posteriors
create_stan_tab <- function(draws) {
  tmp <- draws |>
    janitor::clean_names() |>
    dplyr::select(contains(c("beta", "gamma")))
  mean_ <- apply(tmp, 2, mean)
  q2_5 <- apply(tmp, 2, \(x)(quantile(x, 0.025)))
  q5 <- apply(tmp, 2, \(x)(quantile(x, 0.05)))
  q97_5 <- apply(tmp, 2, \(x)(quantile(x, 0.975)))
  q95 <- apply(tmp, 2, \(x)(quantile(x, 0.9)))
  tibble(para = names(mean_), mean_, q2_5, q5, q95, q97_5)
}

#' @title clean tabs
create_gamma_tab <- function(fit_tab, stan_dat) {
  x_name <- colnames(stan_dat$x)
  u_name <- rownames(stan_dat$u)
  fit_tab |>
    filter(str_detect(para, "gamma")) |>
    mutate(pred_name = rep(x_name,  length(u_name))) |>
    mutate(trait_name = rep(u_name, each = length(x_name)))# |>
}
#' @title clean tabs
create_beta_tab <- function(fit_tab, stan_dat) {
  x_name <- colnames(stan_dat$x)
  sp_name <- paste0("sp_",1:stan_dat$J)
  fit_tab |>
    filter(str_detect(para, "beta")) |>
    mutate(pred_name = rep(x_name, length(sp_name))) |>
    mutate(sp_name = rep(sp_name, each = length(x_name)))
}

coef_pointrange0 <- function(fit_gamma, stan_data, title = "Dry", int = TRUE) {

  gamma_lab <- c(
      "gamma_1_1" = expression(Intercept~(gamma["1,1"])),
      "gamma_2_1" = expression(ln~Height~(gamma["2,1"])),
      "gamma_3_1" = expression(ConS~(gamma["3,1"])),
      "gamma_4_1" = expression(ConT~(gamma["4,1"])),
      "gamma_5_1" = expression(HetS~(gamma["5,1"])),
      "gamma_6_1" = expression(HetT~(gamma["6,1"])),
      "gamma_7_1" = expression(Rainfall~(gamma["7,1"])),
      "gamma_8_1" = expression(ConS%*%Rainfall~(gamma["8,1"])),
      "gamma_9_1" = expression(ConT%*%Rainfall~(gamma["9,1"])),
      "gamma_10_1" = expression(HetS%*%Rainfall~(gamma["10,1"])),
      "gamma_11_1" = expression(HetT%*%Rainfall~(gamma["11,1"])))

  if (!int) {
   gamma_lab <- gamma_lab[1:7]
  }

  fig_dat <- fit_gamma |>
    filter(trait_name == "intercept") |>
    mutate(sig = ifelse(q2_5 * q97_5 > 0, "sig", "ns"))

  fig_dat |>
    mutate(para = factor(para, levels = paste0("gamma_", 1:nrow(fig_dat), "_1") |> rev())) |>
    ggplot(aes(y = para)) +
    geom_linerange(
      aes(xmin = q2_5, xmax = q97_5),
      color = "#3366FF"
    ) +
    geom_linerange(
      aes(xmin = q5, xmax = q95),
      size = 1.5,
      color = "#3366FF"
    ) +
    geom_point(
      aes(x = mean_, fill = sig),
      color = "#3366FF",
      shape = 21,
      size = 3) +
    scale_y_discrete(labels = gamma_lab) +
    scale_fill_manual(
        values = c(
          "sig" = "#33CCFF",
          # "sig95" = "#33CCFF",
          # "sig90" = "grey",
          "ns" = "#FFFFFF"
        )) +
    ylab("") +
    xlab("Standardized coefficients") +
    ggtitle(paste(title))
}

print_summary_tbl <- function(mcmc_summary, mcmc_stan_data, alpha = c(0.05, 0.01, 0.5)) {
  gamma_row <- mcmc_stan_data$x |> colnames()
  gamma_col <- mcmc_stan_data$u |> rownames()

  if (alpha == 0.05) {
    mcmc_summary <- mcmc_summary |>
     mutate(sig = ifelse(q2.5 * q97.5 > 0, "sig", "ns"))
  } else if (alpha == 0.1) {
    mcmc_summary <- mcmc_summary |>
     mutate(sig = ifelse(q5 * q95 > 0, "sig", "ns"))
  } else if (alpha == 0.5) {
    mcmc_summary <- mcmc_summary |>
     mutate(sig = ifelse(q25 * q75 > 0, "sig", "ns"))
  }

  mcmc_summary |>
    filter(str_detect(variable, "gamma")) |>
    filter(sig == "sig") |>
    mutate(gamma_row_num = str_split_fixed(variable,  "\\[|\\]|,", 4)[, 2] |>
      as.numeric()) |>
    mutate(gamma_col_num = str_split_fixed(variable,  "\\[|\\]|,", 4)[, 3] |>
      as.numeric()) |>
    mutate(ind_pred = gamma_row[gamma_row_num]) |>
    mutate(sp_pred = gamma_col[gamma_col_num]) |>
    dplyr::select(variable, ind_pred, sp_pred, q2.5, q5, q25, q50, q75, q95, q97.5) |>
    kbl() |>
    kable_styling(bootstrap_options = c("striped", "HOLD_position"))
}

load_mcmc_summary <- function(loo_tbl, season = "dry", trait = "ab") {
  tmp <- loo_tbl |>
    filter(season == {{season}})
  tmp
  if (trait == "ab") {
    tmp <- tmp |>
      filter(str_detect(model, "ab|ba"))
  } else {
    tmp <- tmp |>
      filter(str_detect(model, "_nlog$"))
  }
  tmp <- tmp |>
    arrange(-elpd) |>
    pull(model)
  tmp_summary <- str_replace(tmp[1], "loo_fit_mcmc", "fit_summary")
  tmp_data <- str_replace(tmp_summary, "fit_summary_logistic_simple_", "")
  tmp_draws <- str_replace(tmp[1], "loo_fit_mcmc", "fit_draws")
  withr::with_dir(rprojroot::find_root('_targets.R'),
    targets::tar_load(tmp_summary))
  withr::with_dir(rprojroot::find_root('_targets.R'),
    targets::tar_load(tmp_data))
  withr::with_dir(rprojroot::find_root('_targets.R'),
    targets::tar_load(tmp_draws))
  list(
    name = tmp[1],
    data = get(tmp_data),
    summary = get(tmp_summary),
    draws = get(tmp_draws)
  )
}


generate_coef_data <- function(draws, data, season = "Dry") {
  ind_pred_tmp <- data$x |> colnames()
  sp_pred_tmp <- data$u |> rownames()

  intervals_data <- mcmc_intervals_data(
    draws,
    regex_pars = "gamma",
    point_est = "median",
    prob = 0.5,
    prob_outer = 0.95) |>
    mutate(ind_pred = rep(ind_pred_tmp, length(sp_pred_tmp))) |>
    mutate(sp_pred = rep(sp_pred_tmp, each = length(ind_pred_tmp))) |>
    mutate(season = season) |>
    mutate(sig = ifelse(ll * hh > 0, "sig", "ns")) |>
    mutate(season_sig = paste0(season, "_", sig)) |>
    filter(str_detect(parameter, "1\\]$")) |>
    filter(parameter != "gamma[1,1]")

    intervals_data |>
      mutate(para = factor(ind_pred,
        levels = c(
          "logh_s",
          "scon_s",
          "shet_s",
          "acon_s_c",
          "ahet_s_c",
          "rain_s",
          "logh_s:rain_s",
          "scon_s:rain_s",
          "shet_s:rain_s",
          "acon_s_c:rain_s",
          "ahet_s_c:rain_s"
        ) |> rev()))
}

generate_loo_tbl <- function(loo_list)  {
  loo_list_ori <- loo_list
  loo_list <- loo_list_ori[str_detect(names(loo_list_ori), "het")]
  loo_list <- loo_list[!str_detect(names(loo_list), "_n$")]
  loo_list <- loo_list[str_detect(names(loo_list), "simple")]
  loo_list <- loo_list[!str_detect(names(loo_list), "wet_het_rain_nlog")]
  loo_list <- loo_list[!str_detect(names(loo_list), "wet_het_norain_nlog")]
  loo_names <- names(loo_list)

  loo_names_split <- str_split_fixed(loo_names, "_", 11)

  loo_tbl <- tibble(model = names(loo_list)) |>
    mutate(elpd = map_dbl(loo_list, \(x)x$elpd_loo)) |>
    mutate(p_loo = map_dbl(loo_list, \(x)x$p_loo)) |>
    mutate(looic = map_dbl(loo_list, \(x)x$looic)) |>
    mutate(season = loo_names_split[, 8]) |>
    mutate(phy = loo_names_split[, 9]) |>
    mutate(rain = loo_names_split[, 10]) |>
    mutate(traits = loo_names_split[, 11])
  loo_tbl
}
