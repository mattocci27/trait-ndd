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

#' @title Coef plot
#' @param fit_gamma_dry summary data of gamma for dry season (e.g., fit1_gamma)
#' @param fit_gamma_wet summary data of gamma for wet season (e.g., fit2_gamma)
#' @param stan_data_dry data for stan model for dry season (e.g., dry_each+int)
#' @param stan_data_wet data for stan model for wet season (e.g., wet_each+int)
coef_pointrange <- function(fit_gamma_dry, fit_gamma_wet, int = TRUE) {
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

  fit_gamma_dry <- fit_gamma_dry |>
    mutate(sig = ifelse(q2_5 * q97_5 > 0, "sig", "ns")) |>
    mutate(season = "Dry") |>
    filter(trait_name == "intercept") |>
    mutate(para = factor(para, levels = paste0("gamma_", 1:nrow(fit_gamma_dry), "_1") |> rev()))

  fit_gamma_wet <- fit_gamma_wet |>
    mutate(sig = ifelse(q2_5 * q97_5 > 0, "sig", "ns")) |>
    mutate(season = "Rainy") |>
    filter(trait_name == "intercept") |>
    mutate(para = factor(para, levels = paste0("gamma_", 1:nrow(fit_gamma_wet), "_1") |> rev()))

  my_col <- brewer.pal(5, "RdBu")

  bind_rows(fit_gamma_dry, fit_gamma_wet) |>
    mutate(sig_season = paste0(season, "_", sig)) |>
    ggplot(aes(y = para)) +
    facet_grid(~season) +
    geom_vline(xintercept = 0, lty = 2, col = "grey40") +
    geom_linerange(
      aes(xmin = q2_5, xmax = q97_5, col = season),
    ) +
    geom_linerange(
      aes(xmin = q5, xmax = q95, col = season),
      size = 1.5,
    ) +
    geom_point(
      aes(x = mean_, fill = sig_season, col = season),
      shape = 21,
      size = 3) +
    scale_y_discrete(labels = gamma_lab) +
    scale_colour_manual(
      values = c(
        "Dry" = my_col[1],
        "Rainy" = my_col[5]
      )
    ) +
    scale_fill_manual(
      values = c(
        "Dry_sig" = my_col[2],
        "Rainy_sig" = my_col[4],
        "Dry_ns" = my_col[3],
        "Rainy_ns" = my_col[3]
      )
    ) +
    ylab("") +
    xlab("Standardized coefficients") +
    theme_bw() +
    theme(
      legend.position = "none")
}

#' @title Coef ridghe plot
#' @param fit_tab summary data (e.g., fit1_tab)
#' @param stan_data data for stan model (e.g., dry_each+int)
#' @param draws mcmc draws (e.g., fit_1_dry_each_int_draws_model_ind)
coef_ridge <- function(fit_tab, stan_dat, draws) {
  n_sp <- ncol(stan_dat$u)
  x_name <- colnames(stan_dat$x)

  fit1_beta <- fit_tab |>
    filter(str_detect(para, "beta")) |>
    mutate(pred_name = rep(x_name,  n_sp)) |>
    mutate(sp_name = rep(paste0("sp_", 1:n_sp), each = length(x_name)))


  tmp <- draws |>
      janitor::clean_names() |>
      dplyr::select(contains(c("beta")))

  tmp2 <- pivot_longer(tmp, 1:ncol(tmp)) |>
    group_by(name) |>
    nest()

  tmp3 <- tmp2 |>
    mutate(value = map(data, unlist)) |>
    mutate(density_= map(value, density)) |>
    mutate(x = map(density_, \(x)x$x)) |>
    mutate(y = map(density_, \(x)x$y)) |>
    mutate(max_y = map(y, max))

  tmp4 <- tmp3 |>
    ungroup() |>
    mutate(pred_name = rep(x_name,  n_sp)) |>
    mutate(sp_name = rep(paste0("sp_", 1:n_sp), each = length(x_name))) |>
    dplyr::select(name, x, y, max_y, pred_name, sp_name) |>
    unnest() |>
    mutate(density = y / max_y)

  ggplot(tmp4, aes(x = x, y = pred_name, gr = sp_name)) +
    geom_ridgeline(aes(height = density), fill = NA, colour = rgb(0, 0, 0, 0.1)) +
    theme_bw()
}

#' @title beta
