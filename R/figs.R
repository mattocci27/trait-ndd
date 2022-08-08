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

#' @title generate beta plot data (takes time)
generate_beta_list <- function(fit_beta, fit_gamma, stan_data, draws, x, y, x_lab, y_lab) {
  tmp_beta <- fit_beta |>
    filter(pred_name == paste(y))
  gamma_slope <- fit_gamma |>
    filter(pred_name == paste(y)) |>
    filter(trait_name == paste(x))
  gamma_int <- fit_gamma |>
    filter(pred_name == paste(y)) |>
    filter(trait_name == "intercept")
  trait <- stan_data$u[paste(x), ]

  tmp_para <- gamma_slope |>
    pull(para) |>
    str_split_fixed("_", 3)

  beta_k <- tmp_beta$para |> str_split_fixed("_", 3)
  k <- beta_k[, 2] |> unique()

  y_lab_parse <- paste0("expression(", y_lab ,"(beta[paste(", k, ",',',", "j)]))")

  beta_trait <- bind_cols(tmp_beta, trait = trait)

  x_steps <- seq(min(beta_trait$trait), max(beta_trait$trait), length = 80)
  new_data <- tibble(
    observation = seq_along(x_steps) |> as.character(),
    x_lt = x_steps)

#  tic()
  tmp <- draws |>
    janitor::clean_names() |>
    dplyr::select(
      gamma_int |> pull(para),
      gamma_slope |> pull(para))
#  toc()

  colnames(tmp) <- c("gamma_int_draw", "gamma_slope_draw")

  tmp2 <- tmp |>
    mutate(rep = paste0("rep", 1:nrow(tmp))) |>
    nest(data = c(gamma_int_draw, gamma_slope_draw)) |>
    mutate(x = list(new_data$x_lt)) |>
    mutate(y = map2(x, data, \(x, data) {data$gamma_int_draw + data$gamma_slope_draw * x}))

  pred_lin <- matrix(unlist(tmp2$y), ncol = nrow(draws), byrow = FALSE) |> t()
  # dim(pred_lin)
  df_pred_lin <- tidy_predictions(pred_lin, new_data)

  list(
    beta_trait = beta_trait,
    df_pred_lin = df_pred_lin,
    x_lab  = x_lab,
    y_lab  = y_lab_parse
    )
}


#' @title beta
  #          x = "chl"
  #          y = "cons_scaled"
  #          targets::tar_load(dry_each_int_s)
  #         stan_data = dry_each_int_s
  # targets::tar_load(fit_9_dry_each_int_s_draws_model_ind)
  # targets::tar_load(fit9_gamma)
  # fit_gamma <- fit9_gamma
  # targets::tar_load(fit9_beta)
  # fit_beta <- fit9_beta
  # draws <- fit_9_dry_each_int_s_draws_model_ind
  #         x_lab = "SDMC"
  #         y_lab  = "ConS~effect~"
beta_plot <- function(list_data) {
  ggplot(list_data$beta_trait) +
    geom_point(aes(x = trait, y = mean_)) +
    geom_errorbar(aes(x = trait, ymin = q2_5, ymax = q97_5)) +
    geom_ribbon(
      data = list_data$df_pred_lin,
      aes(ymin = lower, ymax = upper, x = x_lt),
      alpha = 0.4,
      fill = "grey60"
    ) +
    geom_line(
      aes(y = mean, x = x_lt),
      data = list_data$df_pred_lin,
      # colour = "#3366FF",
      size = 1
    ) +
    xlab(list_data$x_lab) +
    ylab(eval(parse(text = list_data$y_lab_parse))) +
    theme_bw()
}

#' @title predicitons for beta
#' @para mat_pred dataframe for MCMC draws of mean predictions (n.mcmc x 80)
#' @para df_pred dataframe for trait (80 x 2)
#' @ref https://www.tjmahr.com/visualizing-uncertainty-rstanarm/
tidy_predictions <- function(
  mat_pred,
  df_data,
  obs_name = "observation",
  prob_lwr = .025,
  prob_upr = .975
) {
  # Get dataframe with one row per fitted value per posterior sample
  df_pred <- mat_pred |>
    as_tibble() |>
    setNames(seq_len(ncol(mat_pred))) |>
    tibble::rownames_to_column("posterior_sample") |>
    tidyr::pivot_longer(
      cols = c(-posterior_sample),
      names_to = obs_name,
      values_to = "fitted"
    )

  # Helps with joining later
  class(df_pred[[obs_name]]) <- class(df_data[[obs_name]])

  # Summarise prediction interval for each observation
  df_pred |>
    group_by(.data[[obs_name]]) |>
    summarise(
      mean = mean(fitted),
      lower = quantile(fitted, prob_lwr),
      upper = quantile(fitted, prob_upr)
    ) |>
    left_join(df_data, by = obs_name)

}


beta_comb_plot <- function(p1, p2, p3){
  (p1 + p2 + p3 + plot_spacer()) /
  (p4 + p5 + p6 + plot_spacer()) /
  (p7 + p8 + p9 + p10)
}
