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
