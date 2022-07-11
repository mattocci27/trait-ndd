#' @title Create summary table for posteriors
create_stan_tab <- function(draws) {
  tmp <- draws |>
    janitor::clean_names() |>
    dplyr::select(contains(c("beta", "gamma")))
  mean_ <- apply(tmp, 2, mean)
  lwr2_5 <- apply(tmp, 2, \(x)(quantile(x, 0.025)))
  lwr5 <- apply(tmp, 2, \(x)(quantile(x, 0.05)))
  upr97_5 <- apply(tmp, 2, \(x)(quantile(x, 0.975)))
  upr95 <- apply(tmp, 2, \(x)(quantile(x, 0.9)))
  tibble(para = names(mean_), mean_, lwr2_5, lwr5, upr95, upr97_5)
}

#' @title use mean estiamte
#' @param ld use leaf density (TRUE) or LMAd (FALSE)
#' @param add use oringal coef (FALSE) or added coef (TRUE)
coef_pointrange <- function(dry_data, wet_data) {
  dry_data <- create_stan_tab(dyr_data) |>
    mutate(season = "Dry")
  wet_data <- create_stan_tab(wet_data) |>
    mutate(season = "Rainy")
  data <- bind_rows(dry_data, wet_data) |>
    mutate(sig = ifelse(lwr2_5 * upr97_5 > 0, "sig", "ns")) |>
    mutate(ci_sig = case_when(
      lwr2_5 * upr97_5 > 0 ~ "sig95",
      lwr5 * upr95 > 0 ~ "sig90",
      TRUE ~ "ns"
      ))
  data <- data |>
    filter(str_detect(para, "gamma")) |>
    filter(str_detect(para, "_1$")) |>
    filter(para != "gamma_1_1") |>
    #mutate(para = factor(para, levels = rev(para))) |>
    mutate(para_chr = colnames(dry_full$x)[-1] |> rep(2)) |>
    mutate(para_fct = factor(para_chr,
      levels = c(
        "logh_scaled",
        "cons_scaled",
        "cona_scaled_c",
        "hets_scaled",
        "heta_scaled_c",
        "rain_scaled",
        "cons_rain",
        "cona_rain",
        "hets_rain",
        "heta_rain"
        ) |> rev())
    )

  gamma_lab <- c(
    cons_scaled = expression(ConS),
    hets_scaled = expression(HetS),
    cona_scaled_c = expression(ConT),
    heta_scaled_c = expression(HetT),
    cons_rain = expression(ConS%*%Rainfall),
    cona_rain = expression(ConT%*%Rainfall),
    hets_rain = expression(HetS%*%Rainfall),
    heta_rain = expression(HetT%*%Rainfall),
    rain_scaled = expression(Rainfall),
    logh_scaled = expression(Height)
  )
  # plot function without title and y-lab
  plot_fun <- function(data) {
    data |>
      ggplot(aes(y = para_fct)) +
      facet_grid(~ season) +
      geom_vline(xintercept = 0, lty  = 2, color = "grey60") +
      geom_linerange(
        aes(xmin = lwr2_5, xmax = upr97_5),
        color = "#3366FF") +
      geom_linerange(
        aes(xmin = lwr5, xmax = upr95),
        size = 1.5,
        color = "#3366FF") +
      geom_point(
        aes(x = mean_, fill = sig),
        color = "#3366FF",
        shape = 21,
        size = 3) +
      ylab("") +
    #  ggtitle("Effects on mean") +
      scale_fill_manual(
        values = c(
          "sig" = "#33CCFF",
          # "sig95" = "#33CCFF",
          # "sig90" = "grey",
          "ns" = "#FFFFFF"
        )) +
      xlab("Standardized coefficients")
  }
  plot_fun(data) +
    scale_y_discrete(labels = gamma_lab)
}
