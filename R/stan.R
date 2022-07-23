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
  lwr2_5 <- apply(tmp, 2, \(x)(quantile(x, 0.025)))
  lwr5 <- apply(tmp, 2, \(x)(quantile(x, 0.05)))
  upr97_5 <- apply(tmp, 2, \(x)(quantile(x, 0.975)))
  upr95 <- apply(tmp, 2, \(x)(quantile(x, 0.9)))
  tibble(para = names(mean_), mean_, lwr2_5, lwr5, upr95, upr97_5)
}

