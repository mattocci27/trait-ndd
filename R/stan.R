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
