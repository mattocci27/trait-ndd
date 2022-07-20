
targets::tar_load(fit_1_dry_full_draws_model_ind)
targets::tar_load(fit_2_wet_full_draws_model_ind)

dry_data <- create_stan_tab(fit_1_dry_full_draws_model_ind) |>
  mutate(season = "Dry")
wet_data <- create_stan_tab(fit_2_wet_full_draws_model_ind) |>
  mutate(season = "Rainy")


library(tidyverse)

draws <- fit_1_dry_full_draws_model_ind |>
    janitor::clean_names() |>
    dplyr::select(contains(c("beta", "gamma")))


coef_pointrange <- function(data, ld = TRUE, add = FALSE) {

  data <- wet_data |>
    mutate(sig = ifelse(lwr2_5 * upr97_5 > 0, "sig", "ns")) |>
    mutate(ci_sig = case_when(
      lwr2_5 * upr97_5 > 0 ~ "sig95",
      lwr5 * upr95 > 0 ~ "sig90",
      TRUE ~ "ns"
      )) |> filter(sig == "sig") |>
    filter(str_detect(para, "gamma"))


  # beta (mean)
  data1 <- dry_data |>
    filter(str_detect(para, "gamma")) |>
    #filter(para != "beta_1") |>
    mutate(para = factor(para, levels = rev(para)))
  # gamma (sd)

  data2 <- data |>
    filter(str_detect(para, "gamma")) |>
        "hets_scaled",
        "heta_scaled_c",
        "rain_scaled",
        "cons_rain",
        "cona_rain",
        "hets_rain",
        "heta_rain"
        ) |> rev() )
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
  plot_fun(data1) +
    scale_y_discrete(labels = gamma_lab)

  targets::tar_load(dry_full)
  colnames(dry_full$x)


  # plot function without title and y-lab
  plot_fun <- function(data) {
    data |>
      ggplot(aes(y = para_fct)) +
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
      ggtitle("Effects on mean") +
      scale_fill_manual(
        values = c(
          "sig" = "#33CCFF",
          # "sig95" = "#33CCFF",
          # "sig90" = "grey",
          "ns" = "#FFFFFF"
        )) +
      xlab("Standardized coefficients")
  }


data <- create_stan_tab(draws)


targets::tar_load(dry_full)
colnames(dry_full$x)

targets::tar_load(dry_wd)
targets::tar_load(dry_pca)

str(dry_full)
str(dry_wd)
str(dry_pca)


model_code <-
  '
  functions {
    vector test(matrix beta, matrix x) {
      int N = rows(x);
      vector[N] y;
      int sp[N];
      sp[1] = 1;
      sp[2] = 1;
      sp[3] = 2;
      sp[4] = 2;
      for (n in 1:N) {
        y[n] = x[n, ] * beta[, sp[n]];
      }
      return y;
   }
    vector test2(matrix beta, matrix x) {
      int N = rows(x);
      vector[N] y;
      int sp[N];
      sp[1] = 1;
      sp[2] = 1;
      sp[3] = 2;
      sp[4] = 2;
      y = x * beta[,sp];
      return y;
   }
  }
'


expose_stan_functions(stanc(model_code = model_code))

n <- 4
k <- 3
j <- 2
x <- matrix(1:12, n, k)
beta <- matrix(1:6, k, j)
test(beta, x)

x[2,] %*% beta[,1]

sp <- rep(1:2, each = n / 2)

