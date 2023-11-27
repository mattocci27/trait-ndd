my_ggsave <- function(filename, plot, units = c("in", "cm",
        "mm", "px"), height = NA, width = NA, dpi = 300, ...) {

  ggsave(
    filename = paste0(filename, ".png"),
    plot = plot,
    height = height,
    width = width,
    units = units,
    dpi = dpi,
    ...
  )

  ggsave(
    filename = paste0(filename, ".pdf"),
    plot = plot,
    height = height,
    width = width,
    units = units,
    dpi = dpi,
    ...
  )

  paste0(filename, c(".png", ".pdf"))
}

my_theme <- function(){
  theme_bw() %+replace%
  theme(
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.text.align = 0,
    # legend.key.height = unit(0.2, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
  )
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

#' @title Coef plot
coef_pointrange <- function(dry_trait, wet_trait, comb = TRUE) {
    if (comb) {
      data <- bind_rows(
        generate_coef_data(dry_trait$draws, dry_trait$data, season = "Dry"),
        generate_coef_data(wet_trait$draws, wet_trait$data, season = "Rainy"))
      coef_pointrange_wrapper(data)
    } else {
      dry_data <- generate_coef_data(dry_trait$draws, dry_trait$data, season = "Dry")
      wet_data <- generate_coef_data(wet_trait$draws, wet_trait$data, season = "Rainy")
      p1 <- coef_pointrange_wrapper(dry_data)
      p2 <- coef_pointrange_wrapper(wet_data)
      p1 + p2 +
        plot_annotation(tag_levels = "a") &
        theme(
          text = element_text(size = 8),
          plot.tag = element_text(face = "bold")
        )
    }
}

coef_pointrange_wrapper <- function(data) {
  gamma_lab <- c(
    # "(Intercept)" = expression(Intercept~(gamma["1,1"])),
    "logh_s" = expression(ln~Height~(gamma["2,1"])),
    "scon_s" = expression(ConS~(gamma["3,1"])),
    "shet_s" = expression(HetS~(gamma["4,1"])),
    "sphy_s" = expression(PhyS~(gamma["4,1"])),
    "acon_s_c" = expression(ConT~(gamma["5,1"])),
    "ahet_s_c" = expression(HetT~(gamma["6,1"])),
    "aphy_s_c" = expression(PhyT~(gamma["6,1"])),
    "rain_s" = expression(Rainfall~(gamma["7,1"])),
    "logh_s:rain_s" = expression(ln~Height%*%Rainfall~(gamma["8,1"])),
    "scon_s:rain_s" = expression(ConS%*%Rainfall~(gamma["9,1"])),
    "shet_s:rain_s" = expression(HetS%*%Rainfall~(gamma["10,1"])),
    "acon_s_c:rain_s" = expression(ConT%*%Rainfall~(gamma["11,1"])),
    "ahet_s_c:rain_s" = expression(HetT%*%Rainfall~(gamma["12,1"])))

  my_col <- RColorBrewer::brewer.pal(5, "RdBu")

  ggplot(data, aes(y = para)) +
    geom_vline(xintercept = 0, lty = 2, col = "grey60") +
    geom_linerange(aes(xmin = ll, xmax = hh, col = season),
      position = position_dodge(width = 0.5)) +
    geom_linerange(aes(xmin = l, xmax = h, col = season), size = 2,
      position = position_dodge(width = 0.5)) +
    geom_point(aes(x = m, col = season, fill = season_sig), shape = 21, size = 2.5,
      position = position_dodge(width = 0.5)) +
    scale_colour_manual(
      values = c(
        "Dry" = my_col[1],
        "Rainy" = my_col[5]
      )) +
    scale_fill_manual(
      values = c(
        "Dry_sig" = my_col[2],
        "Rainy_sig" = my_col[4],
        "Dry_ns" = my_col[3],
        "Rainy_ns" = my_col[3]
      )
    ) +
    scale_y_discrete(labels = gamma_lab) +
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

# #' @title generate beta plot data (takes time)
# generate_beta_list <- function(fit_beta, fit_gamma, stan_data, draws, x, y, x_lab, y_lab) {
#   tmp_beta <- fit_beta |>
#     filter(pred_name == paste(y))
#   gamma_slope <- fit_gamma |>
#     filter(pred_name == paste(y)) |>
#     filter(trait_name == paste(x))
#   gamma_int <- fit_gamma |>
#     filter(pred_name == paste(y)) |>
#     filter(trait_name == "intercept")
#   trait <- stan_data$u[paste(x), ]

#   tmp_para <- gamma_slope |>
#     pull(para) |>
#     str_split_fixed("_", 3)

#   beta_k <- tmp_beta$para |> str_split_fixed("_", 3)
#   k <- beta_k[, 2] |> unique()

#   y_lab_parse <- paste0("expression(", y_lab ,"(beta[paste(", k, ",',',", "j)]))")

#   beta_trait <- bind_cols(tmp_beta, trait = trait)

#   x_steps <- seq(min(beta_trait$trait), max(beta_trait$trait), length = 80)
#   new_data <- tibble(
#     observation = seq_along(x_steps) |> as.character(),
#     x_lt = x_steps)

# #  tic()
#   tmp <- draws |>
#     janitor::clean_names() |>
#     dplyr::select(
#       gamma_int |> pull(para),
#       gamma_slope |> pull(para))
# #  toc()

#   colnames(tmp) <- c("gamma_int_draw", "gamma_slope_draw")

#   tmp2 <- tmp |>
#     mutate(rep = paste0("rep", 1:nrow(tmp))) |>
#     nest(data = c(gamma_int_draw, gamma_slope_draw)) |>
#     mutate(x = list(new_data$x_lt)) |>
#     mutate(y = map2(x, data, \(x, data) {data$gamma_int_draw + data$gamma_slope_draw * x}))

#   pred_lin <- matrix(unlist(tmp2$y), ncol = nrow(draws), byrow = FALSE) |> t()
#   # dim(pred_lin)
#   df_pred_lin <- tidy_predictions(pred_lin, new_data)

#   list(
#     beta_trait = beta_trait,
#     df_pred_lin = df_pred_lin,
#     x_lab  = x_lab,
#     y_lab  = y_lab_parse
#     )
# }


# #' @title beta
# beta_plot <- function(list_data, rain = FALSE) {
#   my_col <- brewer.pal(5, "RdBu")
#   if (!rain) {
#     my_col1 <- my_col[1]
#     my_col2 <- my_col[2]
#   } else {
#     my_col1 <- my_col[5]
#     my_col2 <- my_col[4]
#   }
#   ggplot(list_data$beta_trait) +
#     geom_ribbon(
#       data = list_data$df_pred_lin,
#       aes(ymin = lower, ymax = upper, x = x_lt),
#       alpha = 0.4,
#       #fill = "grey60"
#       fill = my_col2
#     ) +
#     geom_line(
#       aes(y = mean, x = x_lt),
#       data = list_data$df_pred_lin,
#       colour = my_col1,
#       size = 0.5
#     ) +
#     geom_point(aes(x = trait, y = mean_),
#       colour = my_col1,
#       size = 0.5
#       ) +
#     geom_errorbar(aes(x = trait, ymin = q2_5, ymax = q97_5),
#       colour = my_col1,
#       size = 0.2
#       ) +
#     xlab(list_data$x_lab) +
#     ylab(eval(parse(text = list_data$y_lab))) +
#     theme_bw()
# }

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


beta_dry_comb_plot <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10){
  p1 + p2 + p3 + plot_spacer() +
  p4 + p5 + p6 + plot_spacer() +
  p7 + p8 + p9 + p10 + plot_layout(ncol = 4, nrow = 3) +
   plot_annotation(tag_levels = "A") &
   theme(
    text = element_text(size = 8),
     plot.tag = element_text(face = "bold")
    )
}

beta_wet_comb_plot <- function(p1, p2, p3, p4){
  p1 + p2 + p3 + p4 + plot_layout(ncol = 2, nrow = 2) +
   plot_annotation(tag_levels = "A") &
   theme(
    text = element_text(size = 8),
     plot.tag = element_text(face = "bold")
    )
}


logit <- function(z) 1 / (1 + exp(-z))

prepare_suv_pred <- function(mcmc_summary, mcmc_stan_data, alpha = c(0.05, 0.01, 0.5)) {
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
    # filter(sig == "sig") |>
    mutate(gamma_row_num = str_split_fixed(variable,  "\\[|\\]|,", 4)[, 2] |>
      as.numeric()) |>
    mutate(gamma_col_num = str_split_fixed(variable,  "\\[|\\]|,", 4)[, 3] |>
      as.numeric()) |>
    mutate(ind_pred = gamma_row[gamma_row_num]) |>
    mutate(sp_pred = gamma_col[gamma_col_num]) |>
    dplyr::select(variable, ind_pred, sp_pred, q50, sig)
}
# summary <- wet_trait$summary
# data <- wet_trait$data
# alpha <- 0.05
# trait_no <- 3
generate_suv_pred <- function(summary, data, alpha, trait_no, keep_cons = FALSE) {
  if (keep_cons) {
  tmp <- prepare_suv_pred(summary, data, alpha = alpha) |>
    mutate(sig = ifelse(ind_pred == "scon_s", "sig", sig)) |>
    mutate(q50 = ifelse(sig != "sig", 0, q50))
  } else {
    tmp <- prepare_suv_pred(summary, data, alpha = alpha) |>
      mutate(q50 = ifelse(sig != "sig", 0, q50))
  }
  k <- ncol(data$x)
  l <- nrow(data$u)
  gamma <- matrix(tmp$q50, nrow = k)
  trait_m <- matrix(c(1, rep(0, l - 1)))
  trait_l <- trait_m
  trait_h <- trait_m
  trait_l[trait_no, 1] <- qnorm(0.25)
  trait_h[trait_no, 1] <- qnorm(0.75)
#  dry_trait$data$x %>% str()

  scon_tmp <- data$x[, 3]
  rain_tmp <- data$x[, 7]
  median(scon_tmp)
  quantile(scon_tmp, 0.75)
  quantile(scon_tmp, 0.9)
  quantile(scon_tmp, 0.99)
  # max(scon_tmp)
  scon <- seq(min(scon_tmp), quantile(scon_tmp, 0.9), length = 20)
  rain <- seq(min(rain_tmp), max(rain_tmp), length = 20)

  h <- 0
  shet <- 0
  acon <- 0
  ahet <- 0
  grid_data <- expand_grid(scon, rain)

  beta_l <- gamma %*% trait_l
  beta_m <- gamma %*% trait_m
  beta_h <- gamma %*% trait_h

  grid_data2 <- grid_data |>
    mutate(suv_z_l = beta_l[1] +
      beta_l[3] * scon +
      beta_l[7] * rain +
      beta_l[8] * rain * scon) |>
    mutate(suv_z_m = beta_m[1] +
      beta_m[3] * scon +
      beta_m[7] * rain +
      beta_m[8] * rain * scon) |>
    mutate(suv_z_h = beta_h[1] +
      beta_h[3] * scon +
      beta_h[7] * rain +
      beta_h[8] * rain * scon) |>
    mutate(suv_p_l = logit(suv_z_l)) |>
    mutate(suv_p_m = logit(suv_z_m)) |>
    mutate(suv_p_h = logit(suv_z_h))
  grid_data2
}

subplot_fun <- function(data, low = TRUE) {
  if (low) {
    p <- ggplot(data, aes(x = scon, y = rain, fill = suv_p_l, z = suv_p_l))
  } else {
    p <- ggplot(data, aes(x = scon, y = rain, fill = suv_p_h, z = suv_p_h))
  }
  my_col <- viridis::viridis(3, option = "B")
  p +
    geom_tile() +
    geom_contour(col = "grey40") +
    scale_fill_viridis_c(option = "B")  +
    ylab("Rainfall") +
    xlab("Conspecific density") +
    labs(fill = "Probability of\nsurvival") +
    # scale_fill_gradient2(
    #   mid = my_col[2],
    #   low = my_col[1],
    #   high = my_col[3],
    #   midpoint = 0.87
    # ) +
    my_theme() #+
#    theme(legend.position = "none")
}

generate_suv_pred2 <- function(summary, data, alpha, trait_no, keep_cons = FALSE) {
  if (keep_cons) {
  tmp <- prepare_suv_pred(summary, data, alpha = alpha) |>
    mutate(sig = ifelse(ind_pred == "scon_s", "sig", sig)) |>
    mutate(q50 = ifelse(sig != "sig", 0, q50))
  } else {
    tmp <- prepare_suv_pred(summary, data, alpha = alpha) |>
      mutate(q50 = ifelse(sig != "sig", 0, q50))
  }
  k <- ncol(data$x)
  l <- nrow(data$u)
  gamma <- matrix(tmp$q50, nrow = k)
  trait_m <- matrix(c(1, rep(0, l - 1)))
  trait_l <- trait_m
  trait_h <- trait_m
  trait_l[trait_no, 1] <- qnorm(0.025)
  trait_h[trait_no, 1] <- qnorm(0.975)
#  dry_trait$data$x %>% str()

  scon_tmp <- data$x[, 5]
  rain_tmp <- data$x[, 7]
  median(scon_tmp)
  quantile(scon_tmp, 0.75)
  quantile(scon_tmp, 0.9)
  quantile(scon_tmp, 0.99)
  # max(scon_tmp)
  scon <- seq(min(scon_tmp), quantile(scon_tmp, 0.9), length = 20)
  rain <- seq(min(rain_tmp), max(rain_tmp), length = 20)

  h <- 0
  shet <- 0
  acon <- 0
  ahet <- 0
  grid_data <- expand_grid(scon, rain)

  beta_l <- gamma %*% trait_l
  beta_m <- gamma %*% trait_m
  beta_h <- gamma %*% trait_h

  grid_data2 <- grid_data |>
    mutate(suv_z_l = beta_l[1] +
      beta_l[5] * scon +
      beta_l[7] * rain +
      beta_l[11] * rain * scon) |>
    mutate(suv_z_m = beta_m[1] +
      beta_m[5] * scon +
      beta_m[7] * rain +
      beta_m[11] * rain * scon) |>
    mutate(suv_z_h = beta_h[1] +
      beta_h[5] * scon +
      beta_h[7] * rain +
      beta_h[11] * rain * scon) |>
    mutate(suv_p_l = logit(suv_z_l)) |>
    mutate(suv_p_m = logit(suv_z_m)) |>
    mutate(suv_p_h = logit(suv_z_h))
  grid_data2
}


dry_trait_suv_contour <- function(dry_trait, alpha = 0.05) {
  p1 <- generate_suv_pred(dry_trait$summary, dry_trait$data, alpha = 0.05, 2) |>
    subplot_fun(low = TRUE) +
    ggtitle("Low LDMC species")
  p2 <- generate_suv_pred(dry_trait$summary, dry_trait$data, alpha = 0.05, 2) |>
    subplot_fun(low = FALSE) +
    ggtitle("High LDMC species")
  p3 <- generate_suv_pred(dry_trait$summary, dry_trait$data, alpha = 0.05, 3) |>
    subplot_fun(low = TRUE) +
    ggtitle("Low SDMC species")
  p4 <- generate_suv_pred(dry_trait$summary, dry_trait$data, alpha = 0.05, 3) |>
    subplot_fun(low = FALSE) +
    ggtitle("High SDMC species")
  p5 <- generate_suv_pred(dry_trait$summary, dry_trait$data, alpha = 0.05, 4) |>
    subplot_fun(low = TRUE) +
    ggtitle(expression(Low~delta*C[13]~species))
  p6 <- generate_suv_pred(dry_trait$summary, dry_trait$data, alpha = 0.05, 4) |>
    subplot_fun(low = FALSE) +
    ggtitle(expression(High~delta*C[13]~species))
  p7 <- generate_suv_pred(dry_trait$summary, dry_trait$data, alpha = 0.05, 10) |>
    subplot_fun(low = TRUE) +
    ggtitle("Low LT species")
  p8 <- generate_suv_pred(dry_trait$summary, dry_trait$data, alpha = 0.05, 10) |>
    subplot_fun(low = FALSE) +
    ggtitle("High LT species")
  p9 <- generate_suv_pred(dry_trait$summary, dry_trait$data, alpha = 0.05, 6) |>
    subplot_fun(low = TRUE) +
    ggtitle(expression(Low~pi[tlp]))
  p10 <- generate_suv_pred(dry_trait$summary, dry_trait$data, alpha = 0.05, 6) |>
    subplot_fun(low = FALSE) +
    ggtitle(expression(High~pi[tlp]))

  (p1 + p2) / (p3 + p4) / (p5 + p6) / (p7 + p8) / (p9 + p10) +
     plot_annotation(tag_levels = "a") &
     theme(
      text = element_text(size = 8),
      plot.tag = element_text(face = "bold"))
}

wet_trait_suv_contour <- function(wet_trait, alpha = 0.05, keep_cons = FALSE) {
  p1 <- generate_suv_pred(wet_trait$summary, wet_trait$data, alpha = 0.05, 9, keep_cons) |>
    subplot_fun(low = TRUE) +
    ggtitle("Low N species")
  p2 <- generate_suv_pred(wet_trait$summary, wet_trait$data, alpha = 0.05, 9, keep_cons) |>
    subplot_fun(low = FALSE) +
    ggtitle("High N species")
  p3 <- generate_suv_pred(wet_trait$summary, wet_trait$data, alpha = 0.05, 6, keep_cons) |>
    subplot_fun(low = TRUE) +
    ggtitle(expression(Low~pi[tlp]~species))
  p4 <- generate_suv_pred(wet_trait$summary, wet_trait$data, alpha = 0.05, 6, keep_cons) |>
    subplot_fun(low = FALSE) +
    ggtitle(expression(High~pi[tlp]~species))

  (p1 + p2) / (p3 + p4) +
     plot_annotation(tag_levels = "a") &
     theme(
      text = element_text(size = 8),
      plot.tag = element_text(face = "bold"))
}

generate_beta_list <- function(draws, stan_data, x_lab, y_lab, ind_pred, sp_pred) {
  gamma_row <- stan_data$x |> colnames()
  gamma_col <- stan_data$u |> rownames()
  print(gamma_row[ind_pred])
  print(gamma_col[sp_pred])

  trait <- stan_data$u[sp_pred, ]

  beta <- draws |>
    dplyr::select(contains(str_c("beta[", ind_pred))) |>
    as.matrix()
  gamma <- draws |>
    dplyr::select(contains(str_c("gamma[", ind_pred))) |>
    as.matrix()

  par_res <- NULL
  for (i in 1:4000) {
    res <- beta[i, ] - gamma[i, ] %*% stan_data$u
    res <- as.numeric(res)
    pred <- gamma[i, 1] + gamma[i, sp_pred] * trait
    par_res_tmp <- res + pred
    par_res <- rbind(par_res, par_res_tmp)
  }

  res_data <- tibble(
    beta_m = apply(beta, 2, quantile, 0.5),
    res_m = apply(par_res, 2, quantile, 0.5),
    res_l = apply(par_res, 2, quantile, 0.025),
    res_h = apply(par_res, 2, quantile, 0.975)) |>
    mutate(trait)

  x_steps <- seq(min(trait), max(trait), length = 80)

  tmp2 <- NULL
  for (i in 1:4000) {
    tmp <- gamma[i, 1] + gamma[i, sp_pred] * x_steps
    tmp2 <- rbind(tmp2, tmp)
  }

  pred_data <- data.frame(
    trait = x_steps,
    ll = apply(tmp2, 2, quantile, 0.025),
    l = apply(tmp2, 2, quantile, 0.25),
    m = apply(tmp2, 2, quantile, 0.5),
    h = apply(tmp2, 2, quantile, 0.75),
    hh = apply(tmp2, 2, quantile, 0.975))

  y_lab_parse <- str_c("expression(", y_lab ,"~(beta[paste(", ind_pred, ",',',", "j)]))")

  list(
    pred_data = pred_data,
    res_data = res_data,
    x_lab  = x_lab,
    y_lab  = y_lab_parse
    )
}

beta_plot <- function(beta_list, partial = TRUE ) {
  my_col <- RColorBrewer::brewer.pal(5, "RdBu")
  if (!partial) {
    beta_list$res_data <- beta_list$res_data |>
      mutate(res_m = beta_m)
  }
  ggplot(beta_list$pred_data, aes(x = trait))  +
    geom_line(aes(y = m))  +
    geom_ribbon(aes(ymin = l, ymax = h), alpha = 0.5) +
    geom_ribbon(aes(ymin = ll, ymax = hh), alpha = 0.4) +
    geom_point(data = beta_list$res_data, aes(y = res_m), size = 0.5) +
    # geom_line(aes(y = m), col = my_col[1])  +
    # geom_ribbon(aes(ymin = l, ymax = h), fill = my_col[2], alpha = 0.8) +
    # geom_ribbon(aes(ymin = ll, ymax = hh), fill = my_col[2], alpha = 0.5) +
    # geom_point(data = beta_list$res_data, aes(y = res_m), col = my_col[1]) +
    # geom_errorbar(data = beta_list$res_data, aes(ymin = res_l, ymax = res_h)) +
    ylab(eval(parse(text = beta_list$y_lab))) +
    xlab(beta_list$x_lab) +
    theme_bw()
    # my_theme()
}


cc_line <- function(wet, dry) {
  wet_cc_het <- which(wet$het == max(wet$het)) / nrow(wet)
  dry_cc_het <- which(dry$het == max(dry$het)) / nrow(dry)

  plot_fun <- function(data) {
    ggplot(data, aes(x = cc, y = phy)) +
      geom_line() +
      xlab("c") +
      ylab("Log-likelihood") +
      my_theme() +
      theme(
        axis.title.x = element_text(face = "italic")
      )
  }

  p1 <- plot_fun(dry) +
    geom_vline(xintercept = dry_cc_het, linetype = "dashed") +
    ggtitle("Dry season")
  p2 <- plot_fun(wet) +
    geom_vline(xintercept = wet_cc_het, linetype = "dashed") +
    ggtitle("Rainy season")

  p1 + p2 +
    plot_annotation(tag_levels = "a")
}


# library(tidyverse)
# library(here)
# library(GGally)
# library(tictoc)
# library(targets)

# tar_load(trait_csv)

my_ggpairs <- function(trait_csv) {
  trait <- read_csv(trait_csv) #|>

  trait2 <- trait |>
    mutate(log_LA = log(LA)) |>
    mutate(log_SLA = log(SLA)) |>
    mutate(log_LT = log(LT)) |>
    mutate(log_N = log(N)) |>
    rename(SD = WD) |>
    dplyr::select(
      LDMC,
      SD,
      SDMC,
      C13,
      C,
      log_N,
      CN,
      tlp,
      log_LA,
      log_SLA,
      log_LT)

# Custom function to remove axes lines
  my_lower <- function(data, mapping, ...){
    p <- ggplot(data = data, mapping = mapping) +
      geom_point() +
      theme_bw()  # removes axes lines, labels, and background
      # theme(plot.margin = margin(0, 0, 0, 0)) # removes margins
    return(p)
  }

#Custom function to modify 'cor' plots
  custom_cor <- function(data, mapping, ...){
    # Create the default 'cor' plot
    p <- ggally_cor(data, mapping, size = 3, ...)

    # Modify the theme to remove axes lines, labels, and background
    p <- p +
      theme_void() + # removes axes lines, labels, and background
      theme(plot.margin = margin(0, 0, 0, 0))

    return(p)
  }

  custom_density <- function(data, mapping, ...){
    # Create the default 'density' plot
    p <- ggally_densityDiag(data, mapping, ...)

    # Modify the theme for the diagonal elements
    p <- p +
      theme_bw() + # example theme
      theme(plot.margin = margin(0, 0, 0, 0)) # removes margins

    return(p)
  }

  trait2 |>
    ggpairs(
     upper = list(continuous = custom_cor),
     diag = list(continuous = custom_density),
     lower = list(continuous = my_lower)) +
    theme(
      strip.background = element_blank(),
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 0.8)
    )
}


pca_panel <- function(traits_df) {

  traits <- traits_df |>
    mutate(la = log(la)) |>
    mutate(lt = log(lt)) |>
    mutate(n = log(n)) |>
    mutate(sla = log(sla)) |>
    rename(log_la = la) |>
    rename(log_lt = lt) |>
    rename(log_n = n) |>
    rename(log_sla = sla)

  pca <- PCA(traits[, 2:13], graph = FALSE)

  p_eig <- fviz_eig(pca)

# Convert the variable coordinates to a data frame
  var_coords <- as.data.frame(get_pca_var(pca)$coord)
  rownames(var_coords) <- rownames(get_pca_var(pca)$coord)

  p12 <- fviz_pca_var(pca, geom = "", arrows = TRUE) +
    geom_segment(data = var_coords,
      aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
      arrow = arrow(type = "closed", length = unit(0.05, "inches")),
      linewidth = 0.5,
      color = "black")  +
    geom_text_repel(
      data = var_coords,
      aes(x = Dim.1, y = Dim.2, label = rownames(var_coords)),
      size = 4,    # Adjust text size
      # box.padding = unit(0.35, "lines"),   # Adjust box padding
      # point.padding = unit(0.3, "lines"),  # Adjust point padding
      segment.color = 'grey50',
      seed = 123
    )

  p_eig + p12 +
    plot_annotation(tag_levels = "a")

}


phy_het_points <- function(seedling_df) {
  formula <- y ~ x  # Define the formula for the regression line

  seedling_df <- seedling_df |>
    mutate(season = if_else(season == "dry", "Dry", "Rainy"))

  p1 <- ggplot(seedling_df, aes(x = shet, y = sphy)) +
    geom_pointdensity() +
    scale_color_viridis_c() +
    stat_smooth(method = "lm", formula = formula, se = FALSE, size = 0.5) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 formula = formula, parse = TRUE) +
    facet_grid(~season) +
    xlab("Heterospecific seedling density") +
    ylab("Phylogenetic seedling density")

  p2 <- ggplot(seedling_df, aes(x = ahet, y = aphy)) +
    geom_pointdensity() +
    scale_color_viridis_c() +
    stat_smooth(method = "lm", formula = formula, se = FALSE, size = 0.5) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 formula = formula, parse = TRUE) +
    facet_grid(~season) +
    xlab("Heterospecific tree density") +
    ylab("Phylogenetic tree density")

  p1 / p2 +
    plot_annotation(tag_levels = "a") &
    my_theme() &
    theme(
      legend.position = "none"
    )

}


