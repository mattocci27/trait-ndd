#' @inheritParams readr::write_csv
my_write_csv <- function(x, path, append = FALSE, col_names = !append) {
    write_csv(x, path, append = FALSE, col_names = !append)
    paste(path)
}

clean_data <- function(rda, write_trait = FALSE) {
  load(rda)
  trait2 <- left_join(trait, community) |>
    filter(!is.na(abundance)) |>
    arrange(latin) |>
    rename(LDMC = LDMCsp) |>
    rename(WD = WDsp) |>
    rename(SDMC = SDMCsp) |>
    rename(LA = Lasp) |>
    rename(SLA = SLAsp) |>
    rename(Chl = Chlsp) |>
    rename(LT = LTsp) |>
    rename(CN = CNratio) |>
    rename(ab = abundance) |>
    rename(BA = basal)

  trait2[trait2$latin == "Alchornea_tiliifolia", "SLA"] <- 189.5

  dataseedling2 <- dataseedling |>
    filter(latin %in% trait2$latin) |>
    mutate(aphy = aphy |> as.numeric()) |>
    mutate(sphy = sphy |> as.numeric()) |>
    dplyr::select(quadrat, tag, latin, h1:shet, aphy:sphy, RF) |>
    dplyr::select(-note1, -dat)

  if (write_trait) {
    my_write_csv(trait2, "data/traits.csv")
  } else {
    my_write_csv(dataseedling2, "data/seedling.csv")
  }
}

calc_scale_cc <- function(seedling_csv, wet = TRUE) {
  seedling <- read_csv(seedling_csv) |>
    janitor::clean_names()
  if (wet) {
    seedling <- seedling |>
      filter(season == "rainy")
  } else {
    seedling <- seedling |>
      filter(season == "dry")
  }
  x1 <- seedling$acon
  x2 <- seedling$ahet
  y <- seedling$surv
  n_len <- 100
  lik <- numeric(n_len)
  lik2 <- numeric(n_len)
  for (i in 1:n_len) {
    d1 <- x1^(i / n_len)
    d2 <- x2^(i / n_len)
    fm1 <- glm(y ~ d1 + d2, family = binomial)
    lik[i] <- logLik(fm1)
  }
  cc <- which(lik == max(lik)) / n_len
  cc
}

# This is for visualization
generate_cc_data <- function(seedling_csv, wet = TRUE) {
  seedling <- read_csv(seedling_csv) |>
    janitor::clean_names()
  if (wet) {
    seedling <- seedling |>
      filter(season == "rainy")
  } else {
    seedling <- seedling |>
      filter(season == "dry")
  }
  x1 <- seedling$acon
  x2 <- seedling$ahet
  z1 <- seedling$aphy
  y <- seedling$surv

  n_len <- 100

  lik <- numeric(n_len)
  lik2 <- numeric(n_len)
  for (i in 1:n_len) {
    d1 <- x1^(i / n_len)
    d2 <- x2^(i / n_len)
    d3 <- z1^(i / n_len)
    fm1 <- glm(y ~ d1 + d2, family = binomial)
    fm2 <- glm(y ~ d1 + d3, family = binomial)
    lik[i] <- logLik(fm1)
    lik2[i] <- logLik(fm2)
  }
  tibble(cc = (1:n_len) / n_len, het = lik, phy = lik2)
}

generate_pca_data <- function(traits_df) {
  traits <- traits_df |>
    mutate(
      log_la = log(la),
      log_lt = log(lt),
      log_n = log(n),
      log_sla = log(sla)
    ) |>
    select(-la, -lt, -n, -sla)

  # use FactoMineR
  pca <- traits |>
    dplyr::select(-ab, -ba, -latin) |>
    PCA(graph = FALSE)

  bind_cols(traits, pca$ind$coord) |>
    janitor::clean_names() #|>
    # my_write_csv("data/trait_pca.csv")
}

write_diagnostics_tables <- function(combined_summary, combined_diagnostics, loo_tbl, out) {
  tmp <- combined_diagnostics |>
    group_by(model) |>
    summarize(n_div = sum(divergent__)) |>
    ungroup()

  tmp2 <- combined_summary |>
    mutate(check = ifelse(ess_tail < 400, 1, 0)) |>
    group_by(model) |>
    summarize(n_ess = sum(check, na.rm = TRUE)) |>
    ungroup()

  tmp3 <- full_join(tmp, tmp2)

  pattern <- "loo_fit_mcmc_suv_ind_"
  loo_tbl |>
    mutate(model = str_remove(model, pattern)) |>
    full_join(tmp3) |>
    my_write_csv(out)
}

render_diagnostics_tables <- function(diagnostics_tbl, season, abund = FALSE) {
  if(abund) {
    tmp <- diagnostics_tbl |>
        filter(traits != "nlog" & traits != "pc12" & traits != "pc15")
  } else {
    tmp <- diagnostics_tbl |>
        filter(traits %in% c("nlog", "pc12"))
  }

  tmp <- tmp  |>
    filter(season == {{season}}) |>
    arrange(desc(elpd)) |>
    mutate(traits = case_when(
      traits == "nlog" ~ "All traits",
      traits == "pc12" ~ "PC1 and PC2",
      traits == "ab" ~ "Abundance",
      traits == "ba" ~ "Basal area",
      traits == "ab1ba" ~ "Abundance + Basal area",
      traits == "ab2ba" ~ "Abundance $\\times$ Basal area",
    )) |>
    mutate_if(is.numeric, \(x) round(x, digits = 1)) |>
    mutate(elpd = format(elpd, nsmall = 1, trim = TRUE)) |>
    mutate(looic = format(looic, nsmall = 1, trim = TRUE)) |>
    mutate(rain = case_when(
      rain == "norain" ~ "No Rainfall",
      rain == "rain" ~ "Rainfall",
      rain == "intrain" ~ "All predictors $\\times$ Rainfall",
      rain == "intrain2" ~ "ConS $\\times$ Rainfall",
      rain == "intrain3" ~ "ConT $\\times$ Rainfall",
      rain == "intrain4" ~ "(ConS + ConT) $\\times$ Rainfall",
    ))

  if (abund)  {
    tmp <- tmp |>
      dplyr::select(
        `Rainfall model` = rain,
        `Abundance` = traits,
        `ELPD` = elpd,
        `LOOIC` = looic,
        `N\\_ESS` = n_ess,
        `N\\_Div` = n_div)
  } else {
    tmp <- tmp |>
      dplyr::select(
        `Rainfall model` = rain,
        `Traits` = traits,
        `ELPD` = elpd,
        `LOOIC` = looic,
        `N\\_ESS` = n_ess,
        `N\\_Div` = n_div)
  }

  kable_output <- tmp  |>
    kbl(booktabs = TRUE, escape = FALSE, format = "latex", longtable = TRUE) |>
    kable_styling(latex_options = c("striped", "scale_down", "HOLD_position", "repeat_header"),
      full_width = FALSE)

   for (i in 1:nrow(tmp)) {
    if (tmp$`N\\_ESS`[i] > 0 | tmp$`N\\_Div`[i] >= 4) {
      kable_output <- kable_output %>%
        row_spec(i, color = "gray")
     }
   }
   kable_output
}

gen_si_tab <- function(data) {

  pred_name_map <- c(
    "(Intercept)" = "Intercept",
    "logh_s" = "ln Height",
    "scon_s" = "ConS",
    "acon_s_c" = "ConT",
    "shet_s" = "HetS",
    "ahet_s_c" = "HetT",
    "sphy_s" = "PhyS",
    "aphy_s_c" = "PhyT",
    "rain_s" = "Rainfall",
    "scon_s:rain_s" = "ConS $\\times$ Rainfall",
    "acon_s_c:rain_s" = "ConT $\\times$ Rainfall",
    "shet_s:rain_s" = "HetS $\\times$ Rainfall",
    "ahet_s_c:rain_s" = "HetT $\\times$ Rainfall",
    "sphy_s:rain_s" = "PhyS $\\times$ Rainfall",
    "aphy_s_c:rain_s" = "PhyT $\\times$ Rainfall",
    "logh_s:rain_s" = "ln Height $\\times$ Rainfall"
  )

  trait_name_map <- c(
    "intercept" = "Intercept",
    "ldmc" = "LDMC",
    "sdmc" = "SDMC",
    "chl" = "Chl",
    "c" = "C",
    "c13" = "$\\delta \\mathrm{C_{13}}$",
    "log_n" = "ln N",
    "tlp" = "$\\pi_\\mathrm{{tlp}}$",
    "log_la" = "ln LA",
    "log_sla" = "ln SLA",
    "log_lt" = "ln LT",
    "log_ab" = "ln Abundance",
    "dim_1" = "PC1",
    "dim_2" = "PC2"
  )

  data |>
    mutate_if(is.numeric, round, 3) |>
    mutate(
      pred_name = coalesce(pred_name_map[pred_name], pred_name),
      trait_name = coalesce(trait_name_map[trait_name], trait_name),
      q50 = cell_spec(q50, bold = ifelse(q2.5 * q97.5 > 0, TRUE, FALSE)),
      `95\\% CI` = paste0("[", q2.5, ", ", q97.5, "]"),
      `95\\% CI` = cell_spec(`95\\% CI`, bold = ifelse(q2.5 * q97.5 > 0, TRUE, FALSE)),
      Parameter = paste0("$\\gamma_{", str_split_fixed(variable, ",|\\[|\\]", 4)[, 2],
                         ",", str_split_fixed(variable, ",|\\[|\\]", 4)[, 3], "}$")
    ) |>
    rename(
      Median = q50,
      `Lower 2.5\\% CI` = q2.5,
      `Upper 95\\% CI` = q95,
      `Upper 97.5\\% CI` = q97.5,
      `Individual-level predictor` = pred_name,
      `Species-level predictor` = trait_name,
      Rhat = rhat
    ) |>
    dplyr::select(Parameter, `Individual-level predictor`, `Species-level predictor`, Median, `95\\% CI`, Rhat) |>
    kbl(booktabs = TRUE, escape = FALSE, format = "latex", longtable = TRUE) |>
    kable_styling(latex_options = c("striped", "scale_down", "HOLD_position", "repeat_header"),
      full_width = FALSE)
}

pre_glm <- function(seedling_csv) {
  seedling <- read_csv(seedling_csv) |>
    janitor::clean_names()
  cc <- 0.26
  cc2 <- 0.26
  seedling_data <- seedling |>
    mutate(scon_s = scale(scon) |> as.numeric()) |>
    mutate(shet_s = scale(shet) |> as.numeric()) |>
    mutate(sphy_s = scale(sphy) |> as.numeric()) |>
    mutate(acon_s_c = as.numeric(scale(acon^cc))) |>
    mutate(ahet_s_c = as.numeric(scale(ahet^cc))) |>
    mutate(aphy_s_c = as.numeric(scale(aphy^cc2))) |>
    mutate(logh_s = log(h1) |> scale() |> as.numeric()) |>
    mutate(rain_s = as.numeric(scale(rf)))
  fit <- glm(surv ~ (logh_s +
          scon_s + shet_s +
          acon_s_c + ahet_s_c) * season, family = binomial, data = seedling_data)
  fit
}

# Utility function for scaling traits
scale_traits <- function(data) {
  data |>
    mutate(across(c(la, sla, n, lt, ab, ba), ~log(.x), .names = "log_{.col}"),
           latin = data$latin) |>
    select(-la, -sla, -lt, -ab, -ba) |>
    mutate(across(where(is.numeric), ~scale(.x) %>% as.numeric()))
}

# Utility function for scaling seedling data
scale_seedling <- function(data, cc) {
  data |>
    mutate(scon_s = scale(scon) |> as.numeric(),
           shet_s = scale(shet) |> as.numeric(),
           acon_s_c = as.numeric(scale(acon^cc)),
           ahet_s_c = as.numeric(scale(ahet^cc)),
           logh_s = log(h1) |> scale() |> as.numeric(),
           rain_s = as.numeric(scale(rf)))
}

# Utility function to generate model matrix based on rain
generate_model_matrix <- function(data, rain) {
  formula <- switch(rain,
    norain = surv ~ logh_s + scon_s + shet_s + acon_s_c + ahet_s_c,
    intrain = surv ~ (logh_s + scon_s + shet_s + acon_s_c + ahet_s_c) * rain_s,
    intrain2 = surv ~ logh_s + scon_s + shet_s + acon_s_c + ahet_s_c + rain_s +
      scon_s:rain_s,
    intrain3 = surv ~ logh_s + scon_s + shet_s + acon_s_c + ahet_s_c + rain_s +
      acon_s_c:rain_s,
    intrain4 = surv ~ logh_s + scon_s + shet_s + acon_s_c + ahet_s_c + rain_s +
      scon_s:rain_s + acon_s_c:rain_s,
    rain = surv ~ logh_s + scon_s + shet_s + acon_s_c + ahet_s_c + rain_s
  )
  model.matrix(formula, data = data)
}

# Utility function to get predictors based on sp_pred
get_sp_pred_data <- function(data, pca_data, sp_pred) {
  if (sp_pred %in% c("nlog", "n", "ab", "ba", "ab1ba", "ab2ba")) {
    cols <- switch(sp_pred,
               nlog = setdiff(names(data), c("n", "log_ab", "log_ba")),
               n = setdiff(names(data), c("log_n", "log_ab", "log_ba")),
               ab = c("latin", "log_ab"),
               ba = c("latin", "log_ba"),
               ab1ba = c("latin", "log_ab", "log_ba"),
               ab2ba = c("latin", "log_ab", "log_ba", "ab2ba"))
    data |> dplyr::select(all_of(cols))
  } else if (sp_pred %in% c("pc12", "pc15")) {
    range <- as.integer(str_sub(sp_pred, -1, -1))
    range <- min(range, 5)
    cols_to_select <- paste0("dim_", 1:range)
    pca_data |> dplyr::select(latin, all_of(cols_to_select))
  } else {
    data
  }
}

#' @title Create data list for stan
#' @para scaling_within_seasons Scaling within seasons or across seasons (default = FALSE)
generate_stan_data <- function(seedling_df, traits_df, scale_cc, season = "dry", rain = "norain", sp_pred = "nlog") {
  # Read and clean data
  # add station column
  seedling <- seedling_df |>
    mutate(station = str_extract(quadrat, "\\d+"))
  traits <- traits_df

  # Season-based filtering and scalinu
  cc_season <- ifelse(season == "wet", "wet", "dry")
  cc <- scale_cc[[cc_season]]
  season <- ifelse(season == "dry", "dry", "rainy")
  seedling <- seedling |> filter(season == {{season}})

  # Adjust traits and seedlings based on sp_pred and rain
  pca_df <- generate_pca_data(traits_df)
  traits2 <- traits |>
    dplyr::select(-wd, -cn) |>
    scale_traits() |>
    get_sp_pred_data(pca_df, sp_pred)
  seedling <- scale_seedling(seedling, cc)
  Xd <- generate_model_matrix(seedling, rain)

  # sp-level matrix for the model
  Ud <- rbind(
    intercept = rep(1, length(traits2$latin)),
    traits2 |>
      dplyr::select(-latin) |>
      as.matrix() |> t())

  n_sp_d <- length(unique(seedling$latin))
  n_para_d <- ncol(Xd)
  n_plot_d <- length(unique(seedling$quadrat))
  n_station_d <- length(unique(seedling$station))
  n_tag_d <- length(unique(seedling$tag))
  n_census_d <- length(unique(seedling$census))

  list(
    N = nrow(seedling),
    J = n_sp_d,
    K = n_para_d,
    S = n_plot_d,
    T = n_census_d,
    M = n_tag_d,
    P = n_station_d,
    L = nrow(Ud),
    cc = cc,
    suv = seedling$surv,
    plot = seedling$quadrat |>
      as.character() |> as.factor() |> as.integer(),
    station = seedling$station |>
      as.character() |> as.factor() |> as.integer(),
    census = seedling$census |>
      as.character() |> as.factor() |> as.integer(),
    sp = seedling$latin |>
      as.character() |> as.factor() |> as.integer(),
    tag = seedling$tag |>
      as.character() |> as.factor() |> as.integer(),
    x = Xd,
    u = Ud)
}
