
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
  cc <- which(lik == max(lik)) / n_len
  cc2 <- which(lik2 == max(lik2)) / n_len
  c(het = cc, phy = cc2)
}

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

generate_pca_data <- function(trait_csv) {
  trait <- read_csv(trait_csv) |>
    janitor::clean_names() |>
    mutate(la = log(la)) |>
    mutate(lt = log(lt)) |>
    mutate(n = log(n)) |>
    mutate(sla = log(sla)) |>
    rename(log_la = la) |>
    rename(log_lt = lt) |>
    rename(log_n = n) |>
    rename(log_sla = sla)
  pca <- prcomp(trait[, 3:14], scale = TRUE)

  bind_cols(trait, pca$x[, 1:5]) |>
    janitor::clean_names() #|>
    # my_write_csv("data/trait_pca.csv")
}

#targets::tar_load(data_list)
#' @title Create data list for stan
#' @para scaling_within_seasons Scaling within seasons or across seasons (default = FALSE)
generate_stan_data <- function(
  seedling, traits, scale_cc,
  season = c("dry", "wet"),
  het = c("phy", "het"),
  rain = c("norain", "rain", "intrain", "intrain2", "intrain3", "intrain4"),
  sp_pred = c("nlog", "n", "ab", "ba", "ab1ba", "ab2ba", "pc12", "pc15")) {

  # targets::tar_load(scale_wet)
  # targets::tar_load(scale_dry)
  # targets::tar_load(seedling)
  # targets::tar_load(traits)
  # scale_cc <- list(wet = scale_wet, dry = scale_dry)
  # season <- "dry"
  # het <- "phy"
  # sp_pred <- "ab1ba"
  # rain <- "intrain"
  pca_df <- generate_pca_data(traits)

  seedling <- read_csv(seedling) |>
    janitor::clean_names()
  traits <- read_csv(traits) |>
    janitor::clean_names()

  if (season == "wet") {
    seedling <- seedling |>
      filter(season == "rainy")
      cc <- scale_cc$wet[names(scale_cc$wet) == "het"]
      cc2 <- scale_cc$wet[names(scale_cc$wet) == "phy"]
  } else {
    seedling <- seedling |>
      filter(season == "dry")
      cc <- scale_cc$dry[names(scale_cc$dry) == "het"]
      cc2 <- scale_cc$dry[names(scale_cc$dry) == "phy"]
  }

  traits2 <- traits |>
    mutate(log_la = log(la)) |>
    mutate(log_sla = log(sla)) |>
    mutate(log_lt = log(lt)) |>
    mutate(log_ab = log(ab)) |>
    mutate(log_ba = log(ba)) |>
    # mutate(log_chl = log(chl)) |>
    # mutate(log_c = log(c)) |>
    mutate(log_n = log(n)) |>
    dplyr::select(-la, -sla, -lt, -ab, -ba) |>
    dplyr::select(latin, ldmc, sdmc, log_la, log_sla, log_lt, c13, log_n, n, tlp, log_ab, log_ba)

  if (sp_pred == "nlog") {
    traits2 <- traits2 |>
      dplyr::select(-n, -log_ab, -log_ba)
  } else if (sp_pred == "n") {
    traits2 <- traits2 |>
      dplyr::select(-log_n, -log_ab, -log_ba)
  } else if (sp_pred == "ab") {
    traits2 <- traits2 |>
      dplyr::select(latin, log_ab)
  } else if (sp_pred == "ba") {
    traits2 <- traits2 |>
      dplyr::select(latin, log_ba)
  } else if (sp_pred == "ab1ba") {
    traits2 <- traits2 |>
      dplyr::select(latin, log_ab, log_ba)
  } else if (sp_pred == "ab2ba") {
    traits2 <- traits2 |>
      dplyr::select(latin, log_ab, log_ba) |>
      mutate(ab2ba = log_ab * log_ba)
  } else if (sp_pred == "pc12") {
    traits2 <- pca_df |>
      dplyr::select(latin, pc1:pc2)
  } else if (sp_pred == "pc15") {
    traits2 <- pca_df |>
      dplyr::select(latin, pc1:pc5)
  }

  # dry and wet season have the same sp number
  # so we can scale the trait data here
  traits3 <- traits2 |>
    summarise_if(is.numeric, \(x) scale(x) |> as.numeric()) |>
    mutate(latin = traits2$latin)

  # tweak sp list for trait and seedling data
  trait_sp <- traits3$latin |> unique()
  seedling_sp <- seedling$latin |> unique()
  sp_c0 <- c(trait_sp, seedling_sp)
  sp_c <- sp_c0[duplicated(sp_c0)] |> unique()

  trait_data <- traits3 |>
    filter(latin %in% sp_c)

  seedling_data <- seedling |>
    filter(latin %in% sp_c)

  # sp-level matrix for the model
  Ud <- rbind(
    intercept = rep(1, length(trait_data$latin)),
    trait_data |>
      dplyr::select(-latin) |>
      as.matrix() |> t())

  seedling_data <- seedling_data |>
    mutate(scon_s = scale(scon) |> as.numeric()) |>
    mutate(shet_s = scale(shet) |> as.numeric()) |>
    mutate(sphy_s = scale(sphy) |> as.numeric()) |>
    mutate(acon_s_c = as.numeric(scale(acon^cc))) |>
    mutate(ahet_s_c = as.numeric(scale(ahet^cc))) |>
    mutate(aphy_s_c = as.numeric(scale(aphy^cc2))) |>
    mutate(logh_s = log(h1) |> scale() |> as.numeric()) |>
    mutate(rain_s = as.numeric(scale(rf)))

  if (rain == "norain") {
    if (het == "het") {
      Xd <- model.matrix(surv ~ logh_s +
        scon_s + shet_s +
        acon_s_c + ahet_s_c, data = seedling_data)
    } else {
      Xd <- model.matrix(surv ~ logh_s +
        scon_s + sphy_s +
        acon_s_c + aphy_s_c, data = seedling_data)
    }
   } else if (rain == "intrain") {
    if (het == "het") {
      Xd <- model.matrix(surv ~ (logh_s +
        scon_s + shet_s +
        acon_s_c + ahet_s_c) * rain_s, data = seedling_data)
    } else {
      Xd <- model.matrix(surv ~ (logh_s +
        scon_s + sphy_s +
        acon_s_c + aphy_s_c) * rain_s, data = seedling_data)
    }
   } else if (rain == "intrain2") {
    if (het == "het") {
      Xd <- model.matrix(surv ~ logh_s +
        scon_s + shet_s +
        acon_s_c + ahet_s_c + rain_s + scon_s:rain_s, data = seedling_data)
    } else {
      Xd <- model.matrix(surv ~ logh_s +
        scon_s + sphy_s +
        acon_s_c + aphy_s_c + rain_s + scon_s:rain_s, data = seedling_data)
    }
   } else if (rain == "intrain3") {
    if (het == "het") {
      Xd <- model.matrix(surv ~ logh_s +
        scon_s + shet_s +
        acon_s_c + ahet_s_c + rain_s + acon_s_c:rain_s, data = seedling_data)
    } else {
      Xd <- model.matrix(surv ~ logh_s +
        scon_s + sphy_s +
        acon_s_c + aphy_s_c + rain_s + acon_s_c:rain_s, data = seedling_data)
    }
   } else if (rain == "intrain4") {
    if (het == "het") {
      Xd <- model.matrix(surv ~ logh_s +
        scon_s + shet_s +
        acon_s_c + ahet_s_c + rain_s + scon_s:rain_s + acon_s_c:rain_s, data = seedling_data)
    } else {
      Xd <- model.matrix(surv ~ logh_s +
        scon_s + sphy_s +
        acon_s_c + aphy_s_c + rain_s + scon_s:rain_s + acon_s_c:rain_s, data = seedling_data)
    }
   } else if (rain == "rain") {
    if (het == "het") {
      Xd <- model.matrix(surv ~ logh_s +
        scon_s + shet_s +
        acon_s_c + ahet_s_c + rain_s, data = seedling_data)
    } else {
      Xd <- model.matrix(surv ~ logh_s +
        scon_s + sphy_s +
        acon_s_c + aphy_s_c + rain_s, data = seedling_data)
    }
   }

  n_sp_d <- length(unique(seedling_data$latin))
  n_para_d <- ncol(Xd)
  n_plot_d <- length(unique(seedling_data$quadrat))
  n_tag_d <- length(unique(seedling_data$tag))

  n_census_d <- length(unique(seedling_data$census))

  if (het == "phy") cc <- cc2

  list(N = nrow(seedling_data),
       J = n_sp_d,
       K = n_para_d,
       S = n_plot_d,
       T = n_census_d,
       M = n_tag_d,
       L = nrow(Ud),
       cc = cc,
       suv = seedling_data$surv,
       plot = seedling_data$quadrat |>
         as.character() |> as.factor() |> as.integer(),
       census = seedling_data$census |>
         as.character() |> as.factor() |> as.integer(),
       sp = seedling_data$latin |>
         as.character() |> as.factor() |> as.integer(),
       tag = seedling_data$tag |>
         as.character() |> as.factor() |> as.integer(),
       x = Xd,
       u = Ud)
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

  pattern <- "loo_fit_mcmc_logistic_simple_stan_data_"
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
    mutate(phy = ifelse(phy == "phy", "Phylogenetic", "Non-phylogenetic")) |>
    mutate(rain = case_when(
      rain == "norain" ~ "No rain",
      rain == "rain" ~ "Rain without interactions",
      rain == "intrain" ~ "Rain with all the interactions",
      rain == "intrain2" ~ "Rain with an interaction of cons",
      rain == "intrain3" ~ "Rain with an interaction of cona",
      rain == "intrain4" ~ "Rain with an interaction of cons and cona",
    )) |>
    mutate(
      across(1:10,
      \(x) cell_spec(x, color = ifelse(n_ess > 0 | n_div >= 4, "gray", "black"))
    ))

  if (abund)  {
    tmp <- tmp |>
      dplyr::select(`Seedling densities` = phy,
        `Rainfall` = rain,
        `Abundance` = traits,
        `ELPD` = elpd,
        `LOOIC` = looic,
        `N\\_ESS` = n_ess,
        `N\\_Div` = n_div)
  } else {
    tmp <- tmp |>
      dplyr::select(`Seedling densities` = phy,
        `Rainfall` = rain,
        `Traits` = traits,
        `ELPD` = elpd,
        `LOOIC` = looic,
        `N\\_ESS` = n_ess,
        `N\\_Div` = n_div)
  }

  tmp  |>
    kbl(booktabs = TRUE, escape = FALSE, format = "latex", longtable = FALSE) |>
    kable_styling(latex_options = c("striped", "scale_down", "HOLD_position", "repeat_header"), full_width = FALSE) #k|>
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
    "c13" = "$\\delta \\mathrm{C_{13}}$",
    "log_n" = "ln N",
    "tlp" = "$\\pi_\\mathrm{{tlp}}$",
    "log_la" = "ln LA",
    "log_sla" = "ln SLA",
    "log_lt" = "ln LT",
    "log_ab" = "ln Abundance",
    "pc1" = "PC1",
    "pc2" = "PC2"
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
    dplyr::select(Parameter, Median, `95\\% CI`, `Individual-level predictor`, `Species-level predictor`, Rhat) |>
    kbl(booktabs = TRUE, escape = FALSE, format = "latex", longtable = TRUE) |>
    kable_styling(latex_options = c("striped", "scale_down", "HOLD_position", "repeat_header"))
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
