
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

#targets::tar_load(data_list)
#' @title Create data list for stan
#' @para scaling_within_seasons Scaling within seasons or across seasons (default = FALSE)
generate_stan_data <- function(
  seedling, traits, scale_cc,
  season = c("dry", "wet"),
  het = c("phy", "het"),
  rain = c("norain", "rain", "intrain", "intrain2", "intrain3", "intrain4"),
  sp_pred = c("nlog", "n", "ab", "ba", "ab1ba", "ab2ba")) {

  # targets::tar_load(scale_wet)
  # targets::tar_load(scale_dry)
  # targets::tar_load(seedling)
  # targets::tar_load(traits)
  # scale_cc <- list(wet = scale_wet, dry = scale_dry)
  # season <- "dry"
  # het <- "phy"
  # sp_pred <- "ab1ba"
  # rain <- "intrain"

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
