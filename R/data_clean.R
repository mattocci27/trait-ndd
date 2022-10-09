#' @title write_csv for targets
#' @inheritParams readr::write_csv
my_write_csv <- function(x, path, append = FALSE, col_names = !append) {
    write_csv(x, path, append = FALSE, col_names = !append)
    paste(path)
}

seedling_archive <- function(data_list, dry = TRUE) {
  d <- data_list[[1]]
  if (dry) {
    d <- d |>
      dplyr::filter(season == "dry")
  } else {
    d <- d |>
      dplyr::filter(season == "rainy")
  }
  d2 <- d |>
    dplyr::select(tag, quadrat, sp_code, height, census,
    season, survive, cons, cona, heta,hets, rainfall) |>
    mutate(plot = quadrat |> as.character() |> as.factor() |> as.integer()) |>
    mutate(census = census |> as.character() |> as.factor() |> as.integer()) |>
    # mutate(sp_code = sp |> as.character() |> as.factor() |> as.integer()) |>
    mutate(tag = tag |> as.character() |> as.factor() |> as.integer()) |>
    mutate(cons_scaled = scale(cons) |> as.numeric()) |>
    mutate(hets_scaled = scale(hets) |> as.numeric()) |>
    mutate(logh_scaled = log(height) |> scale() |> as.numeric()) |>
    mutate(rain_scaled = as.numeric(scale(rainfall))) |>
    dplyr::select(plot, census, sp = sp_code, tag, cons_scaled, hets_scaled, cona, heta,
      logh_scaled, rain_scaled)

  if (dry) {
    d2 |>
    write_csv("data/seedling_dry_season.csv")
    paste("data/seedling_dry_season.csv")
  } else  {
    d2 |>
    write_csv("data/seedling_rainy_season.csv")
    paste("data/seedling_rainy_season.csv")
  }
}

trait_archive <- function(data_list) {
  d <- data_list[[2]]
  d |>
    dplyr::select(sp, ldmc, wd, sdmc, chl, c13, c_mass, n_mass,
    cn, tlp, log_la, log_sla, log_lt) |>
  write_csv("data/trait.csv")
  paste("data/trait.csv")
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
    mutate(sphy = sphy |> as.numeric())

  if (write_trait) {
    my_write_csv(trait2, "data/traits.csv")
  } else {
    my_write_csv(dataseedling2, "data/seedling.csv")
  }
}

calc_scale_cc <- function(seedling_csv) {
  seedling <- read_csv(seedling_csv) |>
    janitor::clean_names()
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


#targets::tar_load(data_list)
#' @title Create data list for stan
#' @para scaling_within_seasons Scaling within seasons or across seasons (default = FALSE)
generate_stan_data <- function(seedling, traits, scale_cc,
                        model = c("phy_season", "het_season", "phy_rain", "het_rain"),
                        abund = c("abund", "ba", "both")
                        ) {

  # seedling <- read_csv("data/seedling.csv") |>
  #   janitor::clean_names()
  # traits <- read_csv("data/traits.csv") |>
  #   janitor::clean_names()

  seedling <- read_csv(seedling) |>
    janitor::clean_names()
  traits <- read_csv(traits) |>
    janitor::clean_names()

  traits2 <- traits |>
    mutate(log_la = log(la)) |>
    mutate(log_sla = log(sla)) |>
    mutate(log_lt = log(lt)) |>
    mutate(log_ab = log(ab)) |>
    mutate(log_ba = log(ba)) |>
    mutate(log_chl = log(chl)) |>
    mutate(log_c = log(c)) |>
    mutate(log_n = log(n)) |>
    dplyr::select(-la, -sla, -lt, -ab, -ba, -chl, -c, -n) |>
    dplyr::select(latin, ldmc, sdmc, log_la, log_sla, log_chl, log_lt, c13, log_c, log_n, tlp, log_ab, log_ba)

  traits3 <- traits2 |>
    summarise_if(is.numeric, \(x) scale(x) |> as.numeric()) |>
    mutate(latin = traits2$latin)

  if (abund == "abund") {
    traits4 <- traits3 |>
      dplyr::select(-log_ba)
  } else if (abund == "ba") {
    traits4 <- traits3 |>
      dplyr::select(-log_ab)
  } else {
    traits4 <- traits3
  }

  # tweak sp list for trait and seedling data
  trait_sp <- traits4$latin |> unique()
  seedling_sp <- seedling$latin |> unique()
  sp_c0 <- c(trait_sp, seedling_sp)
  sp_c <- sp_c0[duplicated(sp_c0)] |> unique()

  trait_data <- traits4 |>
    filter(latin %in% sp_c)

  seedling_data <- seedling |>
    filter(latin %in% sp_c)

  # sp-level matrix for the model
  Ud <- rbind(
    intercept = rep(1, length(trait_data$latin)),
    trait_data |>
      dplyr::select(-latin) |>
      as.matrix() |> t())

  ##  Use Detto et al. 2019 Ecology letters -----------------------------
  # x1 <- seedling_data$acon
  # x2 <- seedling_data$ahet
  # z1 <- seedling_data$aphy
  # y <- seedling_data$surv

  # n_len <- 100

  # lik <- numeric(n_len)
  # lik2 <- numeric(n_len)
  # for (i in 1:n_len) {
  #   d1 <- x1^(i / n_len)
  #   d2 <- x2^(i / n_len)
  #   d3 <- z1^(i / n_len)
  #   fm1 <- glm(y ~ d1 + d2, family = binomial)
  #   fm2 <- glm(y ~ d1 + d3, family = binomial)
  #   lik[i] <- logLik(fm1)
  #   lik2[i] <- logLik(fm2)
  # }
  # cc <- which(lik == max(lik)) / n_len
  # cc2 <- which(lik2 == max(lik2)) / n_len

  cc <- scale_cc[names(scale_cc) == "het"]
  cc2 <- scale_cc[names(scale_cc) == "phy"]

  seedling_data <- seedling_data |>
    mutate(scon_s = scale(scon) |> as.numeric()) |>
    mutate(shet_s = scale(shet) |> as.numeric()) |>
    mutate(sphy_s = scale(sphy) |> as.numeric()) |>
    mutate(acon_s_c = as.numeric(scale(acon^cc))) |>
    mutate(ahet_s_c = as.numeric(scale(ahet^cc))) |>
    mutate(aphy_s_c = as.numeric(scale(aphy^cc2))) |>
    mutate(logh_s = log(h1) |> scale() |> as.numeric()) |>
    mutate(rain_s = as.numeric(scale(rf)))

  if (model == "het_season") {
    Xd <- model.matrix(surv ~ (logh_s +
      scon_s + shet_s +
      acon_s_c + ahet_s_c) * season, data = seedling_data)
  } else if (model == "phy_season") {
    Xd <- model.matrix(surv ~ (logh_s +
      scon_s + sphy_s +
      acon_s_c + aphy_s_c) * season, data = seedling_data)
  } else if (model == "het_rain") {
    Xd <- model.matrix(surv ~ (logh_s +
      scon_s + shet_s +
      acon_s_c + ahet_s_c) * season * rain_s, data = seedling_data)
  } else if (model == "phy_rain") {
    Xd <- model.matrix(surv ~ (logh_s +
      scon_s + sphy_s +
      acon_s_c + aphy_s_c) * season * rain_s, data = seedling_data)
  }

  n_sp_d <- length(unique(seedling_data$latin))
  n_para_d <- ncol(Xd)
  n_plot_d <- length(unique(seedling_data$quadrat))
  n_census_d <- length(unique(seedling_data$census))
  n_tag_d <- length(unique(seedling_data$tag))

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
       sp = seedling_data$sp |>
         as.character() |> as.factor() |> as.integer(),
       tag = seedling_data$tag |>
         as.character() |> as.factor() |> as.integer(),
       x = Xd,
       u = Ud)
}
