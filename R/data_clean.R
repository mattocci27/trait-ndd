#' @title write_csv for targets
#' @inheritParams readr::write_csv
my_write_csv <- function(x, path, append = FALSE, col_names = !append) {
    write_csv(x, path, append = FALSE, col_names = !append)
    paste(path)
}


#' @title Just cleaning seedling csv
seedling_clean <- function(seedling_csv) {
  seedling_clean <- read_csv(seedling_csv) |>
    janitor::clean_names()
  write_csv(seedling_clean, "data/seedling_clean.csv")
  paste("data/seedling_clean.csv")
}

#' @title Just cleaning trait csv
trait_clean <- function(trait_csv) {
  trait_clean <- read_csv(trait_csv) |>
    janitor::clean_names() |>
    rename(sp_code = s_pcode)
  write_csv(trait_clean, "data/trait_clean.csv")
  paste("data/trait_clean.csv")
}

gen_seedling <- function(seedling_csv, trait_csv, habitat_csv, n_ab = 50) {
  # seedling_csv <- "data/seedlingmatrix.csv"
  # habitat_csv <- "data/habitat150.csv"
  seedling <- read_csv(seedling_csv) |>
    janitor::clean_names() |>
    rename(sp_code = s_pcode) |>
    mutate(tmp = str_match(sp_code, "BBSP([0-9]+)")[,2]) |>
    mutate(sp =
      case_when(
        str_length(tmp) == 1 ~ str_c("BBSP00", tmp),
        str_length(tmp) == 2 ~ str_c("BBSP0", tmp),
        str_length(tmp) > 2 ~ str_c("BBSP", tmp)
    ))

  trait <- read_csv(trait_csv) |>
      janitor::clean_names() |>
      rename(sp = s_pcode) |>
      rename(c_mass = c) |>
      rename(n_mass = n) |>
    # remove species starting with Sn
      filter(str_detect(sp, "^BBSP")) |>
      mutate(tmp = str_match(sp, "BBSP([0-9]+)")[,2]) |>
      mutate(sp =
        case_when(
          str_length(tmp) == 1 ~ str_c("BBSP00", tmp),
          str_length(tmp) == 2 ~ str_c("BBSP0", tmp),
          str_length(tmp) > 2 ~ str_c("BBSP", tmp)
      )) |>
      dplyr::select(-tmp)

  # trait
  trait2 <- trait |>
    mutate(log_la = log(la)) |>
    mutate(log_sla = log(sla)) |>
    mutate(log_lt = log(lt)) |>
    dplyr::select(-la, -sla, -lt) |>
    arrange(sp) |>
    na.omit()

  tmp10 <- seedling |>
    group_by(sp) |>
    summarize(n = n()) |>
    filter(n >= n_ab)

# Eventually we didn't use habitat
  habitat <- read_csv(habitat_csv) |>
    janitor::clean_names()
# seedling data for abundance >= n_ab
  seedling <- right_join(seedling, tmp10, by = "sp") |>
    full_join(habitat, by = "qua")

# seedling + trait
# data with only traits are available
  full_dat <- right_join(seedling, trait2, by = "sp")

# sp list
  sp_name <- full_dat |>
    pull(sp) |>
    unique()

# drop species from trait data
  trait_clean <- trait2 |>
    filter(sp %in% sp_name)

  pca_res <- prcomp(
    trait_clean |>
    dplyr::select(-sp),
    scale = TRUE, center = TRUE)

  trait_pca <- bind_cols(trait_clean,
            pca_res$x[,1:(ncol(trait_clean) - 2)] |>
            as_tibble()) |>
    janitor::clean_names() |>
    dplyr::select(sp, starts_with("pc"))

# full trait with scaled values
# still need to adjust sp number according to seedling data
  trait_tmp <- full_join(trait_clean, trait_pca, by = "sp")

  trait_re <- trait_tmp |>
    dplyr::select(!starts_with("pc")) |>
    dplyr::select(-sp) |>
    scale() |>
    as_tibble() |>
    bind_cols(trait_tmp |>
    dplyr::select(starts_with("pc"), sp))

  list(seedling = full_dat |>
        dplyr::select(qua:habit5),
        trait = trait_re)

}

#targets::tar_load(data_list)
#' @title Create data list for stan
#' @para scaling_within_seasons Scaling within seasons or across seasons (default = FALSE)
gen_stan_dat <- function(data_list,
                        season = "dry",
                        inter = TRUE,
                        trait_set = c("each", "pca"),
                        one_inter = FALSE,
                        scaling_within_seasons = FALSE
                        ) {

  # targets::tar_load(data_list)
  seedling <- data_list$seedling
  if (scaling_within_seasons) {
    seedling <- data_list$seedling |>
      filter(season == {{season}})
  }

  if (trait_set == "each") {
    trait <- data_list$trait |>
      dplyr::select(!starts_with("pc")) |>
      dplyr::select(-wd) |>
      dplyr::select(-cn) |>
      na.omit()
  } else if (trait_set == "pca") {
    trait <- data_list$trait |>
      dplyr::select(sp, pc1, pc2, pc3) |>
      na.omit()
  } else if (trait_set == "cn") {
    trait <- data_list$trait |>
      dplyr::select(!starts_with("pc")) |>
      dplyr::select(-wd) |>
      dplyr::select(-c_mass) |>
      dplyr::select(-n_mass) |>
      na.omit()
  }

  # tweak sp list for trait and seedling data
  trait_sp <- trait$sp |> unique()
  seedling_sp <- seedling$sp |> unique()
  sp_c0 <- c(trait_sp, seedling_sp)
  sp_c <- sp_c0[duplicated(sp_c0)] |> unique()

  trait_dat <- trait |>
    filter(sp %in% sp_c)

  seedling_dat <- seedling |>
    filter(sp %in% sp_c)

  # sp-level matrix for the model
  Ud <- rbind(
    intercept = rep(1, length(trait_dat$sp)),
    trait_dat |>
      dplyr::select(-sp) |>
      as.matrix() |> t())

  ##  Use Detto et al. 2019 Ecology letters -----------------------------
  x1 <- seedling_dat$cona
  x2 <- seedling_dat$heta
  y <- seedling_dat$survive
  lik <- numeric(100)
  for (i in 1:100) {
    d1 <- x1^(i/100)
    d2 <- x2^(i/100)
    #fm1 <- glm(y ~ scale(d1) + scale(d2), family = binomial)
    fm1 <- glm(y ~ d1 + d2, family = binomial)
    lik[i] <- logLik(fm1)
  }

  cc <- which(lik == max(lik)) / 100

  # Scaling will be done across seasons if (!scale_within_seasons)
  # Scaling will be done within seasons if (scale_within_seasons)
  seedling_dat <- seedling_dat |>
    mutate(cons_scaled = scale(cons) |> as.numeric()) |>
    mutate(hets_scaled = scale(hets) |> as.numeric()) |>
    mutate(logh_scaled = log(height) |> scale() |> as.numeric()) |>
    mutate(heta_scaled_c = as.numeric(scale(heta^cc))) |>
    mutate(cona_scaled_c = as.numeric(scale(cona^cc))) |>
    mutate(rain_scaled = as.numeric(scale(rainfall))) |>
    mutate(heta_rain = heta_scaled_c * rain_scaled) |>
    mutate(hets_rain = hets_scaled * rain_scaled) |>
    mutate(acon_rain = acon_scaled_c * rain_scaled) |>
    mutate(cons_rain = cons_scaled * rain_scaled) |>
    filter(season == {{season}})

  if (inter) {
    Xd <- cbind(rep(1, nrow(seedling_dat)),
                seedling_dat[,c( "logh_scaled",
                                 "cons_scaled",
                                 "cona_scaled_c",
                                 "hets_scaled",
                                 "heta_scaled_c",
                                 "rain_scaled",
                                 "cons_rain",
                                 "cona_rain",
                                 "hets_rain",
                                 "heta_rain"
                                 )])
  } else  {
    Xd <- cbind(rep(1, nrow(seedling_dat)),
                seedling_dat[,c( "logh_scaled",
                                 "cons_scaled",
                                 "acon_scaled_c",
                                 "hets_scaled",
                                 "heta_scaled_c",
                                 "rain_scaled")])
  }

  if (one_inter) {
    Xd <- cbind(rep(1, nrow(seedling_dat)),
                seedling_dat[,c( "logh_scaled",
                                 "cons_scaled",
                                 "acon_scaled_c",
                                 "hets_scaled",
                                 "heta_scaled_c",
                                 "rain_scaled",
                                 "acon_rain"
                                 )])
  }

  colnames(Xd)[1] <- "int"

  n_sp_d <- length(unique(seedling_dat$sp))
  n_para_d <- ncol(Xd)
  n_plot_d <- length(unique(seedling_dat$quadrat))
  n_census_d <- length(unique(seedling_dat$census))
  n_tag_d <- length(unique(seedling_dat$tag))

  list(N = nrow(seedling_dat),
       J = n_sp_d,
       K = n_para_d,
       S = n_plot_d,
       T = n_census_d,
       M = n_tag_d,
       L = nrow(Ud),
       cc = cc,
       suv = seedling_dat$survive,
       plot = seedling_dat$quadrat |>
         as.character() |> as.factor() |> as.integer(),
       census = seedling_dat$census |>
         as.character() |> as.factor() |> as.integer(),
       sp = seedling_dat$sp |>
         as.character() |> as.factor() |> as.integer(),
       tag = seedling_dat$tag |>
         as.character() |> as.factor() |> as.integer(),
       x = Xd |> as.matrix(),
       u = Ud)
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


  # traits2 %>%
  #  pivot_longer(ldmc:log_ba) %>%
  #  ggplot(., aes(x = value)) +
  #  geom_histogram() +
  #  facet_wrap(~ name, scale = "free")

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
    mutate(logh_s = log(h1) |> scale() |> as.numeric()) |>
    mutate(ahet_s_c = as.numeric(scale(ahet^cc))) |>
    mutate(acon_s_c = as.numeric(scale(acon^cc))) |>
    mutate(aphy_s_c = as.numeric(scale(aphy^cc2))) |>
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
       suv = seedling_data$survive,
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
