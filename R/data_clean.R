
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

generate_pca_data <- function(trait_csv) {
  targets::tar_load(trait_csv)
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
    janitor::clean_names() |>
    my_write_csv("data/trait_pca.csv")
}

calculate_lik <- function(seedling_csv) {
  seedling <- read_csv(seedling_csv) |> janitor::clean_names()
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

  return(list(lik = lik, lik2 = lik2))
}

trait <- tar_read(trait_csv) |>
  read_csv() |>
  janitor::clean_names()

pca <- prcomp(trait[, 3:14], scale = TRUE)


calc_scale_cc <- function(seedling_csv) {
  results <- calculate_lik(seedling_csv)
  cc <- which(results$lik == max(results$lik)) / 100
  cc2 <- which(results$lik2 == max(results$lik2)) / 100

  c(het = cc, phy = cc2)
}

generate_cc_data <- function(seedling_csv) {
  results <- calculate_lik(seedling_csv)
  tibble(cc = (1:100) / 100, het = results$lik, phy = results$lik2)
}

#targets::tar_load(data_list)
#' @title Create data list for stan
#' @para scaling_within_seasons Scaling within seasons or across seasons (default = FALSE)
generate_stan_data <- function(
  seedling_csv, trait_pca_csv, scale_cc,
  het = c("phy", "het"),
  rain = c("norain", "rain"),
  sp_pred = c("pc12", "pc15", "ab", "ba", "ab1ba")) {

  # targets::tar_load(scale_cc)
  # targets::tar_load(seedling_csv)
  # targets::tar_load(trait_pca_csv)
  # scale_cc <- list(wet = scale_wet, dry = scale_dry)
  # season <- "dry"
  # het <- "phy"
  # sp_pred <- "ab1ba"
  # rain <- "intrain"

  seedling <- read_csv(seedling_csv) |>
    janitor::clean_names()
  traits <- read_csv(trait_pca_csv) |>
    janitor::clean_names() |>
    mutate(log_ab = log(ab)) |>
    mutate(log_ba = log(ba))

  cc <- scale_cc[names(scale_cc) == "het"]
  cc2 <- scale_cc[names(scale_cc) == "phy"]

  if (sp_pred == "pc12") {
    traits2 <- traits |>
      dplyr::select(latin, pc1, pc2)
  } else if (sp_pred == "pc15") {
    traits2 <- traits |>
      dplyr::select(latin, pc1, pc2, pc3, pc4, pc5)
  } else if (sp_pred == "ab") {
    traits2 <- traits |>
      dplyr::select(latin, log_ab)
  } else if (sp_pred == "ba") {
    traits2 <- traits |>
      dplyr::select(latin, log_ba)
  } else if (sp_pred == "ab1ba") {
    traits2 <- traits |>
      dplyr::select(latin, log_ab, log_ba)
  }

  # scale abundance
  if (sp_pred != "pc12" | sp_pred != "pc15") {
    traits2 <- traits2 |>
      summarise_if(is.numeric, \(x) scale(x) |> as.numeric()) |>
      mutate(latin = traits2$latin)
  }

  # tweak sp list for trait and seedling data
  trait_sp <- traits2$latin |> unique()
  seedling_sp <- seedling$latin |> unique()
  sp_c0 <- c(trait_sp, seedling_sp)
  sp_c <- sp_c0[duplicated(sp_c0)] |> unique()

  trait_data <- traits2 |>
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
      Xd <- model.matrix(surv ~ (logh_s +
        scon_s + shet_s +
        acon_s_c + ahet_s_c) * season, data = seedling_data)
    } else {
      Xd <- model.matrix(surv ~ (logh_s +
        scon_s + sphy_s +
        acon_s_c + aphy_s_c) * season, data = seedling_data)
     }
   } else if (rain == "rain") {
    if (het == "het") {
      Xd <- model.matrix(surv ~ (logh_s +
        scon_s + shet_s +
        acon_s_c + ahet_s_c) * season * rain_s, data = seedling_data)
    } else {
      Xd <- model.matrix(surv ~ (logh_s +
        scon_s + sphy_s +
        acon_s_c + aphy_s_c) * season * rain_s, data = seedling_data)
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
