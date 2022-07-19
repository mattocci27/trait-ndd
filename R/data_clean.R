
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

  habitat <- read_csv(habitat_csv) |>
    janitor::clean_names()
# seedling data for abundance >= n_ab
  seedling <- right_join(seedling, tmp10, by = "sp") |>
    full_join(habitat, by = "qua") |>
    mutate(cons_scaled = scale(cons) |> as.numeric()) |>
    mutate(hets_scaled = scale(hets) |> as.numeric()) |>
    mutate(cona_scaled = scale(cona) |> as.numeric()) |>
    mutate(heta_scaled = scale(heta) |> as.numeric()) |>
    mutate(logh_scaled = log(height) |> scale() |> as.numeric())

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
        dplyr::select(qua:logh_scaled),
        trait = trait_re)

}

#targets::tar_load(data_list)
gen_stan_dat <- function(data_list, season = "dry", habitat = "all",
                        inter = TRUE,
                        trait_set = c("cn", "wd", "pca", "sla_wd_lt")) {
  seedling <- data_list$seedling |>
    filter(season == {{season}})
  if (habitat != "all") {
  seedling<- seedling_dat |>
    filter(habit3 == {{habitat}})
  }

  if (trait_set == "cn") {
    trait <- data_list$trait |>
      dplyr::select(!starts_with("pc")) |>
      dplyr::select(-wd) |>
      dplyr::select(-cn) |>
      na.omit()
  } else if (trait_set == "wd") {
    trait <- data_list$trait |>
      dplyr::select(!starts_with("pc")) |>
      dplyr::select(-wd) |>
      na.omit()
  } else if (trait_set == "pca") {
    trait <- data_list$trait |>
      dplyr::select(sp, pc1, pc2, pc3) |>
      na.omit()
  } else if (trait_set == "sla_wd_lt") {
    trait <- data_list$trait |>
      dplyr::select(sp, log_sla, log_lt, wd) |>
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

  seedling_dat <- seedling_dat |>
    mutate(heta_scaled_c = as.numeric(scale(heta^cc))) |>
    mutate(cona_scaled_c = as.numeric(scale(cona^cc))) |>
    mutate(rain_scaled = as.numeric(scale(rainfall))) |>
    mutate(heta_rain = heta_scaled_c * rain_scaled) |>
    mutate(hets_rain = hets_scaled * rain_scaled) |>
    mutate(cona_rain = cona_scaled_c * rain_scaled) |>
    mutate(cons_rain = cons_scaled * rain_scaled)


  if (inter) {
    Xd <- cbind(rep(1, nrow(seedling_dat)),
                seedling_dat[,c("cons_scaled",
                                 "cona_scaled_c",
                                 "hets_scaled",
                                 "heta_scaled_c",
                                 "rain_scaled",
                                 "cons_rain",
                                 "cona_rain",
                                 "hets_rain",
                                 "heta_rain",
                                 "logh_scaled")])
  } else  {
    Xd <- cbind(rep(1, nrow(seedling_dat)),
                seedling_dat[,c("cons_scaled",
                                 "cona_scaled_c",
                                 "hets_scaled",
                                 "heta_scaled_c",
                                 "rain_scaled",
                                 "logh_scaled")])
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
