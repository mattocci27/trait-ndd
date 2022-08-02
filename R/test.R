library(tidyverse)

targets::tar_load(dry_each_int_s)
targets::tar_load(wet_each_int_s)
wet_each_int_s |> str()
dry_each_int_s |> str()


  seedling <- data_list$seedling

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
  dry_tmp <- seedling_dat |> filter(season == "dry")
  x1 <- dry_tmp$cona
  x2 <- dry_tmp$heta
  y <- dry_tmp$survive
  lik <- numeric(100)
  tic()
  for (i in 1:100) {
    d1 <- x1^(i/100)
    d2 <- x2^(i/100)
    #fm1 <- glm(y ~ scale(d1) + scale(d2), family = binomial)
    fm1 <- glm(y ~ d1 + d2, family = binomial)
    lik[i] <- logLik(fm1)
  }
  dry_cc <- which(lik == max(lik)) / 100

  ##  Use Detto et al. 2019 Ecology letters -----------------------------
  wet_tmp <- seedling_dat |> filter(season == "rainy")
  x1 <- wet_tmp$cona
  x2 <- wet_tmp$heta
  y <- wet_tmp$survive
  lik <- numeric(100)
  tic()
  for (i in 1:100) {
    d1 <- x1^(i/100)
    d2 <- x2^(i/100)
    #fm1 <- glm(y ~ scale(d1) + scale(d2), family = binomial)
    fm1 <- glm(y ~ d1 + d2, family = binomial)
    lik[i] <- logLik(fm1)
  }
  wet_cc <- which(lik == max(lik)) / 100
