rm(list = ls())

# Deps
library(tidyverse)
library(rstan)
#library(rstanarm)
library(bayesplot)
options(mc.cores = parallel::detectCores())

set.seed(123)

argv <- commandArgs(trailingOnly = TRUE)
model_name <- as.character(argv[1])
n_iter <- as.numeric(argv[2])
n_warm <- as.numeric(argv[3])
n_thin <- as.numeric(argv[4])
n_chains <- as.numeric(argv[5])
a_delta <- as.numeric(argv[6])
n_ab <- as.numeric(argv[7])
dry <- as.character(argv[8])
ng_data <- as.character(argv[9])
trait_data <- as.character(argv[10])
hab <- as.character(argv[11])


#model_name <- "model_ind"
#n_iter <- 10
#n_warm <- 5
#n_thin <- 1
#n_chains <- 1
#a_delta <-  0.8
#n_ab <- 50
#dry <- "dry"
#ng_data <- "full" 
#trait_data <- "Full"
#hab <- "ridge"

model_path <- str_c("./model/", model_name, ".stan")

print(paste("Model ", model_name))
print(paste("Model for ", dry, "season"))
print(paste("Use", ng_data))
print(paste("Habitat =", hab))
print(paste("n_iter =", n_iter))
print(paste("n_warm =", n_warm))
print(paste("n_thin =", n_thin))
print(paste("n_chains =", n_chains))
print(paste("adapt_delta =", a_delta))
print(paste("minimum sp abund =", n_ab))

# Data
seedling_all <- read_csv("./data/seedlingmatrix.csv") |>
  dplyr::select(!starts_with("PC")) |>
  mutate(tmp = str_match(SPcode, "BBSP([0-9]+)")[,2]) |>
  mutate(sp =
    case_when(
      str_length(tmp) == 1 ~ str_c("BBSP00", tmp),
      str_length(tmp) == 2 ~ str_c("BBSP0", tmp),
      str_length(tmp) > 2 ~ str_c("BBSP", tmp)
  ))

hab_dat <- read_csv("./data/habitat150.csv")

trait <- read_csv("./data/trait SPcode.csv") |>
  rename(sp = SPcode) |>
  #dplyr::select(-Cname) |>
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
  #dplyr::select(-Species) |>
  mutate(logLA = log(LA)) |>
  mutate(logSLA = log(SLA)) |>
  mutate(logLT = log(LT)) |>
  dplyr::select(-LA, -SLA, -LT) |>
  arrange(sp) |>
  na.omit()

tmp10 <- seedling_all |>
  group_by(sp) |>
  summarize(n = n()) |>
  filter(n >= n_ab)

# seedling data for abundance >= n_ab
seedling <- right_join(seedling_all, tmp10, by = "sp") |>
  full_join(hab_dat, by = "qua") |>
  mutate(CONS_scaled = scale(CONS) |> as.numeric()) |>
  mutate(HETS_scaled = scale(HETS) |> as.numeric()) |>
  mutate(CONA_scaled = scale(CONA) |> as.numeric()) |>
  mutate(HETA_scaled = scale(HETA) |> as.numeric()) |>
  mutate(logH_scaled = log(height) |> scale() |> as.numeric()) 

# seedling + trait
# data with only traits are available
full_dat <- right_join(seedling, trait2, by = "sp")

#full_dat |>
#  filter(is.na(logSLA)) |>
#  pull(SPcode) |>
#  unique()

# sp list
sp_name <- full_dat |>
  pull(sp) |>
  unique()

# drop species from trait data
trait3 <- trait2 |>
  filter(sp %in% sp_name)

message("No. of species")
message(length(sp_name))

pca_res <- prcomp(
  trait3 |>
  dplyr::select(-sp),
  scale = TRUE, center = TRUE)

trait_pca <- bind_cols(trait3,
          pca_res$x[,1:(ncol(trait3) - 2)] |> 
          as_tibble()) |>
  dplyr::select(sp, starts_with("PC"))

trait4 <- full_join(trait3, trait_pca, by = "sp")

# full trait with scaled values
# still need to adjust sp number according to seedling data
trait5 <- trait4 |>
  dplyr::select(!starts_with("PC")) |>
  dplyr::select(-sp) |>
  scale() |>
  as_tibble() |>
  bind_cols(trait4 |>
  dplyr::select(starts_with("PC"), sp))

# seedling data with traits for analysis
seedling_dat <- seedling |>
  filter(season == {{dry}})

if (hab == "valley" | hab == "ridge" | hab == "slope") {
  seedling_dat <- seedling_dat |>
    filter(habit3 == {{hab}})
}

# different trait sets
# trait data
#- Full
#- except for WD
#- use PC1-3

if (trait_data == "Full") {
  print("Sp-level: 1 + all the traits")
  trait6 <- trait5 |>
    dplyr::select(!starts_with("PC")) |>
    na.omit()
} else if (trait_data == "WD") {
  print("Sp-level: except for WD")
  trait6 <- trait5 |>
    dplyr::select(!starts_with("PC")) |>
    dplyr::select(-WD) |>
    na.omit()
} else if (trait_data == "PCA") {
  print("Sp-level: 1 + PC1 + PC2 + PC3")
  trait6 <- trait5 |>
    dplyr::select(sp, PC1, PC2, PC3) |>
    na.omit()
} else if (trait_data == "SLA_WD_LT") {
  print("Sp-level: 1 + SLA + WD + LT")
  trait6 <- trait5 |>
    dplyr::select(sp, logSLA, logLT, WD) |>
    na.omit()
}

# tweak sp list for trait and seedling data
trait_sp <- trait6$sp |> unique()
seedling_sp <- seedling_dat$sp |> unique()
sp_c0 <- c(trait_sp, seedling_sp)
sp_c <- sp_c0[duplicated(sp_c0)] |> unique()

trait_dat <- trait6 |>
  filter(sp %in% sp_c)

seedling_dat2 <- seedling_dat |>
  filter(sp %in% sp_c)

str_c("sp number in seedling data: ",
  seedling_dat2$sp |> unique() |> length()) |> print()
str_c("sp number in trait data: ",
  trait_dat$sp |> unique() |> length()) |> print()

# sp-level matrix for the model

Ud <- cbind(
  intercept = rep(1, length(trait_dat$sp)),
  trait_dat |>
    dplyr::select(-sp) |>
    as.matrix())

# Model

writeLines(readLines(model_path))

##  Use Detto et al. 2019 Ecology letters -----------------------------

x1 <- seedling_dat2$CONA
x2 <- seedling_dat2$HETA
y <- seedling_dat2$survive
lik <- numeric(100)
for (i in 1:100) {
  d1 <- x1^(i/100)
  d2 <- x2^(i/100)
  #fm1 <- glm(y ~ scale(d1) + scale(d2), family = binomial)
  fm1 <- glm(y ~ d1 + d2, family = binomial)
  lik[i] <- logLik(fm1)
}

#plot(seq(0.01, 1, length = 100), lik, type = "l")

cc <- which(lik == max(lik)) / 100
print(str_c("use c = ", cc, " as a scaling parameter for the distance effect"))

seedling_dat2 <- seedling_dat2 |>
  mutate(HETA_scaled_c = as.numeric(scale(HETA^cc))) |>
  mutate(CONA_scaled_c = as.numeric(scale(CONA^cc))) |>
  mutate(rain_scaled = as.numeric(scale(Rainfall))) |>
  mutate(HETA_rain = HETA_scaled_c * rain_scaled) |>
  mutate(HETS_rain = HETS_scaled * rain_scaled) |>
  mutate(CONA_rain = CONA_scaled_c * rain_scaled) |>
  mutate(CONS_rain = CONS_scaled * rain_scaled) 

# -------------------------------------------------------------------------


if (ng_data == "full") {
  #Xd <- cbind(rep(1, nrow(seedling_dat2)),
  #            seedling_dat2[,c("CONS_scaled",
  #                             "CONA_scaled_c",
  #                             "HETS_scaled",
  #                             "HETA_scaled_c",
  #                             "logH_scaled")])
  Xd <- cbind(rep(1, nrow(seedling_dat2)),
              seedling_dat2[,c("CONS_scaled",
                               "CONA_scaled_c",
                               "HETS_scaled",
                               "HETA_scaled_c",
                               "rain_scaled",
                               "CONS_rain",
                               "CONA_rain",
                               "HETS_rain",
                               "HETA_rain",
                               "logH_scaled")])
} else if (ng_data == "seedling") {
  Xd <- cbind(rep(1, nrow(seedling_dat2)),
              seedling_dat2[,c("CONS_scaled",
                               "HETS_scaled",
                               "rain_scaled",
                               "CONS_rain",
                               "HETS_rain",
                               "logH_scaled")])
} else if (ng_data == "adult") {
  Xd <- cbind(rep(1, nrow(seedling_dat2)),
              seedling_dat2[,c(
                               "CONA_scaled_c",
                               "HETA_scaled_c",
                               "rain_scaled",
                               "CONA_rain",
                               "HETA_rain",
                               "logH_scaled")])
}

colnames(Xd)[1] <- "Int"

n_sp_d <- length(unique(seedling_dat2$sp))
n_para_d <- ncol(Xd)
n_plot_d <- length(unique(seedling_dat2$quadrat))
n_census_d <- length(unique(seedling_dat2$census))
n_tag_d <- length(unique(seedling_dat2$tag))

list_dat_d <- list(N = nrow(seedling_dat2),
                   J = n_sp_d,
                   K = n_para_d,
                   S = n_plot_d,
                   T = n_census_d,
                   M = n_tag_d,
                   L = ncol(Ud),
                   suv = seedling_dat2$survive,
                   plot = seedling_dat2$quadrat |>
                     as.character() |> as.factor() |> as.integer(),
                   census = seedling_dat2$census |>
                     as.character() |> as.factor() |> as.integer(),
                   sp = seedling_dat2$sp |>
                     as.character() |> as.factor() |> as.integer(),
                   tag = seedling_dat2$tag |>
                     as.character() |> as.factor() |> as.integer(),
                   x = Xd |> as.matrix(),
                   u = Ud)

print(str_c("n_sp = J =", n_sp_d))
print(str_c("n_para = K = ", n_para_d))
print(str_c("n_plot = S = ", n_plot_d))
print(str_c("n_census = T = ", n_census_d))
print(str_c("n_tag = M = ", n_tag_d))

model <- stan_model(model_path) 

fit <- sampling(model,
            data = list_dat_d,
            verbose = TRUE,
            iter = n_iter,
            warmup = n_warm,
            thin = n_thin,
            chains =  n_chains,
            refresh = 200,
            control = list(adapt_delta = a_delta, max_treedepth = 20))

print(fit, pars = c("gamma", "sig", "lp__"))

save_name <- str_c("./rda/", dry, "_spab_", n_ab, "_", model_name, "_", ng_data, "_", trait_data, "_", hab, ".rda")
print(save_name)

save.image(save_name)

print("MCMC done!!")
