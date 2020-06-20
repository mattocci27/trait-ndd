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

model_path <- str_c("./model/", model_name, ".stan")

print(paste("Model ", model_name))
print(paste("Model for ", dry, "season"))
print(paste("Use", ng_data))
print(paste("n_iter =", n_iter))
print(paste("n_warm =", n_warm))
print(paste("n_thin =", n_thin))
print(paste("n_chains =", n_chains))
print(paste("adapt_delta =", a_delta))
print(paste("minimum sp abund =", n_ab))

# Data
seedling_all <- read_csv("./data/seedling_for_drought.csv")
#wp <- read_csv("./data/wp.csv")
tlp <- read_csv("./data/tlp.csv")
trait <- read_csv("./data/BB_SeedlingTrait.csv")


# trait
trait2 <- trait %>%
  #dplyr::select(-Species) %>%
  mutate(logLA = log(LA)) %>%
  mutate(logSLA = log(SLA)) %>%
  mutate(logLT = log(LT)) %>%
  dplyr::select(-LA, -SLA, -LT)

colnames(seedling_all)
colnames(tlp)

tmp10 <- seedling_all %>%
  group_by(sp) %>%
  summarize(n = n()) %>%
  filter(n >= n_ab)

# 10
seedling <- right_join(seedling_all, tmp10, by = "sp")

full_dat <- full_join(seedling, trait, by = c("sp" = "Species"))

if (trait_data == "C13") {
  trait_select <- "C13"

} else {
  trait_select <- "SLA"
}

sp_name <- full_dat %>%
  dplyr::select(sp, {{trait_select}}, surv) %>%
  na.omit() %>%
  select(sp) %>%
  unique %>%
  unlist

names(sp_name) <- NULL

seedling_tlp <- seedling %>%
    filter(sp %in% sp_name)

trait3 <- trait2 %>%
  filter(Species %in% sp_name) %>%
  arrange(Species)

trait4 <- trait3 %>%
  filter(Species %in% sp_name) %>%
  arrange(Species) %>%
  dplyr::select(-Species) %>%
  scale %>%
  as_tibble %>%
  mutate(Species = trait3$Species)

dim(trait3)
dim(trait4)


# Model

writeLines(readLines(model_path))

# Dry or wet season

if (dry == "dry") {
  seedling_tlpd <- seedling_tlp %>%
    filter(season == "dry")
} else {
  seedling_tlpd <- seedling_tlp %>%
    filter(season != "dry")
}

##  Use Detto et al. 2019 Ecology letters -----------------------------

x1 <- seedling_tlpd$acon
x2 <- seedling_tlpd$ahet
y <- seedling_tlpd$surv
lik <- numeric(100)
for (i in 1:100) {
  d1 <- x1^(i/100)
  d2 <- x2^(i/100)
  fm1 <- glm(y ~ scale(d1) + scale(d2), family=binomial)
  lik[i] <- logLik(fm1)
}

#plot(1:100, lik, type = "l")
cc <- which(lik == max(lik)) / 100
print(str_c("use c = ", cc, " as a scaling parameter for distance effect"))

seedling_tlpd <- seedling_tlpd %>%
  mutate(SC_ahet = as.numeric(scale(ahet^cc))) %>%
  mutate(SC_acon = as.numeric(scale(acon^cc)))

# -------------------------------------------------------------------------


if (ng_data == "full") {
  Xd <- cbind(rep(1, nrow(seedling_tlpd)),
              seedling_tlpd[,c("S_scon",
                               "SC_acon",
                               "S_shet",
                               "SC_ahet",
                               "S_log_h1")])
} else if (ng_data == "seedling") {
  Xd <- cbind(rep(1, nrow(seedling_tlpd)),
              seedling_tlpd[,c("S_scon",
                               "S_shet",
                               "S_log_h1")])
} else if (ng_data == "adult") {
  Xd <- cbind(rep(1, nrow(seedling_tlpd)),
              seedling_tlpd[,c("SC_acon",
                               "SC_ahet",
                               "S_log_h1")])
}



colnames(Xd)[1] <- "Int"

trait5 <- trait4 %>%
  filter(Species %in% seedling_tlpd$sp)

dim(Xd)

intercept <- rep(1, length(unique(trait5$Species)))

if (trait_data == "C13") {
  print("Sp-level: 1 + StemD + logSLA + C13")
  Ud <- cbind(
    intercept,
    trait5$StemD,
    trait5$logSLA,
    trait5$C13)
} else if (trait_data == "WP") {
  print("Sp-level: 1 + WP")
  Ud <- cbind(
    intercept,
    trait5$TLP)
} else if (trait_data == "LT") {
  print("Sp-level: 1 + StemD + logSLA + logLT")
  Ud <- cbind(
    intercept,
    trait5$StemD,
    trait5$logSLA,
    trait5$logLT)
} else if (trait_data == "PCA") {
  print("Sp-level: 1 + PC1 + PC2")

  trait_pca <- prcomp(
    trait_pca <- trait5 %>%
      dplyr::select(-Species, -C13),
      scale = TRUE, center = TRUE)

  Ud <- cbind(
    intercept,
    trait_pca$x[,1],
    trait_pca$x[,2])
}



n_sp_d <- length(unique(seedling_tlpd$sp))
n_para_d <- ncol(Xd)
n_plot_d <- length(unique(seedling_tlpd$quadrat))
n_census_d <- length(unique(seedling_tlpd$census))
n_tag_d <- length(unique(seedling_tlpd$tag))

list_dat_d <- list(N = nrow(seedling_tlpd),
                   J = n_sp_d,
                   K = n_para_d,
                   S = n_plot_d,
                   T = n_census_d,
                   M = n_tag_d,
                   L = ncol(Ud),
                   suv = seedling_tlpd$surv,
                   plot = seedling_tlpd$quadrat %>%
                     as.character %>% as.factor %>% as.integer,
                   census = seedling_tlpd$census %>%
                     as.character %>% as.factor %>% as.integer,
                   sp = seedling_tlpd$sp %>%
                     as.character %>% as.factor %>% as.integer,
                   tag = seedling_tlpd$tag %>%
                     as.character %>% as.factor %>% as.integer,
                   x = Xd %>% as.matrix,
                   u = Ud)

print(str_c("n_sp = J =", n_sp_d))
print(str_c("n_para = K = ", n_para_d))
print(str_c("n_plot = S = ", n_plot_d))
print(str_c("n_census = T = ", n_census_d))
print(str_c("n_tag = M = ", n_tag_d))

fit <- stan(file = model_path,
            data = list_dat_d,
            verbose = TRUE,
            iter = n_iter,
            warmup = n_warm,
            thin = n_thin,
            chains =  n_chains,
            refresh = 200,
            control = list(adapt_delta = a_delta, max_treedepth = 20))

print(fit, pars = c("gamma", "sig", "lp__"))

save_name <- str_c("./data/", dry, "_spab_", n_ab, "_", model_name, "_", ng_data, "_", trait_data, ".rda")
print(save_name)

save.image(save_name)

print("MCMC done!!")
