rm(list = ls())

# Deps
library(tidyverse)
library(rstan)
library(rstanarm)
library(bayesplot)
options(mc.cores = parallel::detectCores())

set.seed(123)

argv <- commandArgs(trailingOnly = TRUE)
n_iter <- as.numeric(argv[1])
n_warm <- as.numeric(argv[2])
n_thin <- as.numeric(argv[3])
n_chains <- as.numeric(argv[4])
a_delta <- as.numeric(argv[5])
n_ab <- as.numeric(argv[6])
dry <- as.character(argv[7])

print(paste("Model for ", dry, "season"))
print(paste("n_iter =", n_iter))
print(paste("n_warm =", n_warm))
print(paste("n_thin =", n_thin))
print(paste("n_chains =", n_chains))
print(paste("adapt_delta =", a_delta))
print(paste("minimum sp abund =", n_ab))

# Data
seedling_all <- read_csv("./data/seedling_for_drought.csv")
wp <- read_csv("./data/wp.csv")

colnames(seedling_all)
colnames(wp)

tmp10 <- seedling_all %>% 
  group_by(sp) %>%
  summarize(n = n()) %>%
  filter(n >= n_ab)

# 10
seedling <- right_join(seedling_all, tmp10, by = "sp") 

seedling_wp <- seedling[seedling$sp%in%wp$Species,]
shared_wp <- wp[wp$Species%in%seedling$sp,]

seedling$sp %>% unique %>% length
seedling_wp$sp %>% unique %>% length

nrow(seedling_all)
nrow(seedling)
nrow(seedling_wp)

length(unique(seedling_all$sp))
length(unique(seedling$sp))
length(unique(seedling_wp$sp))

# Model

writeLines(readLines("./model/model_ind.stan"))

# Dry or wet season

if (dry == "dry") {
  seedling_wpd <- seedling_wp %>%
    filter(season == "dry")
} else {
  seedling_wpd <- seedling_wp %>%
    filter(season != "dry")
}

Xd <- cbind(rep(1, nrow(seedling_wpd)),
            seedling_wpd[,c("S_scon",
                            "S_acon",
                            "S_ahet",
                            "S_shet",
                            "S_soilpc1",
                            "S_soilpc2",
                            "S_soilpc3",
                            "S_log_h1")])

colnames(Xd)[1] <- "Int"

shared_wpd2 <- shared_wp %>%
  filter(Species %in% seedling_wpd$sp)

dim(Xd)

Ud <- cbind(rep(1,length(unique(seedling_wpd$sp))),
            shared_wpd2$avg_spwp)

Ud[, 2] <- scale(Ud[,2])

dim(Ud)

n_sp_d <- length(unique(seedling_wpd$sp))
n_para_d <- ncol(Xd)
n_plot_d <- length(unique(seedling_wpd$quadrat))
n_census_d <- length(unique(seedling_wpd$census))
list_dat_d <- list(N = nrow(seedling_wpd),
                   J = n_sp_d,
                   K = n_para_d,
                   S = n_plot_d,
                   T = n_census_d,
                   L = ncol(Ud),
                   suv = seedling_wpd$surv,
                   plot = seedling_wpd$quadrat %>%
                     as.character %>% as.factor %>% as.integer,
                   census = seedling_wpd$census %>%
                     as.character %>% as.factor %>% as.integer,
                   sp = seedling_wpd$sp %>%
                     as.character %>% as.factor %>% as.integer,
                   x = Xd %>% as.matrix,
                   # x = t(X),
                   u = Ud)
str(list_dat_d)

fit_dry <- stan(file = "./model/model_ind.stan",
                data = list_dat_d,
                verbose = TRUE,
                iter = n_iter,
                warmup = n_warm,
                thin = n_thin,
                chains =  n_chains,
                refresh = 200,
                control = list(adapt_delta = a_delta, max_treedepth = 20))

print(fit_dry, pars = c("gamma", "sigma"))

save_name <- str_c("./data/", dry, "_spab_", n_ab, ".rda")

save.image(save_name)

print("MCMC done!!")
