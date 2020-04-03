rm(list = ls())

# Deps
library(tdyverse)
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

model_path <- str_c("./model/", model_name, ".stan")

print(paste("Model ", model_name))
print(paste("Model for ", dry, "season"))
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

colnames(seedling_all)
colnames(tlp)

tmp10 <- seedling_all %>% 
  group_by(sp) %>%
  summarize(n = n()) %>%
  filter(n >= n_ab)

# 10
seedling <- right_join(seedling_all, tmp10, by = "sp") 

seedling_tlp <- seedling[seedling$sp%in%tlp$Species,]
shared_tlp <- tlp[tlp$Species%in%seedling$sp,]

seedling$sp %>% unique %>% length
seedling_tlp$sp %>% unique %>% length

nrow(seedling_all)
nrow(seedling)
nrow(seedling_tlp)

length(unique(seedling_all$sp))
length(unique(seedling$sp))
length(unique(seedling_tlp$sp))
length(unique(shared_tlp$Species))

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

#Xd <- cbind(rep(1, nrow(seedling_wpd)),
#            seedling_wpd[,c("S_scon",
#                            "SC_acon",
#                            "S_shet",
#                            "SC_ahet",
#                            "S_soilpc1",
#                            "S_soilpc2",
#                            "S_soilpc3",
#                            "S_log_h1")])

Xd <- cbind(rep(1, nrow(seedling_tlpd)),
            seedling_tlpd[,c("S_scon",
                             "SC_acon",
                             "S_shet",
                             "SC_ahet",
                             "S_log_h1")])

colnames(Xd)[1] <- "Int"

shared_tlpd2 <- shared_tlp %>%
  filter(Species %in% seedling_tlpd$sp)

dim(Xd)

Ud <- cbind(rep(1,length(unique(seedling_tlpd$sp))),
            shared_tlpd2$tlp)

Ud[, 2] <- scale(Ud[,2])

dim(Ud)

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
                   # x = t(X),
                   u = Ud)
str(list_dat_d)

fit <- stan(file = model_path,
            data = list_dat_d,
            verbose = TRUE,
            iter = n_iter,
            warmup = n_warm,
            thin = n_thin,
            chains =  n_chains,
            refresh = 200,
            control = list(adapt_delta = a_delta, max_treedepth = 20))

print(fit, pars = c("gamma", "sigma", "L_sigma", "lp__"))

save_name <- str_c("./data/", dry, "_spab_", n_ab, "_", model_name, ".rda")
print(save_name)

save.image(save_name)

print("MCMC done!!")
