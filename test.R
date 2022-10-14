library(tidyverse)
library(targets)
hoge <- tar_manifest()


hoge |>
  filter(str_detect(name, "loo"))


mcmc_names <- expand_grid(a1 = c("phy", "het"), a2 = c("season", "rain"), a3 = c("ab", "ba", "both")) |>
  mutate(model = str_c("fit_mcmc_multilevel_logistic", a1, a2, a3, sep = "_")) |>
  pull(model)
