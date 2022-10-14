library(tidyverse)
library(targets)
hoge <- tar_manifest()


hoge |>
  filter(str_detect(name, "loo"))


mcmc_names <- expand_grid(a1 = c("phy", "het"), a2 = c("season", "rain"), a3 = c("ab", "ba", "both")) |>
  mutate(model = str_c("fit_mcmc_multilevel_logistic", a1, a2, a3, sep = "_")) |>
  pull(model)

mcmc_names <- expand_grid(a1 = c("phy", "het"), a2 = c("season", "rain"), a3 = c("ab", "ba", "both")) |>
  mutate(model = str_c("fit_mcmc_multilevel_logistic", a1, a2, a3, sep = "_")) |>
  pull(model)

data_names <- expand_grid(a1 = c("phy", "het"), a2 = c("season", "rain"), a3 = c("ab", "ba", "both")) |>
  mutate(data = str_c(a1, a2, a3, sep = "_")) |>
  pull(data)

library(tidyverse)
  seedling <- read_csv("data/seedling.csv") |>
    janitor::clean_names()
  traits <- read_csv("data/traits.csv") |>
    janitor::clean_names()

  d <- read_csv("data/seedling.csv")

  d |>
    dplyr::select(plot, tag, latin, h1:shet, aphy:sphy, RF) |>

dplyr::select(-note1, -dat)

targets::tar_read(phy_season_both) |> str()

hoge <- targets::tar_read(phy_season_both)
hoge$u
