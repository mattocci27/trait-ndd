library(targets)
library(tidyverse)

# Fetch the model ID from the command-line arguments
args <- commandArgs(trailingOnly = TRUE)
model_id <- as.numeric(args[1])

values <- expand_grid(
  season = c("dry", "wet"),
  rain = c("norain", "rain", "intrain", "intrain2", "intrain3", "intrain4"),
  sp_pred = c("nlog", "ab", "ba", "pc12")
  )

model_name <- values |>
      mutate(model_name = paste(season, rain, sp_pred, sep = "_")) |>
      pull(model_name)

tar_name <- paste0("fit_mcmc_suv_ind_", model_name[model_id])

print(tar_name)

model_id <- ifelse(model_id %% 9 == 0, 9, model_id %% 9)
model_id <- model_id - 1

Sys.sleep(model_id * 20)

# Run the specific target corresponding to the model ID
targets::tar_make(callr_function = NULL,
  names = tar_name)

