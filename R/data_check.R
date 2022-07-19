
trait_pair <- function(trait_log, full = TRUE) {
  trait_log <- trait_log |>
    dplyr::select(-sp_code) |>
    rename(LDMC = ldmc) |>
    rename(WD = wd) |>
    rename(SDMC = sdmc) |>
    rename(Chl = chl) |>
    rename(C13 = c13) |>
    rename(C = c) |>
    rename(N = n) |>
    rename(CN = cn) |>
    rename(`log[LA]` = log_la) |>
    rename(`log[SLA]` = log_sla) |>
    rename(`log[LT]` = log_lt)

  if(!full) {
    trait_log <- trait_log  |>
      dplyr::select(-CN, -WD)
  }

  trait_log  |>
    ggpairs(
      aes(alpha = 0.6),
      upper = list(continuous = wrap("cor", size = 2.5))
    ) +
    theme(text = element_text(family = "Arial", size = 6),
      strip.text = element_text(size = 8)
    )
}
