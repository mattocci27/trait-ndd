
trait_pair <- function(trait_log, full = TRUE) {
  trait_log <- trait_log |>
    dplyr::select(-sp_code) |>
    rename(LDMC = ldmc) |>
    rename(WD = wd) |>
    rename(SDMC = sdmc) |>
    rename(Chl = chl) |>
    rename("\u03b4C[13]" = c13) |>
    rename("\u03c0[tlp]" = tlp) |>
    rename(C = c) |>
    rename(N = n) |>
    rename(CN = cn) |>
    rename(`ln~LA` = log_la) |>
    rename(`ln~SLA` = log_sla) |>
    rename(`ln~LT` = log_lt)

  if(!full) {
    trait_log <- trait_log  |>
      dplyr::select(-CN, -WD)
  }

  trait_log  |>
    ggpairs(
      aes(alpha = 0.6),
      upper = list(continuous = wrap("cor", size = 2.5)),
      labeller = "label_parsed"
    ) +
    theme_bw() +
    theme(text = element_text(family = "Arial", size = 6),
      strip.text = element_text(size = 8)
    )
}
