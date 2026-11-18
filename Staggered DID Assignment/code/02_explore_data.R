pacman::p_load(here, fixest, data.table, ggplot2, collapse, bacondecomp)

dt <- fread(here("temp_data", "divorce_clean.csv"))

summary_table <- data.table(
  statistic = c(
    "Number of states",
    "Number of years",
    "Number of observations",
    "First year",
    "Last year",
    "Mean female suicide rate",
    "SD female suicide rate",
    "Minimum female suicide rate",
    "Maximum female suicide rate"
  ),
  value = c(
    uniqueN(dt$state),
    uniqueN(dt$year),
    nrow(dt),
    min(dt$year),
    max(dt$year),
    mean(dt$Y, na.rm = TRUE),
    sd(dt$Y, na.rm = TRUE),
    min(dt$Y, na.rm = TRUE),
    max(dt$Y, na.rm = TRUE)
  )
)

print(summary_table)

fwrite(summary_table,
       here("output", "01_summary_statistics.csv"))

dt[, E := if (any(D == 1)) {
  min(year[D == 1])
} else {
  Inf
}, by = state]

dt[, always_treated := all(D == 1),
   by = state]

dt[, never_treated := all(D == 0),
   by = state]

adoption <- unique(dt[, .(
  state, E, always_treated, never_treated
)])

adoption_table <- adoption[, .(
  N_states = .N
), by = .(E, always_treated, never_treated)]

setorder(adoption_table, E)

print(adoption_table)

fwrite(adoption_table,
       here("output", "02_adoption_cohorts.csv"))

dt <- dt[always_treated == FALSE]

saveRDS(dt, here("temp_data", "panel.rds"))
