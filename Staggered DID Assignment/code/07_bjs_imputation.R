library(data.table)
library(fixest)
library(here)
library(ggplot2)


# PART 2.1: DEFINE THE SETS

dt <- readRDS(here("temp_data", "panel.rds"))

A <- readRDS(here("temp_data", "cohorts.rds"))

k <- 5

dt[, E := as.numeric(E)]

dt[never_treated == TRUE, E := Inf]

dt <- dt[always_treated == FALSE]

stopifnot(!anyNA(dt$E))

dt[, Ei := E]

dt[, Kit := year - Ei]

dt[, Dit := as.integer(year >= Ei)]

stopifnot(all(dt$Dit == dt$D))

omega0 <- dt[Dit == 0]

omega1 <- dt[Dit == 1]

N0 <- nrow(omega0)

N1 <- nrow(omega1)

N_never <- uniqueN(dt[is.infinite(Ei), state])

stopifnot(N0 + N1 == nrow(dt))

sets_table <- data.table(
  statistic = c(
    "Total states",
    "Total observations",
    "Untreated observations",
    "Treated observations",
    "Never-treated states"
  ),
  value = c(
    uniqueN(dt$state),
    nrow(dt),
    N0,
    N1,
    N_never
  )
)

print(sets_table)

fwrite(
  sets_table,
  here("output", "07_imputation_sets.csv")
)


# PART 2.2: FIXED EFFECTS (STEP 1)

bjs_model <- feols(
  Y ~ 1 | state + year,
  data = omega0,
  cluster = ~state
)

summary(bjs_model)


# PRE-TRENDS MODEL

for (h in 1:5) {
  
  dt[, (paste0("lead", h)) :=
       as.integer(Kit == -h)]
  
}

stopifnot(
  all(paste0("lead", 1:5) %in% names(dt))
)

bjs_pretrend <- feols(
  Y ~ lead5 + lead4 + lead3 + lead2 + lead1 |
    state + year,
  data = dt[Dit == 0],
  cluster = ~state
)

summary(bjs_pretrend)

pretrend_test <- wald(
  bjs_pretrend,
  keep = "lead"
)

print(pretrend_test)

pretrend_table <- as.data.table(
  coeftable(bjs_pretrend),
  keep.rownames = "term"
)

print(pretrend_table)

fwrite(
  pretrend_table,
  here("output", "13_BJS_pretrend.csv")
)


# PART 2.3: IMPUTE UNTREATED OUTCOMES (STEP 2)

missing_states <- setdiff(
  unique(dt$state),
  unique(omega0$state)
)

missing_years <- setdiff(
  unique(dt$year),
  unique(omega0$year)
)

print(missing_states)
print(missing_years)

stopifnot(length(missing_states) == 0)
stopifnot(length(missing_years) == 0)

dt[, Y0_hat := predict(
  bjs_model,
  newdata = dt
)]

print(sum(is.na(dt$Y0_hat)))

print(
  dt[is.na(Y0_hat), .(state, year, Ei, Kit, Dit)]
)

stopifnot(!anyNA(dt$Y0_hat))

dt[, tau_hat := ifelse(
  Dit == 1,
  Y - Y0_hat,
  NA_real_
)]

stopifnot(!anyNA(dt[Dit == 1, tau_hat]))

summary(dt$Y0_hat)
summary(dt$tau_hat)


# IMPUTATION USING THE PRE-TRENDS MODEL

prediction_data <- copy(dt)

for (h in 1:5) {
  
  prediction_data[, (paste0("lead", h)) := 0L]
  
}

stopifnot(
  all(paste0("lead", 1:5) %in%
        names(prediction_data))
)

dt[, Y0_pretrend := predict(
  bjs_pretrend,
  newdata = prediction_data
)]

print(sum(is.na(dt$Y0_pretrend)))

stopifnot(!anyNA(dt$Y0_pretrend))

dt[, tau_pretrend := ifelse(
  Dit == 1,
  Y - Y0_pretrend,
  NA_real_
)]

stopifnot(!anyNA(dt[Dit == 1, tau_pretrend]))

summary(dt$Y0_pretrend)
summary(dt$tau_pretrend)


# PART 2.4: AVERAGE TREATMENT EFFECTS (STEP 3)


# PART 2.4(a): OVERALL ATT

ATT_overall <- dt[
  Dit == 1,
  mean(tau_hat)
]

print(ATT_overall)


# PART 2.4(b): EVENT-TIME ATT

event_ATT <- dt[
  Dit == 1 &
    Ei %in% A &
    Kit %in% 0:k,
  .(
    ATT = mean(tau_hat),
    N = .N
  ),
  by = Kit
]

setorder(event_ATT, Kit)

stopifnot(nrow(event_ATT) == k + 1)
stopifnot(all(event_ATT$Kit == 0:k))
stopifnot(uniqueN(event_ATT$N) == 1)

print(event_ATT)

ATT_trimmed <- mean(event_ATT$ATT)

print(ATT_trimmed)

fwrite(
  event_ATT,
  here("output", "08_BJS_event_ATT.csv")
)


# PART 2.4(c): CALENDAR-TIME ATT

calendar_ATT <- dt[
  Dit == 1,
  .(
    ATT = mean(tau_hat),
    N_treated = .N
  ),
  by = year
]

setorder(calendar_ATT, year)

print(calendar_ATT)

fwrite(
  calendar_ATT,
  here("output", "09_BJS_calendar_ATT.csv")
)


# PLOT CALENDAR-TIME ATT

p2 <- ggplot(
  calendar_ATT,
  aes(x = year, y = ATT)
) +
  geom_line() +
  geom_point() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Calendar-Time Average Treatment Effects",
    subtitle = "BJS imputation estimator",
    x = "Calendar year",
    y = "Effect per 100,000 women"
  ) +
  theme_minimal()

print(p2)

ggsave(
  here("output", "figure_02_calendar_ATT.png"),
  plot = p2,
  width = 8,
  height = 5,
  dpi = 300
)


# SUMMARY TABLE OF ATT ESTIMATES

ATT_table <- data.table(
  estimator = c(
    "BJS overall ATT",
    "BJS trimmed ATT"
  ),
  ATT = c(
    ATT_overall,
    ATT_trimmed
  )
)

print(ATT_table)

fwrite(
  ATT_table,
  here("output", "10_BJS_average_ATT.csv")
)


# SAVE INTERMEDIATE DATA AND MODELS

saveRDS(
  dt,
  here("temp_data", "imputed_panel.rds")
)

saveRDS(
  bjs_model,
  here("temp_data", "bjs_model.rds")
)

saveRDS(
  bjs_pretrend,
  here("temp_data", "bjs_pretrend.rds")
)
