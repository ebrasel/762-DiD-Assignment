library(data.table)
library(fixest)
library(here)
library(ggplot2)


# LOAD DATA AND MODELS

dt <- readRDS(
  here("temp_data", "imputed_panel.rds")
)

A <- readRDS(
  here("temp_data", "cohorts.rds")
)

event_results <- fread(
  here("output", "05_stacked_event_results.csv")
)

event_ATT <- fread(
  here("output", "08_BJS_event_ATT.csv")
)

bjs_pretrend <- readRDS(
  here("temp_data", "bjs_pretrend.rds")
)

k <- 5


# EXTRACT STACKED DID ESTIMATES

stacked_results <- event_results[
  model == "Q-weighted FE",
  .(e, estimate, lower, upper, model)
]

stacked_results[, model := "Stacked DID"]

setorder(stacked_results, e)

print(stacked_results)


# EXTRACT BJS POST-TREATMENT ESTIMATES

bjs_post <- event_ATT[, .(
  e = Kit,
  estimate = ATT,
  lower = NA_real_,
  upper = NA_real_,
  model = "BJS imputation"
)]

print(bjs_post)


# EXTRACT BJS PRE-TREATMENT ESTIMATES

b <- coef(bjs_pretrend)

V <- vcov(bjs_pretrend)

lead_names <- paste0("lead", k:1)

stopifnot(all(lead_names %in% names(b)))

bjs_pre <- data.table(
  e = -k:-1,
  estimate = as.numeric(b[lead_names]),
  se = sqrt(diag(V)[lead_names])
)

critical <- qt(
  0.975,
  df = uniqueN(dt$state) - 1
)

bjs_pre[, lower := estimate - critical * se]

bjs_pre[, upper := estimate + critical * se]

bjs_pre[, model := "BJS pre-trends"]

bjs_pre[, se := NULL]

print(bjs_pre)


# COMBINE EVENT STUDY ESTIMATES

comparison <- rbindlist(
  list(
    stacked_results,
    bjs_pre,
    bjs_post
  ),
  use.names = TRUE
)

setorder(comparison, model, e)

print(comparison)

fwrite(
  comparison,
  here("output", "11_comparison.csv")
)


# PLOT EVENT STUDY COMPARISON

p3 <- ggplot(
  comparison,
  aes(
    x = e,
    y = estimate,
    color = model,
    group = model
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = -0.5,
    linetype = "dotted"
  ) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = -k:k) +
  labs(
    title = "Comparison of Staggered DID Estimators",
    subtitle = "Unilateral divorce laws and female suicide",
    x = "Years relative to adoption",
    y = "Effect per 100,000 women",
    color = "Estimator"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p3)

ggsave(
  here("output", "figure_03_comparison.png"),
  plot = p3,
  width = 8,
  height = 5,
  dpi = 300
)


# COMPARE POST-TREATMENT EVENT-TIME ESTIMATES

post_comparison <- merge(
  stacked_results[
    e >= 0,
    .(e, stacked_ATT = estimate)
  ],
  bjs_post[
    ,
    .(e, BJS_ATT = estimate)
  ],
  by = "e"
)

post_comparison[, difference := BJS_ATT - stacked_ATT]

setorder(post_comparison, e)

print(post_comparison)

fwrite(
  post_comparison,
  here("output", "14_post_comparison.csv")
)


# COMPARE TRIMMED AVERAGE TREATMENT EFFECTS

ATT_comparison <- data.table(
  estimator = c(
    "Q-weighted stacked DID",
    "BJS imputation"
  ),
  trimmed_ATT = c(
    mean(post_comparison$stacked_ATT),
    mean(post_comparison$BJS_ATT)
  )
)

print(ATT_comparison)

fwrite(
  ATT_comparison,
  here("output", "15_ATT_comparison.csv")
)
