library(data.table)
library(fixest)
library(here)
library(ggplot2)


# LOAD DATA

dt <- readRDS(
  here("temp_data", "imputed_panel.rds")
)

A <- readRDS(
  here("temp_data", "cohorts.rds")
)

event_ATT <- fread(
  here("output", "08_BJS_event_ATT.csv")
)

k <- 5

stopifnot(!anyNA(dt$Ei))
stopifnot(!anyNA(dt$Dit))
stopifnot(!anyNA(dt$Y0_hat))


# SET UP STATE-CLUSTERED BOOTSTRAP

set.seed(123)

B <- 999

states <- unique(dt$state)

S <- length(states)

boot_results <- matrix(
  NA_real_,
  nrow = B,
  ncol = k + 1
)

colnames(boot_results) <- 0:k


# STATE-CLUSTERED BOOTSTRAP

for (b in 1:B) {
  
  sampled_states <- sample(
    states,
    size = S,
    replace = TRUE
  )
  
  boot_list <- vector("list", S)
  
  for (j in seq_len(S)) {
    
    temp <- copy(
      dt[state == sampled_states[j]]
    )
    
    temp[, boot_state := j]
    
    boot_list[[j]] <- temp
    
  }
  
  boot_dt <- rbindlist(boot_list)
  
  boot_omega0 <- boot_dt[Dit == 0]
  
  boot_omega1 <- boot_dt[Dit == 1]
  
  if (nrow(boot_omega0) == 0 ||
      nrow(boot_omega1) == 0) next
  
  missing_states <- setdiff(
    unique(boot_dt$boot_state),
    unique(boot_omega0$boot_state)
  )
  
  missing_years <- setdiff(
    unique(boot_dt$year),
    unique(boot_omega0$year)
  )
  
  if (length(missing_states) > 0 ||
      length(missing_years) > 0) next
  
  model_b <- tryCatch(
    feols(
      Y ~ 1 | boot_state + year,
      data = boot_omega0,
      notes = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(model_b)) next
  
  predictions <- tryCatch(
    predict(
      model_b,
      newdata = boot_dt
    ),
    error = function(e) NULL
  )
  
  if (is.null(predictions)) next
  
  if (anyNA(predictions)) next
  
  boot_dt[, Y0_hat := predictions]
  
  boot_dt[, tau_hat := ifelse(
    Dit == 1,
    Y - Y0_hat,
    NA_real_
  )]
  
  for (h in 0:k) {
    
    effects <- boot_dt[
      Dit == 1 &
        Ei %in% A &
        Kit == h,
      tau_hat
    ]
    
    if (length(effects) > 0 &&
        all(is.finite(effects))) {
      
      boot_results[b, h + 1] <- mean(effects)
      
    }
    
  }
  
  if (b %% 100 == 0) {
    
    print(paste("Bootstrap replication:", b))
    
  }
  
}


# CHECK BOOTSTRAP REPLICATIONS

valid <- complete.cases(boot_results)

B_valid <- sum(valid)

B_invalid <- B - B_valid

print(B_valid)
print(B_invalid)

stopifnot(B_valid > 1)

if (B_valid < 0.95 * B) {
  
  warning("More than 5% of bootstrap replications failed.")
  
}

boot_valid <- boot_results[
  valid,
  ,
  drop = FALSE
]

bootstrap_summary <- data.table(
  statistic = c(
    "Requested replications",
    "Successful replications",
    "Unsuccessful replications"
  ),
  value = c(
    B,
    B_valid,
    B_invalid
  )
)

print(bootstrap_summary)

fwrite(
  bootstrap_summary,
  here("output", "16_bootstrap_summary.csv")
)


# CALCULATE BOOTSTRAP STANDARD ERRORS

boot_se <- apply(
  boot_valid,
  2,
  sd
)

print(boot_se)

event_ATT[, SE := boot_se[
  match(Kit, 0:k)
]]

stopifnot(!anyNA(event_ATT$SE))


# CALCULATE 95% CONFIDENCE INTERVALS

event_ATT[, lower := ATT - 1.96 * SE]

event_ATT[, upper := ATT + 1.96 * SE]

event_ATT[, z_stat := ATT / SE]

event_ATT[, p_value := 2 * pnorm(
  -abs(z_stat)
)]

print(event_ATT)

fwrite(
  event_ATT,
  here("output", "12_BJS_inference.csv")
)


# INFERENCE FOR THE TRIMMED ATT

ATT_trimmed <- mean(event_ATT$ATT)

boot_trimmed <- rowMeans(boot_valid)

SE_trimmed <- sd(boot_trimmed)

trimmed_inference <- data.table(
  estimator = "BJS trimmed ATT",
  ATT = ATT_trimmed,
  SE = SE_trimmed
)

trimmed_inference[, lower := ATT - 1.96 * SE]

trimmed_inference[, upper := ATT + 1.96 * SE]

trimmed_inference[, z_stat := ATT / SE]

trimmed_inference[, p_value := 2 * pnorm(
  -abs(z_stat)
)]

print(trimmed_inference)

fwrite(
  trimmed_inference,
  here("output", "17_BJS_trimmed_inference.csv")
)


# PLOT BJS EVENT STUDY WITH CONFIDENCE INTERVALS

p4 <- ggplot(
  event_ATT,
  aes(x = Kit, y = ATT)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line() +
  geom_point() +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.15
  ) +
  scale_x_continuous(breaks = 0:k) +
  labs(
    title = "BJS Imputation Event Study",
    subtitle = "95% state-clustered bootstrap confidence intervals",
    x = "Years since adoption",
    y = "Effect per 100,000 women"
  ) +
  theme_minimal()

print(p4)

ggsave(
  here("output", "figure_04_BJS_inference.png"),
  plot = p4,
  width = 8,
  height = 5,
  dpi = 300
)


# LOAD STACKED DID AND BJS PRE-TREND ESTIMATES

comparison <- fread(
  here("output", "11_comparison.csv")
)

comparison <- comparison[
  model != "BJS imputation"
]

bjs_inference <- event_ATT[, .(
  e = Kit,
  estimate = ATT,
  lower,
  upper,
  model = "BJS imputation"
)]


# COMBINE ESTIMATES WITH CONFIDENCE INTERVALS

comparison_inference <- rbindlist(
  list(
    comparison,
    bjs_inference
  ),
  use.names = TRUE
)

setorder(comparison_inference, model, e)

print(comparison_inference)

fwrite(
  comparison_inference,
  here("output", "18_comparison_inference.csv")
)


# PLOT COMPARATIVE EVENT STUDY WITH CONFIDENCE INTERVALS

p5 <- ggplot(
  comparison_inference,
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
  geom_line(
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    position = position_dodge(width = 0.4)
  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.15,
    position = position_dodge(width = 0.4)
  ) +
  scale_x_continuous(breaks = -k:k) +
  labs(
    title = "Comparison of Staggered DID Estimators",
    subtitle = "95% confidence intervals",
    x = "Years relative to adoption",
    y = "Effect per 100,000 women",
    color = "Estimator"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p5)

ggsave(
  here("output", "figure_05_comparison_inference.png"),
  plot = p5,
  width = 9,
  height = 6,
  dpi = 300
)


# SAVE BOOTSTRAP RESULTS

saveRDS(
  boot_results,
  here("temp_data", "BJS_bootstrap.rds")
)
