#Building stacked regressions

stack <- readRDS(
  here("temp_data", "stack_weighted.rds")
)

# Model 1

m1 <- feols(
  Y ~ Dsa + i(e, ref = -1) + i(e, Dsa, ref = -1),
  data = stack,
  cluster = ~state
)

# Model 2

m2 <- feols(
  Y ~ Dsa + i(e, ref = -1) + i(e, Dsa, ref = -1),
  data = stack,
  weights = ~Q,
  cluster = ~state
)

# Model 3

m3 <- feols(
  Y ~ i(e, Dsa, ref = -1) | state_exp + exp_event,
  data = stack,
  cluster = ~state
)

# Model 4

m4 <- feols(
  Y ~ i(e, Dsa, ref = -1) | state_exp + exp_event,
  data = stack,
  weights = ~Q,
  cluster = ~state
)

etable(m1, m2, m3, m4)

saveRDS(list(m1, m2, m3, m4),
        here("temp_data", "stack_models.rds"))

#Event coefficients
extract_event <- function(model, model_name) {
  
  ct <- as.data.table(
    coeftable(model),
    keep.rownames = "term"
  )
  
  ct <- ct[grepl(":Dsa$", term)]
  
  results <- data.table(
    e = as.numeric(
      sub("e::(-?[0-9]+):Dsa", "\\1", ct$term)
    ),
    estimate = ct[[2]],
    se = ct[[3]]
  )
  
  critical <- qt(
    0.975, df = uniqueN(stack$state) - 1
  )
  
  results[, lower := estimate - critical * se]
  results[, upper := estimate + critical * se]
  
  results <- rbind(
    results,
    data.table(
      e = -1,
      estimate = 0,
      se = 0,
      lower = 0,
      upper = 0
    )
  )
  
  results[, model := model_name]
  
  setorder(results, e)
  
  return(results)
}

event_results <- rbindlist(list(
  extract_event(m1, "Unweighted saturated"),
  extract_event(m2, "Q-weighted saturated"),
  extract_event(m3, "Unweighted FE"),
  extract_event(m4, "Q-weighted FE")
))

print(event_results)

fwrite(event_results,
       here("output", "05_stacked_event_results.csv"))

#plot event studies

p1 <- ggplot(
  event_results,
  aes(x = e, y = estimate,
      color = model, group = model)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = -0.5, linetype = "dotted") +
  geom_line(position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.15,
    position = position_dodge(width = 0.5)
  ) +
  scale_x_continuous(breaks = -5:5) +
  labs(
    title = "Effect of Unilateral Divorce on Female Suicide",
    subtitle = "Four stacked difference-in-differences specifications",
    x = "Years relative to adoption",
    y = "Effect per 100,000 women",
    color = "Specification"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)

ggsave(
  here("output", "figure_01_stacked_DID.png"),
  p1, width = 9, height = 6, dpi = 300
)

#Post period average

post_average <- event_results[e >= 0, .(
  post_ATT = mean(estimate)
), by = model]

print(post_average)

fwrite(post_average,
       here("output", "06_stacked_post_ATT.csv"))