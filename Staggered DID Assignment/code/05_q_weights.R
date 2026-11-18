stack <- readRDS(here("temp_data", "stack.rds"))

Q_table <- stack[, .(
  ND_a = uniqueN(state[Dsa == 1]),
  NC_a = uniqueN(state[Dsa == 0])
), by = experiment]

ND <- sum(Q_table$ND_a)
NC <- sum(Q_table$NC_a)

Q_table[, treated_share := ND_a / ND]

Q_table[, control_share := NC_a / NC]

Q_table[, Q_control := treated_share / control_share]

print(Q_table)

stack <- merge(
  stack,
  Q_table[, .(experiment, Q_control)],
  by = "experiment"
)

stack[, Q := ifelse(Dsa == 1, 1, Q_control)]

fwrite(Q_table, here("output", "04_Q_weights.csv"))

saveRDS(stack, here("temp_data", "stack_weighted.rds"))

#checking if weights exceed 1

Q_table[Q_control > 1]

summary(stack$Q)