dt <- readRDS(here("temp_data", "panel.rds"))

A <- readRDS(here("temp_data", "cohorts.rds"))

k <- 5

dt[is.na(E) & never_treated == TRUE, E := Inf]

sum(is.na(dt$E))

stack_list <- vector("list", length(A))

for (i in seq_along(A)) {
  
  a <- A[i]
  
  temp <- copy(dt[
    (E == a | E > a + k) &
      year >= a - k &
      year <= a + k
  ])
  
  temp[, experiment := a]
  
  temp[, Dsa := as.integer(E == a)]
  
  temp[, e := year - a]
  
  stack_list[[i]] <- temp
}

stack <- rbindlist(stack_list)

setorder(stack, experiment, state, year)

check <- stack[, .N, by = .(state, experiment)]

table(check$N)

stopifnot(all(check$N == 2 * k + 1))

state_frequency <- stack[, .(
  N_experiments = uniqueN(experiment)
), by = state]

print(state_frequency)

table(state_frequency$N_experiments)

nrow(state_frequency[N_experiments > 1])

stack[, state_exp := interaction(
  state, experiment, drop = TRUE
)]

stack[, exp_event := interaction(
  experiment, e, drop = TRUE
)]

saveRDS(stack, here("temp_data", "stack.rds"))