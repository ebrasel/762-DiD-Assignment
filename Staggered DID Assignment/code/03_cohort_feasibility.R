k <- 5

units <- unique(dt[, .(state, E)])

cohorts <- sort(unique(units$E[is.finite(units$E)]))

feasibility <- function(k) {
  
  results <- data.table(
    a = cohorts,
    window = k
  )
  
  results[, IC1a :=
            a - k >= min(dt$year) &
            a + k <= max(dt$year)]
  
  results[, ND_a := sapply(a, function(x) {
    sum(units$E == x, na.rm = TRUE)
  })]
  
  results[, NC_a := sapply(a, function(x) {
    sum(units$E > x + k, na.rm = TRUE)
  })]
  
  results[, IC2a := NC_a > 0]
  
  results[, feasible := IC1a & IC2a]
  
  return(results)
}

feasible_5 <- feasibility(5)
feasible_10 <- feasibility(10)

print(feasible_5)
print(feasible_10)

saveRDS(dt, here("temp_data", "panel.rds"))

saveRDS(A, here("temp_data", "cohorts.rds"))

fwrite(
  rbind(feasible_5, feasible_10),
  here("output", "03_feasibility.csv")
)