pacman::p_load(here, fixest, data.table, ggplot2, collapse, bacondecomp)


data(divorce, package = "bacondecomp")

dt <- data.table(divorce)

dim(dt)
names(dt)
str(dt)

table(dt$sex)
summary(dt$suicide)

grep("div|unil|reform|adopt|post",
     names(dt), value = TRUE)

dt <- dt[sex == 2]

dt[, Y := suicide * 100000]

dt[, state := st]
dt[, D := unilateral]

dt <- dt[, .(state, year, Y, D)]

setorder(dt, state, year)

stopifnot(!anyDuplicated(dt[, .(state, year)]))

summary(dt)

fwrite(dt, here("temp_data", "divorce_clean.csv"))
