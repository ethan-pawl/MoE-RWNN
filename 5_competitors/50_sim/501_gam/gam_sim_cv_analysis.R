cvscore_arr <- readRDS(file.path("5_competitors", "50_sim", "501_gam", "cv-summary.RDS"))

cvscore_arr

min_gamma <- 0.1
max_gamma <- 10000

gamma_grid <- flowmix::logspace(min_gamma, max_gamma, nrow(cvscore_arr))

matplot(gamma_grid, log(cvscore_arr + abs(min(cvscore_arr)) + 1), type = "l")
lines(gamma_grid, log(rowMeans(cvscore_arr) + abs(min(cvscore_arr)) + 1), col = "blue")

gamma_grid[which.min(rowMeans(cvscore_arr))]
# [1] 774.2637