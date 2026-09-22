library(flowmix) 

# Load the data
datobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper.RDS"))
datobj |> list2env(envir = .GlobalEnv) |> invisible()

# Estimation settings
numclust <- 10
maxdev <- 0.5

X_dir <- file.path(
  "data", 
  "X_variations"
)

X_pc <- readRDS(file.path(X_dir, "X_pc_9_nh_NA_seed_NA_ofold_NA_ifold_NA.RDS"))
load(file.path("2_application", "21_all_data", "linear_settings.Rdata"))

res <- flowmix_once(
  ylist = ylist, 
  X = X_pc, 
  countslist = countslist, 
  numclust = numclust, 
  prob_lambda = prob_lambda, 
  mean_lambda = mean_lambda, 
  verbose = TRUE, 
  maxdev = maxdev, 
  seed = seed
)

save(res, file = file.path("2_application", "21_all_data", "linear_fit.Rdata"))
