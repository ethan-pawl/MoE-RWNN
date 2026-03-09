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

load_X <- function(readdir, n_PCs = "NA", n_h = "NA", seed = "NA", ofold = "NA", ifold = "NA") {

  file_name <- paste("X_pc", n_PCs, "nh", n_h, "seed", seed, "ofold", ofold, "ifold", ifold, sep = "_")
  file_name <- paste0(file_name, ".RDS")
  file_name <- file.path(readdir, file_name)

  X <- readRDS(file_name)
  return(X)
}

X_pc <- load_X(X_dir, 9)
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
