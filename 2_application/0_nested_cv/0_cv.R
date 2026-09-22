library(flowmix)
library(magrittr)
library(parallel)

# Load all data
datobj <- readRDS(file.path("data", "MGL1704-hourly-paper.RDS"))
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()

# Get indices of out-of-sample data
X_dir <- file.path(
  "data", 
  "X_variations"
)

load(file.path("data", "ofolds__ifolds__ifolds_inner_inds.Rdata"))

load_X <- function(readdir, n_PCs = "NA", n_h = "NA", seed = "NA", ofold = "NA", ifold = "NA") {

  file_name <- paste("X_pc", n_PCs, "nh", n_h, "seed", seed, "ofold", ofold, "ifold", ifold, sep = "_")
  file_name <- paste0(file_name, ".RDS")
  file_name <- file.path(readdir, file_name)

  X <- readRDS(file_name)
}

############################

# Calculate out-of-sample negative log-likelihoods (OOS NLLs)
job_grid <- expand.grid(
  list(
    ofold = 1:5,
    n_h = c("70", "NA")
  )
)

nlls <- mclapply(1:nrow(job_grid), function(ijob) {
  n_h <- job_grid[ijob,"n_h"]
  ofold <- job_grid[ijob,"ofold"]
  seed <- if(n_h == "NA") "NA" else "1"

  X <- load_X(X_dir, n_PCs = 9, n_h = n_h, seed = seed, ofold = ofold)

  outsample_inds <- ofolds[[ofold]]
  insample_inds <- unlist(ofolds[-ofold])

  X_in <- X[insample_inds,]
  X_out <- X[outsample_inds,]
  ylist_in <- ylist[insample_inds]
  ylist_out <- ylist[outsample_inds]
  countslist_in <- countslist[insample_inds]
  countslist_out <- countslist[outsample_inds]

  load(file.path("2_application", "20_nested_cv", paste("nh", n_h, "seed", seed, "ofold", ofold, "settings.Rdata", sep = "_")))

  res <- flowmix_once(
    ylist = ylist_in, 
    X = X_in, 
    countslist = countslist_in, 
    numclust = 10, 
    prob_lambda = prob_lambda, 
    mean_lambda = mean_lambda, 
    maxdev = 0.5, 
    seed = seed
  )

  pred <- predict(res, newx = X_out)
  nll <- objective_newdat(pred, ylist_out, countslist_out)

  return(nll)
}, mc.cores = 6)
traceback()

nl_oos_nlls <- nlls[1:5]
l_oos_nlls <- nlls[6:10]

mean(unlist(l_oos_nlls)) # 3.486
mean(unlist(nl_oos_nlls)) # 3.414
