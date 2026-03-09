library(flowmix) 
library(parallel)
library(parallelly)
library(RhpcBLASctl)

blas_set_num_threads(1)
omp_set_num_threads(1)

options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
arrnum <- args[1]

mc.cores <- availableCores()

datobj <- file.path("data", "MGL1704-hourly-paper.RDS") |> 
  readRDS()

ylist <- datobj$ylist 
countslist <- datobj$countslist

# Estimation settings
numclust <- 10
maxdev <- 0.5

# Cross-validation settings
cv_gridsize <- 10
nfold <- 5
blocksize <- 20
nrep <- 10
max_mean_lambda <- 40
max_prob_lambda <- 2 

# The one_job function expects the data point indices 
# from the original data set (not subsetted, like we are doing here), 
# so we need to use the rownames to store the original indices.
names(ylist) <- 1:length(ylist)
names(countslist) <- 1:length(countslist)

jobs <- list(
  n_PCs = c(18, 27, 37), 
  seed = 2:5
)

job_grid <- expand.grid(jobs)
n_PCs <- job_grid[arrnum,"n_PCs"]
seed <- job_grid[arrnum,"seed"]

experiment <- paste0("nl_", n_PCs, "_70_", seed)
destin <- file.path("2_application", "23_PCs", "results", experiment)
if(!dir.exists(destin)) dir.create(destin, recursive = TRUE)

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

seedtab <- read.csv(
  file.path(
    "2_application", "23_PCs", "seedtabs", paste0("nl_", n_PCs, "_70_", seed, ".csv")
  )
)

X <- load_X(X_dir, n_PCs = n_PCs, n_h = 70, seed = seed)

maxres <- get_max_lambda(
  destin, 
  "maxres.Rdata",
  ylist, 
  countslist, 
  X, 
  numclust, 
  maxdev, 
  max_prob_lambda,
  max_mean_lambda
)

prob_lambdas <- logspace(min = 1e-4, max = maxres$alpha, length = cv_gridsize)
mean_lambdas <- logspace(min = 1e-4, max = maxres$beta, length = cv_gridsize)

save(prob_lambdas, file = file.path(destin, "prob_lambdas.Rdata"))
save(mean_lambdas, file = file.path(destin, "mean_lambdas.Rdata"))

save(
  prob_lambdas,
  mean_lambdas, 
  nfold, 
  nrep, 
  cv_gridsize, 
  file = file.path(destin, "meta.Rdata")
)

print(paste0("wrote meta data to ", file.path(destin, 'meta.Rdata')))

# Cross-validation
iimat <- make_iimat(
  cv_gridsize = cv_gridsize, 
  nfold = nfold, 
  nrep = nrep
)

folds <- make_cv_folds(ylist, nfold, blocksize)

empty <- mclapply(
  1:nrow(iimat), 
  function(ii) {
    ialpha <- iimat[ii,"ialpha"]
    ibeta <- iimat[ii,"ibeta"]
    ifold <- iimat[ii,"ifold"]
    irep <- iimat[ii,"irep"]

    cat("\r", ii, "out of", nrow(iimat), "cross-validation jobs.")

    cur_X <- load_X(X_dir, n_PCs = n_PCs, n_h = 70, seed = seed, ofold = ifold)

    one_job(
      ialpha = ialpha,
      ibeta = ibeta, 
      ifold = ifold, 
      irep = irep, 
      folds = folds, 
      destin = destin, 
      mean_lambdas = mean_lambdas, 
      prob_lambdas = prob_lambdas, 
      seedtab = seedtab, 
      ylist = ylist, 
      countslist = countslist, 
      X = cur_X,
      maxdev = maxdev, 
      numclust = numclust
    )
  }, 
  mc.cores = mc.cores, 
  mc.preschedule = FALSE
)

cv.flowmix(
  ylist = ylist, 
  countslist = countslist, 
  X = X, 
  destin = destin, 
  mean_lambdas = mean_lambdas, 
  prob_lambdas = prob_lambdas,
  iimat = NULL, 
  maxdev = maxdev, 
  numclust = numclust, 
  nfold = nfold, 
  nrep = nrep, 
  verbose = TRUE, 
  refit = TRUE, 
  save_meta = FALSE, 
  mc.cores = mc.cores,
  blocksize = blocksize, 
  folds = folds, 
  seedtab = seedtab
)

cv_summary(destin = destin, save = TRUE, filename = paste0(experiment, "_pcr_summary.RDS"))