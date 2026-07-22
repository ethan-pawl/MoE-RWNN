library(flowmix)
library(parallel)
library(parallelly)

# Parallelization controls
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

mc.cores <- availableCores()

# Warnings print as they occur
options(warn = 1)

datobj <- readRDS(
  file.path(
    "4_3dsim", 
    "3dsimdata", 
    "simdata_imean_2_iprob_1_iint_5_dataSeed_1.RDS"
  )
)

ybin_list <- datobj$ybin_list
countslist <- datobj$countslist

# Estimation settings
maxdev <- diff(range(datobj$mean_spec[,1])) / 2 
numclust <- 2
nrep <- 10

# PCA 
n_PCs <- 9

# NN Settings
n_h <- "NA"
seed <- "NA"

# K-fold cross-validation settings
nfold <- 5 
cv_gridsize <- 10
blocksize <- 20
max_prob_lambda <- 24
max_mean_lambda <- 24

seedtab <- read.csv(
  file.path(
    "4_3dsim", 
    "3dsim_seedtab.csv"
  )
)

destin <- file.path(
  "5_competitors", 
  "50_sim", 
  "507_linear", 
  "results"
)

if(!dir.exists(destin)) dir.create(destin)

X_dir <- file.path("data", "X_variations")
load_X <- function(readdir, n_PCs = "NA", n_h = "NA", seed = "NA", ofold = "NA", ifold = "NA") {

  file_name <- paste("X_pc", n_PCs, "nh", n_h, "seed", seed, "ofold", ofold, "ifold", ifold, sep = "_")
  file_name <- paste0(file_name, ".RDS")
  file_name <- file.path(readdir, file_name)

  X <- readRDS(file_name)
}

######################################

X <- load_X(X_dir, n_PCs = n_PCs, n_h = n_h, seed = seed)
maxres <- get_max_lambda(
  destin, 
  "maxres.Rdata", 
  ybin_list, 
  countslist, 
  X, 
  numclust, 
  maxdev, 
  max_prob_lambda, 
  max_mean_lambda
)

prob_lambdas <- logspace(
  min = 0.0001, max = maxres$alpha, length = cv_gridsize
)

mean_lambdas <- logspace(
  min = 0.0001, max = maxres$beta, length = cv_gridsize
)

save(prob_lambdas, file = file.path(destin, "prob_lambdas.Rdata"))
save(mean_lambdas, file = file.path(destin, "mean_lambdas.Rdata"))

save(
  nfold, 
  nrep, 
  cv_gridsize, 
  mean_lambdas, 
  prob_lambdas, 
  file = file.path(destin, "meta.Rdata")
)

print(paste0("wrote meta data to ", file.path(destin, 'meta.Rdata')))

############################################

# Cross-validation

folds <- make_cv_folds(ybin_list, nfold, blocksize)

iimat <- make_iimat(
  cv_gridsize = cv_gridsize, 
  nfold = nfold, 
  nrep = nrep
)

empty <- mclapply(
  1:nrow(iimat), 
  function(ii) {
    ialpha <- iimat[ii,"ialpha"]
    ibeta <- iimat[ii,"ibeta"]
    ifold <- iimat[ii,"ifold"]
    irep <- iimat[ii,"irep"]

    cat("\r", ii, "out of", nrow(iimat), "cross-validation jobs.")

    cur_X <- load_X(X_dir, n_PCs = n_PCs, n_h = n_h, seed = seed, ofold = ifold)

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
      ylist = ybin_list, 
      countslist = countslist, 
      X = cur_X,
      maxdev = maxdev, 
      numclust = numclust
    )
  }, 
  mc.cores = mc.cores, 
  mc.preschedule = FALSE
)

##############################
# Refit

cv.flowmix(
  ylist = ybin_list, 
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

cv_summary(
  destin = destin, 
  save = TRUE, 
  filename = "3dsim_summary_linear.RDS"
)