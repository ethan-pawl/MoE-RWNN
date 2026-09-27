library(flowmix) 
library(parallel)
library(parallelly)
library(RhpcBLASctl)

blas_set_num_threads(1)
omp_set_num_threads(1)

options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
arrnum <- as.integer(args[1])

mc.cores <- availableCores()

datobj <- file.path("data", "MGL1704-hourly-paper.RDS") |> 
  readRDS()

ylist <- datobj$ylist 
countslist <- datobj$countslist

# Estimation settings
numclust <- 10
maxdev <- 0.5

# Cross-validation settings
nfold <- 5
blocksize <- 20
nrep <- 10

num_PCs <- c("9", "18", "27", "37")
NNseeds <- c(4, 1, 2, 1)

num_PC <- num_PCs[arrnum]
NNseed <- NNseeds[arrnum]

experiment <- paste0("nl_", num_PC, "_70_", NNseed)
destin <- file.path("2_application", "3_PCs", "results", experiment)
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
    "2_application", "3_PCs", "seedtabs", paste0("nl_", num_PC, "_70_", NNseed, ".csv")
  )
)

X <- load_X(X_dir, n_PCs = num_PC, n_h = 70, seed = NNseed)

if(num_PC == "9") {
  load(file.path("2_application", "21_all_data", "nl_settings.Rdata"))
} else {
  load(file.path("2_application", "3_PCs", paste0(num_PC, "_settings.Rdata")))
}
prob_lambdas <- prob_lambda
mean_lambdas <- mean_lambda
cv_gridsize <- 1

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

    cur_X <- load_X(X_dir, n_PCs = num_PC, n_h = 70, seed = NNseed, ofold = ifold)

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
