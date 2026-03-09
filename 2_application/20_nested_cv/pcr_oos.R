library(flowmix) 
library(magrittr)

# Take cross-validation settings from SLURM index
args <- commandArgs(trailingOnly = TRUE)
cv_step <- args[1]
arraynum <- as.integer(args[2])
stopifnot(cv_step %in% c("maxres", "cv", "refit", "summary"))

# Load the data
datobj <- readRDS(file.path("data", "MGL1704-hourly-paper.RDS"))
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()

# Estimation settings
n.h <- 70
numclust <- 10
maxdev <- 0.5

# Cross-validation settings
cv_gridsize <- 10
nfold <- 5
blocksize <- 20
nrep <- 10
max_mean_lambda <- 40
max_prob_lambda <- 2 

##################################

# Model Cross-Validation

# The one_job function expects the data point indices 
# from the original data set (not subsetted, like we are doing here), 
# so we need to use the rownames to store the original indices.
names(ylist) <- 1:length(ylist)
names(countslist) <- 1:length(countslist)

jobs <- list()
jobs$n_PCs <- 9
jobs$nh <- c("NA", "70")
jobs$ofold <- 1:5

simple_jobs <- expand.grid(jobs)
simple_jobs$seed <- ifelse(simple_jobs$nh == "NA", "NA", "1")

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

# Perform the selected cross-validation step
if(cv_step == "maxres") {
  job <- simple_jobs[arraynum,]

  # Folder to save results
  destin <- file.path("2_application", 
                      "20_nested_cv", 
                      "results", 
                      paste0("nh_", job$nh, "_seed_", job$seed, "_ofold_", job$ofold))

  if(!dir.exists(destin)) dir.create(destin, recursive = TRUE)
  
  # get_max_lambda with PC learned on all the data except the outer fold
  X <- load_X(X_dir, n_PCs = job$n_PCs, n_h = job$nh, seed = job$seed, job$ofold)

  # The one_job function expects the data point indices 
  # from the original data set (not subsetted, like we are doing here), 
  # so we need to use the rownames to store the original indices.
  rownames(X) <- 1:nrow(X) 

  insample_inds <- unlist(ofolds[-job$ofold])

  maxres <- get_max_lambda(destin,
                          maxres_file = "maxres.Rdata",
                          ylist = ylist[insample_inds],
                          countslist = countslist[insample_inds],
                          X = X[insample_inds,],
                          numclust = numclust,
                          maxdev = maxdev,
                          max_mean_lambda = max_mean_lambda,
                          max_prob_lambda = max_prob_lambda)
    
  prob_lambdas <- logspace(min = 0.0001, max = maxres$alpha, 
                          length = cv_gridsize)
  mean_lambdas <-logspace(min = 0.0001, max = maxres$beta, 
                            length = cv_gridsize)
                            
  save(prob_lambdas, file = file.path(destin, "prob_lambdas.RData"))
  save(mean_lambdas, file = file.path(destin, "mean_lambdas.RData"))

  save(nfold,
       nrep,
       cv_gridsize,
       mean_lambdas,
       prob_lambdas,
       file = file.path(destin, 'meta.Rdata'))

} else if(cv_step == "cv") {
  jobs$ialpha <- jobs$ibeta <- 1:cv_gridsize
  jobs$ifold <- 1:nfold
  jobs$irep <- 1:nrep
  job_grid <- expand.grid(jobs)
  job_grid$seed <- ifelse(job_grid$nh == "NA", "NA", "1")

  orig_ylist <- ylist
  orig_countslist <- countslist
  # The below for loop assigns multiple cross-validation jobs 
  # to each SLURM array index. This is useful if you only have 
  # access to a limited number of CPUs
  for(j in -9:0 + 10 * arraynum) { 
    job <- job_grid[j,]

    X <- load_X(X_dir, n_PCs = 9, n_h = job$nh, seed = job$seed, ofold = job$ofold, ifold = job$ifold)

    # The one_job function expects the data point indices 
    # from the original data set (not subsetted, like we are doing here), 
    # so we need to use the rownames to store the original indices.
    rownames(X) <- 1:nrow(X) 

    destin <- file.path("2_application", 
                      "20_nested_cv", 
                      "results", 
                      paste0("nh_", job$nh, "_seed_", job$seed, "_ofold_", job$ofold))

    load(file.path(destin, "prob_lambdas.RData"))
    load(file.path(destin, "mean_lambdas.RData"))

    seedfile <- file.path("2_application", 
                          "20_nested_cv", 
                          "seedtabs", 
                          paste0("seedtab_", job$nh, "_", job$ofold, ".csv"))
    seedtab <- read.csv(seedfile)

    # Subsample the response data based on indices from the original data
    insample_inds <- unlist(ifolds[[job$ofold]])
    ylist <- orig_ylist[insample_inds]
    countslist <- orig_countslist[insample_inds]

    # The previous indexing reshuffles the data, so we need to reshuffle X as well
    # use the indices mapped to the subsampled dataset
    folds <- ifolds_inner_inds[[job$ofold]]
    X  <- X[unlist(folds),]

    # The previous indexing steps reshuffle the data into blocked cross-validation folds,
    # so we can now just index sequentially
    ends <- cumsum(sapply(folds, length))
    starts <- c(1, head(ends, -1) + 1)
    folds <- Map(":", starts, ends)

    one_job(ialpha = job$ialpha, 
            ibeta = job$ibeta, 
            ifold = job$ifold, 
            irep = job$irep, 
            folds = folds, 
            destin = destin, 
            mean_lambdas = mean_lambdas, 
            prob_lambdas = prob_lambdas, 
            numclust = numclust,
            maxdev = maxdev,
            sim = FALSE, 
            seedtab = seedtab, 
            ylist = ylist, 
            countslist = countslist, 
            X = X,  
            verbose = TRUE)
  }
} else if(cv_step == "refit") {
  jobs$ialpha <- jobs$ibeta <- 1:cv_gridsize
  job_grid <- expand.grid(jobs)
  job_grid$seed <- ifelse(job_grid$nh == "NA", "NA", "1")
  
  job <- job_grid[arraynum,]

  # Either way, loads in an object called X
  X <- load_X(X_dir, n_PCs = 9, n_h = job$nh, seed = job$seed, ofold = job$ofold)

  # The one_job function expects the data point indices 
  # from the original data set (not subsetted, like we are doing here), 
  # so we need to use the rownames to store the original indices.
  rownames(X) <- 1:nrow(X) 

  destin <- file.path("2_application", 
                      "20_nested_cv", 
                      "results", 
                      paste0("nh_", job$nh, "_seed_", job$seed, "_ofold_", job$ofold))

  load(file.path(destin, "prob_lambdas.RData"))
  load(file.path(destin, "mean_lambdas.RData"))

  seedfile <- file.path("2_application", 
                        "20_nested_cv", 
                        "seedtabs", 
                        paste0("seedtab_", job$nh, "_", job$ofold, ".csv"))
  seedtab <- read.csv(seedfile)

  insample_inds <- unlist(ofolds[-job$ofold])

  one_job_refit(ialpha = job$ialpha,
                ibeta = job$ibeta, 
                destin = destin, 
                mean_lambdas = mean_lambdas, 
                prob_lambdas = prob_lambdas, 
                nrep = nrep, 
                numclust = numclust,
                maxdev = maxdev,
                seedtab = seedtab, 
                ylist = ylist[insample_inds], 
                countslist = countslist[insample_inds], 
                X = X[insample_inds,], 
                verbose = TRUE)
} else if(cv_step == "summary") {
  job <- simple_jobs[arraynum,]

  destin <- file.path("2_application", 
                      "20_nested_cv", 
                      "results", 
                      paste0("nh_", job$nh, "_seed_", job$seed, "_ofold_", job$ofold))

  cv_summary(destin = destin,
             save = TRUE,
             filename = paste0("nh_", job$nh, "_seed_", job$seed, "_ofold_", job$ofold, ".RDS"))
}