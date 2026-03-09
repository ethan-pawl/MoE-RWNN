library(flowmix)

args <- commandArgs(trailingOnly = TRUE)
cv_step <- args[1]
arraynum <- as.integer(args[2])
stopifnot(cv_step %in% c("maxres", "cv", "refit", "summary"))

# Load the response data
imean <- 2
iprob <- 1
iint <- 5
dataSeed <- 1

simdata <- readRDS(
  file.path(
      "1_simulation", 
      "simdata",
      paste0(
          "simdata", 
          "_imean_", 
          imean, 
          "_iprob_", 
          iprob, 
          "_iint_", 
          iint,
          "_dataSeed_", 
          dataSeed, 
          ".RDS"
      )
  )
)

# Data
ylist <- simdata$ybin_list
countslist <- simdata$countslist

# Estimation settings
maxdev <- diff(range(simdata$mean_spec)) / 2
numclust <- 2

# Cross-validation settings
cv_gridsize <- 10
nfold <- 5
blocksize <- 20
nrep <- 10
max_prob_lambda <- 24
max_mean_lambda <- 24

##################################

# Model Cross-Validation

# The one_job function expects the data point indices 
# from the original data set (not subsetted, like we are doing here), 
# so we need to use the rownames to store the original indices.
names(ylist) <- 1:length(ylist)
names(countslist) <- 1:length(countslist)

jobs <- list(
  n_PCs = 9, 
  n_h = 70, 
  seed = 1:5
)

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
if(cv_step == "maxres") { # 1-4
  # model <- models[arraynum]
  job_grid <- expand.grid(jobs)
  job <- job_grid[arraynum,]

  destin <- file.path("3_seedsim",  
                      "results", 
                      job$seed)
                      
  if(!dir.exists(destin)) dir.create(destin, recursive = TRUE)

  # n_h <- if(model == "linear") "NA" else "70"
  # seed <- if(model == "linear") "NA" else "1"
  X <- load_X(X_dir, n_PCs = job$n_PCs, n_h = job$n_h, seed = job$seed)

  # The one_job function expects the data point indices 
  # from the original data set (not subsetted, like we are doing here), 
  # so we need to use the rownames to store the original indices.
  rownames(X) <- 1:nrow(X) 

  maxres <- get_max_lambda(destin,
                          maxres_file = "maxres.Rdata",
                          ylist = ylist,
                          countslist = countslist,
                          X = X,
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
} else if(cv_step == "cv") {
  jobs$ialpha <- jobs$ibeta <- 1:cv_gridsize
  jobs$ifold <- 1:nfold
  jobs$irep <- 1:nrep


  job_grid <- expand.grid(jobs)

  folds <- make_cv_folds(ylist, nfold, blocksize)

  # The below for loop assigns multiple cross-validation jobs 
  # to each SLURM array index. This is useful if you only have 
  # access to a limited number of CPUs
  for(j in -4:0 + 5 * arraynum) { # 1-5000 SLURM array --> 1-25000 job_grid rows
    job <- job_grid[j,]

    X <- load_X(X_dir, n_PCs = job$n_PCs, n_h = job$n_h, seed = job$seed, ofold = job$ifold)
    
    # The one_job function expects the data point indices 
    # from the original data set (not subsetted, like we are doing here), 
    # so we need to use the rownames to store the original indices.
    rownames(X) <- 1:nrow(X) 

    destin <- file.path("3_seedsim",  
                      "results", 
                      job$seed)

    load(file.path(destin, "prob_lambdas.RData"))
    load(file.path(destin, "mean_lambdas.RData"))

    seedfile <- file.path("3_seedsim", 
                          "seedtabs", 
                          paste0(job$seed, "_seedtab.csv"))

    seedtab <- read.csv(seedfile)

    # save meta file
    save(nfold,
         nrep,
         cv_gridsize,
         mean_lambdas,
         prob_lambdas,
         file = file.path(destin, 'meta.Rdata'))

    one_job(job$ialpha, 
            job$ibeta, 
            job$ifold, 
            job$irep, 
            folds,
            destin, 
            mean_lambdas, 
            prob_lambdas, 
            numclust = numclust,
            maxdev = maxdev,
            sim = FALSE, 
            seedtab = seedtab, 
            ylist = ylist, 
            countslist = countslist, 
            X = X)
  }
} else if(cv_step == "refit") { # 1-500
  jobs$ialpha <- jobs$ibeta <- 1:cv_gridsize
  job_grid <- expand.grid(jobs)
  job <- job_grid[arraynum,]

  destin <- file.path("3_seedsim",  
                      "results", 
                      job$seed)

  load(file.path(destin, "prob_lambdas.RData"))
  load(file.path(destin, "mean_lambdas.RData"))

  X <- load_X(X_dir, n_PCs = job$n_PCs, n_h = job$n_h, seed = job$seed)

  # The one_job function expects the data point indices 
  # from the original data set (not subsetted, like we are doing here), 
  # so we need to use the rownames to store the original indices.
  rownames(X) <- 1:nrow(X) 

  seedfile <- file.path("3_seedsim", 
                        "seedtabs", 
                        paste0(job$seed, "_seedtab.csv"))
  seedtab <- read.csv(seedfile)

  one_job_refit(job$ialpha,
                job$ibeta, 
                destin, 
                mean_lambdas, 
                prob_lambdas,
	              nrep = nrep,	
                numclust = numclust,
                maxdev = maxdev,
                seedtab = seedtab, 
                ylist = ylist, 
                countslist = countslist, 
                X = X)
} else if(cv_step == "summary") { # 1-5
  job_grid <- expand.grid(jobs)
  job <- job_grid[arraynum,]

  destin <- file.path("3_seedsim",  
                      "results", 
                      job$seed)
  
  cv_summary(destin = destin,
             save = TRUE,
             filename = paste0(job$seed, "_summary", ".RDS"))
}
