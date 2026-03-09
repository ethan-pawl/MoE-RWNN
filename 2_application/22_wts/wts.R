library(flowmix)

args <- commandArgs(trailingOnly = TRUE)
cv_step <- args[1]
arraynum <- as.integer(args[2])
stopifnot(cv_step %in% c("maxres", "cv", "refit", "summary"))

# Load the data
datobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper-counts-and-biomass.RDS"))
datobj |> list2env(envir = .GlobalEnv) |> invisible()

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

##################################

# Model Cross-Validation

# The one_job function expects the data point indices 
# from the original data set (not subsetted, like we are doing here), 
# so we need to use the rownames to store the original indices.
names(ybin_list) <- 1:length(ybin_list)
names(countslist) <- 1:length(countslist)

jobs_nl_no_wts <- list(
  model = "nl_no_wts",
  n_PCs = 9,
  n_h = 70, 
  seed = 4
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
if(cv_step == "maxres") { # 1
  # model <- models[arraynum]
  job_grid <- expand.grid(jobs_nl_no_wts, stringsAsFactors = FALSE)

  job <- job_grid[arraynum,]

  destin <- file.path("2_application", 
                      "22_wts", 
                      "results", 
                       paste(job$model, job$n_PCs, job$n_h, job$seed, sep = "_"))
                      
  if(!dir.exists(destin)) dir.create(destin, recursive = TRUE)

  X <- load_X(X_dir, n_PCs = job$n_PCs, n_h = job$n_h, seed = job$seed)

  # The one_job function expects the data point indices 
  # from the original data set (not subsetted, like we are doing here), 
  # so we need to use the rownames to store the original indices.
  rownames(X) <- 1:nrow(X) 

  maxres <- get_max_lambda(destin,
                          maxres_file = "maxres.Rdata",
                          ylist = ybin_list,
                          countslist = countslist, # uses bin counts, not biomasses
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

  # save meta file
  save(
    nfold,
    nrep,
    cv_gridsize,
    mean_lambdas,
    prob_lambdas,
    file = file.path(destin, 'meta.Rdata')
  )
} else if(cv_step == "cv") {
  jobs_nl_no_wts$ialpha <- jobs_nl_no_wts$ibeta <- 1:cv_gridsize
  jobs_nl_no_wts$ifold <- 1:nfold
  jobs_nl_no_wts$irep <- 1:nrep

  job_grid <- expand.grid(jobs_nl_no_wts, stringsAsFactors = FALSE)

  folds <- make_cv_folds(ybin_list, nfold, blocksize)

  # The below for loop assigns multiple cross-validation jobs 
  # to each SLURM array index. This is useful if you only have 
  # access to a limited number of CPUs
  for(j in 1 * arraynum) { # 1-5000 SLURM array
    job <- job_grid[j,]

    X <- load_X(X_dir, n_PCs = job$n_PCs, n_h = job$n_h, seed = job$seed, ofold = job$ifold)
    
    # The one_job function expects the data point indices 
    # from the original data set (not subsetted, like we are doing here), 
    # so we need to use the rownames to store the original indices.
    rownames(X) <- 1:nrow(X) 

    destin <- file.path("2_application", 
                      "22_wts", 
                      "results", 
                       paste(job$model, job$n_PCs, job$n_h, job$seed, sep = "_"))

    load(file.path(destin, "prob_lambdas.RData"))
    load(file.path(destin, "mean_lambdas.RData"))

    seedfile <- file.path(
        "2_application", 
        "22_wts", 
        "seedtabs", 
        paste0(
            paste(
                job$model, 
                job$n_PCs, 
                job$n_h, 
                job$seed, 
                "seedtab",
                sep = "_"
            ), 
            ".csv"
        )
    )
    seedtab <- read.csv(seedfile)

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
            ylist = ybin_list, 
            countslist = countslist, 
            X = X,
            verbose = TRUE)
  }
} else if(cv_step == "refit") { # 1-100
  jobs_nl_no_wts$ialpha <- jobs_nl_no_wts$ibeta <- 1:cv_gridsize

  job_grid <- expand.grid(jobs_nl_no_wts, stringsAsFactors = FALSE)
  job <- job_grid[arraynum,]

  destin <- file.path("2_application", 
                      "22_wts", 
                      "results", 
                       paste(job$model, job$n_PCs, job$n_h, job$seed, sep = "_"))

  load(file.path(destin, "prob_lambdas.RData"))
  load(file.path(destin, "mean_lambdas.RData"))

  X <- load_X(X_dir, n_PCs = job$n_PCs, n_h = job$n_h, seed = job$seed)

  # The one_job function expects the data point indices 
  # from the original data set (not subsetted, like we are doing here), 
  # so we need to use the rownames to store the original indices.
  rownames(X) <- 1:nrow(X) 

  seedfile <- file.path(
      "2_application", 
      "22_wts", 
      "seedtabs", 
      paste0(
          paste(
              job$model, 
              job$n_PCs, 
              job$n_h, 
              job$seed, 
              "seedtab",
              sep = "_"
          ), 
          ".csv"
      )
  )
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
                ylist = ybin_list, 
                countslist = countslist, 
                X = X,
                verbose = TRUE)
} else if(cv_step == "summary") { # 1
  job_grid <- expand.grid(jobs_nl_no_wts, stringsAsFactors = FALSE)
  job <- job_grid[arraynum,]

  destin <- file.path("2_application", 
                      "22_wts", 
                      "results", 
                       paste(job$model, job$n_PCs, job$n_h, job$seed, sep = "_"))
  
  cv_summary(destin = destin,
             save = TRUE,
             filename = paste0(paste(job$model, job$n_PCs, job$n_h, job$seed, "wts", "summary", sep = "_"), ".RDS"))
}
