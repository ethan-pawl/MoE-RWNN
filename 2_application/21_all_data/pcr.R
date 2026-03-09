proj_name <- "pcr"

library(flowmix) 
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
cv_step <- args[1]
arraynum <- as.integer(args[2])
stopifnot(cv_step %in% c("maxres", "cv", "refit", "summary"))

# Load the data
datobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper.RDS"))
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()

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
names(ylist) <- 1:length(ylist)
names(countslist) <- 1:length(countslist)

jobs_nl_PCs <- list(
  model = "nl", 
  n_PCs = c(9, 18, 27, 37),
  n_h = 70,
  seed = 1
) 

jobs_nl_seed <- list(
  model = "nl", 
  n_PCs = 9,
  n_h = 70,
  seed = 2:5
) 

jobs_linear <- list(
  model = "linear",
  n_PCs = 9,
  n_h = "NA", 
  seed = "NA"
)

jobs_nl_no_wts <- list(
  model = "nl_no_wts",
  n_PCs = 9,
  n_h = 70, 
  seed = 1
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
if(cv_step == "maxres") { # 1-10
  # model <- models[arraynum]
  job_grid <- rbind(
    expand.grid(jobs_nl_PCs, stringsAsFactors = FALSE), 
    expand.grid(jobs_nl_seed, stringsAsFactors = FALSE),
    expand.grid(jobs_linear, stringsAsFactors = FALSE), 
    expand.grid(jobs_nl_no_wts, stringsAsFactors = FALSE)
  )

  job <- job_grid[arraynum,]

  destin <- file.path("2_application", 
                      "21_all_data", 
                      "results", 
                       paste(job$model, job$n_PCs, job$n_h, job$seed, sep = "_"))
                      
  if(!dir.exists(destin)) dir.create(destin, recursive = TRUE)

  # n_h <- if(model == "linear") "NA" else "70"
  # seed <- if(model == "linear") "NA" else "1"
  X <- load_X(X_dir, n_PCs = job$n_PCs, n_h = job$n_h, seed = job$seed)

  # The one_job function expects the data point indices 
  # from the original data set (not subsetted, like we are doing here), 
  # so we need to use the rownames to store the original indices.
  rownames(X) <- 1:nrow(X) 

  # Remove wts if not using pseudolikelihood (using likelihood-based model)
  countslist <- if(job$model == "nl_no_wts") NULL else countslist

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
  jobs_nl_PCs$ialpha <- jobs_nl_PCs$ibeta <- 1:cv_gridsize
  jobs_nl_PCs$ifold <- 1:nfold
  jobs_nl_PCs$irep <- 1:nrep

  jobs_nl_seed$ialpha <- jobs_nl_seed$ibeta <- 1:cv_gridsize
  jobs_nl_seed$ifold <- 1:nfold
  jobs_nl_seed$irep <- 1:nrep

  jobs_linear$ialpha <- jobs_linear$ibeta <- 1:cv_gridsize
  jobs_linear$ifold <- 1:nfold
  jobs_linear$irep <- 1:nrep

  jobs_nl_no_wts$ialpha <- jobs_nl_no_wts$ibeta <- 1:cv_gridsize
  jobs_nl_no_wts$ifold <- 1:nfold
  jobs_nl_no_wts$irep <- 1:nrep

  job_grid <- rbind(
    expand.grid(jobs_nl_PCs, stringsAsFactors = FALSE), 
    expand.grid(jobs_nl_seed, stringsAsFactors = FALSE),
    expand.grid(jobs_linear, stringsAsFactors = FALSE), 
    expand.grid(jobs_nl_no_wts, stringsAsFactors = FALSE)
  )

  folds <- make_cv_folds(ylist, nfold, blocksize)

  # The below for loop assigns multiple cross-validation jobs 
  # to each SLURM array index. This is useful if you only have 
  # access to a limited number of CPUs
  for(j in -9:0 + 10 * arraynum) { # 1-5000 SLURM array --> 1-50000 job_grid rows
    job <- job_grid[j,]

    # n_h <- if(job$model == "linear") "NA" else "70"
    # seed <- if(job$model == "linear") "NA" else "1"
    X <- load_X(X_dir, n_PCs = job$n_PCs, n_h = job$n_h, seed = job$seed, ofold = job$ifold)
    
    # The one_job function expects the data point indices 
    # from the original data set (not subsetted, like we are doing here), 
    # so we need to use the rownames to store the original indices.
    rownames(X) <- 1:nrow(X) 

    destin <- file.path("2_application", 
                      "21_all_data", 
                      "results", 
                       paste(job$model, job$n_PCs, job$n_h, job$seed, sep = "_"))

    load(file.path(destin, "prob_lambdas.RData"))
    load(file.path(destin, "mean_lambdas.RData"))

    seedfile <- file.path(
        "2_application", 
        "21_all_data", 
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

    # Remove wts if not using pseudolikelihood (using likelihood-based model)
    countslist <- if(job$model == "nl_no_wts") NULL else countslist

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
            X = X,
            verbose = TRUE)
  }
} else if(cv_step == "refit") { # 1-1000
  jobs_nl_PCs$ialpha <- jobs_nl_PCs$ibeta <- 1:cv_gridsize
  jobs_nl_seed$ialpha <- jobs_nl_seed$ibeta <- 1:cv_gridsize
  jobs_linear$ialpha <- jobs_linear$ibeta <- 1:cv_gridsize
  jobs_nl_no_wts$ialpha <- jobs_nl_no_wts$ibeta <- 1:cv_gridsize

  job_grid <- rbind(
    expand.grid(jobs_nl_PCs, stringsAsFactors = FALSE), 
    expand.grid(jobs_nl_seed, stringsAsFactors = FALSE),
    expand.grid(jobs_linear, stringsAsFactors = FALSE), 
    expand.grid(jobs_nl_no_wts, stringsAsFactors = FALSE)
  )

  # jobs$ialpha <- jobs$ibeta <- 1:cv_gridsize
  # job_grid <- expand.grid(jobs)
  job <- job_grid[arraynum,]

  destin <- file.path("2_application", 
                      "21_all_data", 
                      "results", 
                       paste(job$model, job$n_PCs, job$n_h, job$seed, sep = "_"))

  load(file.path(destin, "prob_lambdas.RData"))
  load(file.path(destin, "mean_lambdas.RData"))

  # n_h <- if(job$model == "linear") "NA" else "70"
  # seed <- if(job$model == "linear") "NA" else "1"
  X <- load_X(X_dir, n_PCs = job$n_PCs, n_h = job$n_h, seed = job$seed)

  # The one_job function expects the data point indices 
  # from the original data set (not subsetted, like we are doing here), 
  # so we need to use the rownames to store the original indices.
  rownames(X) <- 1:nrow(X) 

  seedfile <- file.path(
      "2_application", 
      "21_all_data", 
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

  # Remove wts if not using pseudolikelihood (using likelihood-based model)
  countslist <- if(job$model == "nl_no_wts") NULL else countslist

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
                X = X,
                verbose = TRUE)
} else if(cv_step == "summary") { # 1-10
  job_grid <- rbind(
    expand.grid(jobs_nl_PCs, stringsAsFactors = FALSE), 
    expand.grid(jobs_nl_seed, stringsAsFactors = FALSE),
    expand.grid(jobs_linear, stringsAsFactors = FALSE), 
    expand.grid(jobs_nl_no_wts, stringsAsFactors = FALSE)
  )

  job <- job_grid[arraynum,]
  # model <- models[arraynum]

  destin <- file.path("2_application", 
                      "21_all_data", 
                      "results", 
                       paste(job$model, job$n_PCs, job$n_h, job$seed, sep = "_"))
  
  cv_summary(destin = destin,
             save = TRUE,
             filename = paste0(paste(job$model, job$n_PCs, job$n_h, job$seed, "pcr", "summary", sep = "_"), ".RDS"))
}
