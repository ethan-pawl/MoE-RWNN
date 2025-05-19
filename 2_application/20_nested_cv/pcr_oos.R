proj_name <- "pcr_oos"

library(flowmix) 
library(tidyverse)
library(conflicted)

# Load training data
load(file.path("2_application", 
               "20_nested_cv", 
               "folds", 
               "in_sample_folds.RData"))

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

# Select in-sample data for the current "held out" dataset
subset_data <- function(in_sample_folds, fold_ind, X, ylist, countslist) {
  in_sample <- in_sample_folds[[fold_ind]]
  in_ind <- unlist(in_sample)

  res <- list()
  res$in_sample <- in_sample
  res$X <- X[in_ind,]
  res$ylist <- ylist[in_ind]
  res$countslist <- countslist[in_ind]
  
  return(res)
}

jobs <- list()
jobs$model <- c("linear", "nl")
jobs$outer_fold <- 1:5

simple_jobs <- expand.grid(jobs)

# Perform the selected cross-validation step
if(cv_step == "maxres") {
  job <- simple_jobs[arraynum,]

  # Either way, loads in an object called X
  load(file.path("data", 
                 switch(job$model, 
                        linear = "X_pc.Rdata", 
                        nl = "X_nl.Rdata")))
  
  # The one_job function expects the data point indices 
  # from the original data set (not subsetted, like we are doing here), 
  # so we need to use the rownames to store the original indices.
  rownames(X) <- 1:nrow(X) 

  # Folder to save results
  destin <- file.path("2_application", 
                      "20_nested_cv", 
                      "results", 
                      paste0(job$model, "_", job$outer_fold))
    
  if(!dir.exists(destin)) dir.create(destin, recursive = TRUE)

  data_in <- subset_data(in_sample_folds, job$outer_fold, X, ylist, 
                        countslist)

  maxres <- get_max_lambda(destin,
                          maxres_file = "maxres.Rdata",
                          ylist = data_in$ylist,
                          countslist = data_in$countslist,
                          X = data_in$X,
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

  orig_ylist <- ylist
  orig_countslist <- countslist
  # The below for loop assigns multiple cross-validation jobs 
  # to each SLURM array index. This is useful if you only have 
  # access to a limited number of CPUs
  for(j in -9:0 + 10 * arraynum) { 
    job <- job_grid[j,]

     # Either way, loads in an object called X
    load(file.path("data", 
                   switch(job$model, 
                          linear = "X_pc.Rdata", 
                          nl = "X_nl.Rdata")))

    # The one_job function expects the data point indices 
    # from the original data set (not subsetted, like we are doing here), 
    # so we need to use the rownames to store the original indices.
    rownames(X) <- 1:nrow(X) 

    destin <- file.path("2_application", 
                      "20_nested_cv", 
                      "results", 
                      paste0(job$model, "_", job$outer_fold))

    load(file.path(destin, "prob_lambdas.RData"))
    load(file.path(destin, "mean_lambdas.RData"))

    data_in <- subset_data(in_sample_folds, job$outer_fold, X, orig_ylist, 
                          orig_countslist)

    seedfile <- file.path("2_application", 
                          "20_nested_cv", 
                          "seedtabs", 
                          paste0(proj_name, "_seedtab_", job$model, "_", job$outer_fold, ".csv"))
    seedtab <- read.csv(seedfile)

    # name objects properly before saving meta file
    folds <- data_in$in_sample
    ylist <- data_in$ylist
    countslist <- data_in$countslist
    X <- data_in$X

    save(folds,
         nfold,
         nrep,
         cv_gridsize,
         mean_lambdas,
         prob_lambdas,
         ylist, 
         countslist, 
         X,
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
            X = X, 
            subsampled = TRUE, 
            verbose = TRUE)
  }
} else if(cv_step == "refit") {
  jobs$ialpha <- jobs$ibeta <- 1:cv_gridsize
  job_grid <- expand.grid(jobs)
  job <- job_grid[arraynum,]

  # Either way, loads in an object called X
  load(file.path("data", 
                 switch(job$model, 
                        linear = "X_pc.Rdata", 
                        nl = "X_nl.Rdata")))

  # The one_job function expects the data point indices 
  # from the original data set (not subsetted, like we are doing here), 
  # so we need to use the rownames to store the original indices.
  rownames(X) <- 1:nrow(X) 

  destin <- file.path("2_application", 
                      "20_nested_cv", 
                      "results", 
                      paste0(job$model, "_", job$outer_fold))

  load(file.path(destin, "prob_lambdas.RData"))
  load(file.path(destin, "mean_lambdas.RData"))

  data_in <- subset_data(in_sample_folds, job$outer_fold, X, ylist, 
                        countslist)

  seedfile <- file.path("2_application", 
                        "20_nested_cv", 
                        "seedtabs", 
                        paste0(proj_name, "_seedtab_", job$model, "_", job$outer_fold, ".csv"))
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
                ylist = data_in$ylist, 
                countslist = data_in$countslist, 
                X = data_in$X, 
                verbose = TRUE)
} else if(cv_step == "summary") {
  job <- simple_jobs[arraynum,]

  destin <- file.path("2_application", 
                      "20_nested_cv", 
                      "results", 
                      paste0(job$model, "_", job$outer_fold))

  cv_summary(destin = destin,
             save = TRUE,
             filename = paste0(proj_name, "_", job$model, "_", job$outer_fold, ".RDS"))
}