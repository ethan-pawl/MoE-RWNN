# IMPORTANT: this also requires my fork of flowmix (ethan-pawl/flowmix,
# branch summary_no_meta)

proj_name <- "pcr"

library(flowmix) 
library(tidyverse)
library(conflicted)

args <- commandArgs(trailingOnly = TRUE)
cv_step <- args[1]
arraynum <- as.integer(args[2])
stopifnot(cv_step %in% c("maxres", "cv", "refit", "summary"))

# Load the data
datobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper.RDS"))
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()

p <- ncol(X)

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

models <- c("linear", "nl")
jobs$model <- models

# Perform the selected cross-validation step
if(cv_step == "maxres") {
  model <- models[arraynum]

  # TODO: check this
  destin <- file.path("02_application", 
                      "all_data", 
                      "results", 
                       model)
                      
  if(!dir.exists(destin)) dir.create(destin, recursive = TRUE)

  # Either way, loads in an object called X
  load(file.path("data", 
                 switch(model, 
                        linear = "X_pc.Rdata", 
                        nl = "X_nl.Rdata")))

  # The one_job function expects the data point indices 
  # from the original data set (not subsetted, like we are doing here), 
  # so we need to use the rownames to store the original indices.
  rownames(X_pc) <- 1:nrow(X_pc) 

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

  for(j in -1:0 + 2 * arraynum) { # 1-5000 SLURM array --> 1-10000 job_grid rows
    job <- job_grid[j,]

    # Either way, loads in an object called X
    load(file.path("data", 
                   switch(job$model, 
                          linear = "X_pc.Rdata", 
                          nl = "X_nl.Rdata")))

    destin <- file.path("02_application", 
                        "all_data", 
                        "results", 
                         job$model)

    load(file.path(destin, "prob_lambdas.RData"))
    load(file.path(destin, "mean_lambdas.RData"))

    seedfile <- file.path("02_application", "all_data", "seedtabs", paste0(proj_name, "_seedtab_", job$model, ".csv"))
    seedtab <- read.csv(seedfile)

    # save meta file
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
            X = X)
  }
} else if(cv_step == "refit") {
  jobs$ialpha <- jobs$ibeta <- 1:cv_gridsize
  job_grid <- expand.grid(jobs)
  job <- job_grid[arraynum,]

  destin <- file.path("02_application", 
                      "all_data", 
                      "results", 
                      job$model)

  load(file.path(destin, "prob_lambdas.RData"))
  load(file.path(destin, "mean_lambdas.RData"))

  # Either way, loads in an object called X
  load(file.path("data", 
                 switch(job$model, 
                        linear = "X_pc.Rdata", 
                        nl = "X_nl.Rdata")))

  seedfile <- file.path("02_application", "all_data", "seedtabs", paste0(proj_name, "_seedtab_", job$model, ".csv"))
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
} else if(cv_step == "summary") {
  model <- models[arraynum]

  destin <- file.path("02_application", 
                      "all_data", 
                      "results", 
                      model)
  
  cv_summary(destin = destin,
             save = TRUE,
             filename = paste0(proj_name, "_summary_", model, ".RDS"))
}
