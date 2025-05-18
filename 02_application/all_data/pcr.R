# Requires MGL1704-hourly-paper.RDS in the project data folder
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
datobj <- readRDS(file = "../../data/paper-data-v2/MGL1704-hourly-paper.RDS")
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()

p <- ncol(X)

# Estimation settings
n.h <- 70
numclust <- 10

# Cross-validation settings
cv_gridsize <- 10
nfold <- 5
blocksize <- 20
nrep <- 10
maxdev <- 0.5
max_mean_lambda <- 40
max_prob_lambda <- 2 

# ELM transformation functions

make_W_unif <- function(p, n.h, a) {
    return(matrix(runif((p + 1) * n.h, - a, a), nrow = p + 1))
}

activate <- function(V) {
    return(1 / (1 + exp(-V)))
}

X_hidden <- function(X, p, n.h, a) {
  W <- make_W_unif(p, n.h, a)
  return(activate(cbind(1, X) %*% W))
}

###############################

#  Project X into the principal components space
X_no_cp <- X[,3:ncol(X)]

# X has already been centered and scaled
X_pca <- prcomp(X_no_cp, center = FALSE, scale = FALSE)
pca_var <- X_pca$sdev**2
cpve <- cumsum(pca_var / sum(pca_var)) 

d <- which(cpve > 0.95)[1]
# The first d principal components explain most of the variation

X_pc <- X_no_cp %*% X_pca$rotation[,1:d]

##################################

# The one_job function expects the data point indices 
# from the original data set (not subsetted, like we are doing here), 
# so we need to use the rownames to store the original indices.
rownames(X_pc) <- 1:nrow(X_pc) 
names(ylist) <- 1:length(ylist)
names(countslist) <- 1:length(countslist)

jobs <- list()

models <- c("linear", "half")
jobs$model <- models

# Perform the selected cross-validation step
if(cv_step == "maxres") {
  model <- models[arraynum]

  destin <- paste0(proj_name, "_", model)
  if(!dir.exists(destin)) {
    dir.create(destin)
  }

  maxres_file <- paste0(proj_name, "_maxres_", model,
                        ".Rdata")

  if(model == "half") {
    set.seed(0)
    use_X <- X_hidden(X_pc, d, n.h, 0.5)
  } else {
    use_X <- X_pc
  }

  maxres <- get_max_lambda(destin,
                          maxres_file = maxres_file,
                          ylist = ylist,
                          countslist = countslist,
                          X = use_X,
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

  orig_X <- X_pc # this ensures proper naming of "X" in the meta file, 
              #   without overwriting original X 

  for(j in -1:0 + 2 * arraynum) { # 1-5000 SLURM array --> 1-10000 job_grid rows
    job <- job_grid[j,]
    destin <- paste0(proj_name, "_", job$model)
    load(file.path(destin, "prob_lambdas.RData"))
    load(file.path(destin, "mean_lambdas.RData"))

    if(job$model == "half") {
      set.seed(0)
      X <- X_hidden(orig_X, d, n.h, 0.5)
    } else {
      X <- orig_X
    }

    seedfile <- paste0(proj_name, "_seedtab_", job$model, ".csv")
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

  destin <- paste0(proj_name, "_", job$model)
  load(file.path(destin, "prob_lambdas.RData"))
  load(file.path(destin, "mean_lambdas.RData"))

  if(job$model == "half") {
    set.seed(0)
    use_X <- X_hidden(X_pc, d, n.h, 0.5)
  } else {
    use_X <- X_pc
  }

  seedfile <- paste0(proj_name, "_seedtab_", job$model, ".csv")
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
                X = use_X)
} else if(cv_step == "summary") {
  model <- models[arraynum]

  destin <- paste0(proj_name, "_", model)
  
  cv_summary(destin = destin,
             save = TRUE,
             filename = paste0(proj_name, "_summary_", model, ".RDS"))
}
