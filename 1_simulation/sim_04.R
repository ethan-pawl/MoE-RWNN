library(flowmix)
library(parallel)
library(parallelly)
library(dplyr)

# To prevent multi-threaded BLAS which conflicts with mclapply
RhpcBLASctl::blas_set_num_threads(1)

# Input SLURM array index
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

#################

# Make a table to map the SLURM array index to a simulation scenario 
sims <- rbind(expand.grid(1, 1:4, 1:10, c("l", "n"), stringsAsFactors = FALSE), 
              expand.grid(2:4, 1, 1:10, c("l", "n"), stringsAsFactors = FALSE))              
colnames(sims) <- c("imean", "iprob", "iint", "modelFit")

############################

# For each simulation scenario, make a table 
# where each row is a seed to randomly initialize 
# cluster means in the EM algorithm at each point in
#  the cross-validation procedure (for reproducibility).

#  If you want to reproduce our results, don't run this function. Instead, 
# use the given seedtabs in the GitHub repository.
make_seedtab <- function(imean, iprob, iint, modelFit) {
    RNGkind("L'Ecuyer-CMRG")
    # 10 x 10 5-fold cross-validation with 30 EM restarts
    nalpha <- 10
    nbeta <- 10
    nfold <- 5
    nrep <- 30
    nrows <- nalpha * nbeta * (nfold + 1) * nrep

    seed_destin <- file.path("1_simulation", "seedtabs")
    if(!dir.exists(seed_destin)) dir.create(seed_destin)

    seedfile <- file.path(seed_destin,
                          paste0(imean, "-", iprob, "-", iint, "-", modelFit, "_seedtab.csv"))
    
    if(file.exists(seedfile)) {
        print("seedtab already exists.")
    } else {
        ## Generate the random number /states/ (7 integers each)
        set.seed(NULL)
        s <- list(.Random.seed)
        for (ii in 2:nrows){
            s[[ii]] <- nextRNGStream(s[[ii-1]])
        }
        s = s %>% do.call(rbind, .) %>% as_tibble()
        colnames(s) = paste0("seed", 1:7)
        # ifold == 0 corresponds to the refit step
        tab = expand.grid(ialpha = 1:nalpha, 
                            ibeta = 1:nbeta, 
                            ifold = 0:nfold, 
                            irep = 1:nrep) %>% as_tibble()
        seedtab = tab %>% bind_cols(s)
        ## Write the table to file
        write.csv(seedtab, file = seedfile, row.names = FALSE)
        print(paste0("Wrote table containing seeds to", seedfile))
    }
    RNGkind("default") # Reset the RNG
}

# Returns TRUE if the model has already been fit to the 
# i-th simulation scenario 
sim_done <- function(i, sims) {
    # Simulation scenario settings
    imean <- sims[i, "imean"]
    iprob <- sims[i, "iprob"]
    iint <- sims[i, "iint"]
    modelFit <- sims[i, "modelFit"]

    fname <- file.path("1_simulation",
                       "results", 
                       paste0(imean, "-", iprob, "-", iint, "-", modelFit),
                       paste0(imean, "-", iprob, "-", iint, "-", modelFit, "_summary.RDS"))
    
    
    is_done <- file.exists(fname)
    if(is_done) {
        cat("Simulation ", imean, "-", iprob, "-", iint, "-", modelFit, " is done.\n", sep = "")
    }

    return(is_done)
}

# If you want to reproduce our results, set new_seedtab = FALSE.
# If you want to run on a new dataset or have some other reason 
# for wanting new results, set new_seedtab = TRUE
run_sim <- function(i, sims, new_seedtab = FALSE) {
    # Simulation scenario settings
    imean <- sims[i, "imean"]
    iprob <- sims[i, "iprob"]
    iint <- sims[i, "iint"]
    modelFit <- sims[i, "modelFit"]

    if(new_seedtab) {
        make_seedtab(imean, iprob, iint, modelFit)
    }

    cat("Running simulation ", imean, "-", iprob, "-", iint, "-", modelFit, "\n", sep = "")

    load(file.path("1_simulation", 
                   "simdata", 
                   paste0("simdata-", imean, "-", iprob, "-", iint, "-", "1", ".Rdata")))

    # Data
    ylist <- simdata$ybin_list
    countslist <- simdata$countslist

    # Estimation settings
    maxdev <- diff(range(simdata$mean_spec)) / 2
    numclust <- 2

    # K-fold cross-validation settings
    nfold <- 5 
    cv_gridsize <- 10
    blocksize <- 20

    # Number of EM restarts
    nrep <- 30

    seedtab <- read.csv(file.path("1_simulation", 
                                  "seedtabs",
                                  paste0(imean, "-", iprob, "-", iint, "-", modelFit, "_seedtab.csv")))

    destin <- file.path("1_simulation",
                        "results", 
                        paste0(imean, "-", iprob, "-", iint, "-", modelFit))

    if(!dir.exists(destin)) {
        dir.create(destin, recursive = TRUE)
    }

    # either way, loads in an object named X
    load(file.path("data",
                   switch(modelFit, 
                          l = "X_pc.Rdata", 
                          n = "X_nl.Rdata")))

    # Define the cross-validation hyperparameter grid
    max_prob_lambda <- 24
    max_mean_lambda <- 24

    maxres <- get_max_lambda(destin, 
                             "maxres.Rdata",
                             ylist, 
                             countslist, 
                             X, 
                             numclust, 
                             maxdev, 
                             max_prob_lambda,
                             max_mean_lambda)

    prob_lambdas <- logspace(min = 0.0001, max = maxres$alpha, 
                            length = cv_gridsize) 
    mean_lambdas <- logspace(min = 0.0001, max = maxres$beta, 
                            length = cv_gridsize)

    save(prob_lambdas, file = file.path(destin, "prob_lambdas.Rdata"))
    save(mean_lambdas, file = file.path(destin, "mean_lambdas.Rdata"))

    # Define cross-validation folds
    folds <- make_cv_folds(ylist, nfold, blocksize)
    
    # Perform k-fold cross-validation
    cv.flowmix(ylist, 
               countslist, 
               X, 
               destin, 
               mean_lambdas, 
               prob_lambdas,
               NULL, 
               maxdev, 
               numclust, 
               nfold, 
               nrep, 
               TRUE, 
               FALSE, 
               TRUE, 
               availableCores() - 2, 
               blocksize, 
               folds, 
               seedtab)

    # Refit model on entire dataset
    cv.flowmix(ylist, 
               countslist, 
               X, 
               destin, 
               mean_lambdas, 
               prob_lambdas,
               NULL, 
               maxdev, 
               numclust, 
               nfold, 
               nrep, 
               TRUE, 
               TRUE, 
               FALSE, 
               availableCores() - 2, 
               blocksize, 
               folds, 
               seedtab)

    # Summarize k-fold cross-validation and refitting results
    cv_summary(destin = destin, 
               save = TRUE, 
               filename = paste0(imean, "-", iprob, "-", iint, "-", modelFit, "_summary.RDS"))
}

# Call run_sim(i, sims, TRUE) if you don't want to reproduce our results 
if(!sim_done(i, sims)) {
    run_sim(i, sims)
}
