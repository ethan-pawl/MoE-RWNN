library(flowmix)
library(parallel)
library(parallelly)
library(dplyr)

# To prevent multi-threaded BLAS which conflicts with mclapply
RhpcBLASctl::blas_set_num_threads(1)

# Input SLURM array index
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

# Set to TRUE if you want new results
# Set to FALSE to reproduce our results
new_seedtab <- as.logical(args[2])

#################

# Make a table to map the SLURM array index to a simulation scenario 
# sims <- rbind(expand.grid(1, 1:4, 1:10, c("l", "n"), stringsAsFactors = FALSE), 
#               expand.grid(2:4, 1, 1:10, c("l", "n"), stringsAsFactors = FALSE))   
sims <- rbind(
    expand.grid(1, 1:4, 1:10, c("NA", "70"), stringsAsFactors = FALSE), 
    expand.grid(2:4, 1, 1:10, c("NA", "70"), stringsAsFactors = FALSE),
    expand.grid(2, 1, 5, c("35", "105", "140", "175"), stringsAsFactors = FALSE) # Robustness against hidden layer width
)
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
    nrep <- 10
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
    nrep <- 10

    seedtab <- read.csv(file.path("1_simulation", 
                                  "seedtabs",
                                  paste0(imean, "-", iprob, "-", iint, "-", modelFit, "_seedtab.csv")))

    destin <- file.path("1_simulation",
                        "results", 
                        paste0(imean, "-", iprob, "-", iint, "-", modelFit))

    if(!dir.exists(destin)) {
        dir.create(destin, recursive = TRUE)
    }

    X_dir <- file.path("data", "X_variations")
    X <- readRDS(file.path(X_dir, paste0("X_pc_9_nh_", modelFit, "_seed_NA_ofold_NA_ifold_NA.RDS")))

    # load(file.path("data",
    #                switch(modelFit, 
    #                       l = "X_pc.Rdata", 
    #                       n = "X_nl.Rdata")))

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

    # Make a grid to index individual CV jobs
    iimat <- make_iimat(
        cv_gridsize = cv_gridsize, 
        nfold = nfold, 
        nrep = nrep
    )

    # saving all metadata except X since X 
    # changes across folds
    save(
        folds,
        nfold,
        nrep,
        cv_gridsize,
        mean_lambdas,
        prob_lambdas,
        ylist, 
        countslist,
        file = file.path(destin, 'meta.Rdata')
    )

    print(paste0("wrote meta data to ", file.path(destin, 'meta.Rdata')))

    # Parallelize over all CV jobs
    empty <- mclapply(
        1:nrow(iimat), 
        function(ii) {
            ialpha <- iimat[,"ialpha"]
            ibeta <- iimat[,"ibeta"]
            ifold <- iimat[,"ifold"]
            irep <- iimat[,"irep"]

            cat("\r", ii "out of", nrow(iimat) "cross-validation jobs.")
            
            # Load the correct version of the 
            # dataset where PCA has been learned 
            # on the subset of the data with the 
            # current fold held out.
            cur_X <- readRDS(
                file.path(
                    X_dir, 
                    paste0(
                        "X_pc_9_nh_", 
                        modelFit, 
                        "_seed_NA_ofold_", 
                        ifold, 
                        "_ifold_NA.RDS"
                    )
                ) 
            )

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
        mc.cores = availableCores() - 1, 
        mc.preschedule = FALSE
    )

    # Now we are refitting on the 
    # entire dataset so it is OK to learn PCA 
    # from all the data.
    cv.flowmix(
        ylist = ylist, 
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
        mc.cores = availableCores() - 1, 
        blocksize = blocksize, 
        folds = folds, 
        seedtab = seedtab
    )

    # Summarize k-fold cross-validation and refitting results
    cv_summary(destin = destin, 
               save = TRUE, 
               filename = paste0(imean, "-", iprob, "-", iint, "-", modelFit, "_summary.RDS"))
}

# Call run_sim(i, sims, TRUE) if you don't want to reproduce our results 
if(!sim_done(i, sims)) {
    run_sim(i, sims, new_seedtab)
}
