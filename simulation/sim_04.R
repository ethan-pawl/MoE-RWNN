library(flowmix)
library(parallel)
library(parallelly)
library(dplyr)

RhpcBLASctl::blas_set_num_threads(1)

args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

#################

sims <- rbind(expand.grid(1, 1:4, 1:10, c("l", "n"), stringsAsFactors = FALSE), 
              expand.grid(2:4, 1, 1:10, c("l", "n"), stringsAsFactors = FALSE))              
colnames(sims) <- c("imean", "iprob", "iint", "modelFit")

############################

make_seedtab <- function(imean, iprob, iint, modelFit) {
    RNGkind("L'Ecuyer-CMRG")
    nalpha <- 10
    nbeta <- 10
    nfold <- 5
    nrep <- 30
    nrows <- nalpha * nbeta * (nfold + 1) * nrep

    seed_destin <- "seedtabs"

    if(!dir.exists(seed_destin)) {
        dir.create(seed_destin, recursive = TRUE)
    }

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
    RNGkind("default")
}

sim_done <- function(i, sims) {
    imean <- sims[i, "imean"]
    iprob <- sims[i, "iprob"]
    iint <- sims[i, "iint"]
    modelFit <- sims[i, "modelFit"]

    fname <- file.path("results", 
                        paste0(imean, "-", iprob, "-", iint, "-", modelFit),
                        paste0(imean, "-", iprob, "-", iint, "-", modelFit, "_summary.RDS"))
    
    
    is_done <- file.exists(fname)
    if(is_done) {
        cat("Simulation ", imean, "-", iprob, "-", iint, "-", modelFit, " is done.\n", sep = "")
    }

    return(is_done)
}

run_sim <- function(i, sims) {
    imean <- sims[i, "imean"]
    iprob <- sims[i, "iprob"]
    iint <- sims[i, "iint"]
    modelFit <- sims[i, "modelFit"]

    make_seedtab(imean, iprob, iint, modelFit)

    cat("Running simulation ", imean, "-", iprob, "-", iint, "-", modelFit, "\n", sep = "")

    load(file.path("~", 
                   "00_Cyto", 
                   "data", 
                   "simdata", 
                   "simdata_04_01",
                   paste0("simdata-", imean, "-", iprob, "-", iint, "-", "1", ".Rdata")))

    ylist <- res$ybin_list
    countslist <- res$countslist
    maxdev <- diff(range(res$mean_spec)) / 2

    numclust <- 2
    nfold <- 5 
    nrep <- 30
    cv_gridsize <- 10
    blocksize <- 20
    seedtab <- read.csv(file.path("seedtabs",
                                  paste0(imean, "-", iprob, "-", iint, "-", modelFit, "_seedtab.csv")))

    destin <- file.path("results", 
                        paste0(imean, "-", iprob, "-", iint, "-", modelFit))

    if(!dir.exists(destin)) {
        dir.create(destin, recursive = TRUE)
    }

    load(file.path("~", 
                   "00_Cyto", 
                   "data",
                   "X_data", 
                    switch(modelFit, 
                           l = "X_pc.Rdata", 
                           n = "X_nl_70.Rdata")))
    # either way, loads in an object named X

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

    folds <- make_cv_folds(ylist, nfold, blocksize)
    
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

    cv_summary(destin = destin, 
               save = TRUE, 
               filename = paste0(imean, "-", iprob, "-", iint, "-", modelFit, "_summary.RDS"))
}

if(!sim_done(i, sims)) {
    run_sim(i, sims)
}
