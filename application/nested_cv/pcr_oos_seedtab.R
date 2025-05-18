proj_name <- "pcr_oos"

library(magrittr)
library(tibble)
library(dplyr)
RNGkind("L'Ecuyer-CMRG")

# CV settings
nalpha <- 10
nbeta <- 10
nfold <- 5
nrep <- 10

nrows <- nalpha * nbeta * (nfold + 1) * nrep

seeddir <- "seedtabs"
if(!dir.exists(seeddir)) dir.create(seeddir)
## Generate all the seeds

## Form destination folder
## and make file names for saving the seeds
## Form destination folder
for(model in c("linear", "half")) {
    for(outer_fold in 1:nfold) {
        seedfile <- paste0(proj_name, "_seedtab_", model, "_", outer_fold, ".csv")

        ## Make sure the seeds have not already been generated.
        if(file.exists(seedfile)) {
            stop("seedtab already exists.")
        }

        ## Generate the random number /states/ (7 integers each)
        set.seed(NULL)
        s <- list(.Random.seed)
        for (ii in 2:nrows){
            s[[ii]] <- parallel::nextRNGStream(s[[ii-1]])
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
        write.csv(seedtab, file = file.path(seeddir, seedfile), row.names = FALSE)
        print(paste0("Wrote table containing seeds to", seedfile))
    }
}
