proj_name <- "pcr"

library(magrittr)
library(tibble)
library(dplyr)
RNGkind("L'Ecuyer-CMRG")

# Cross-validation settings
nalpha <- 10
nbeta <- 10
nfold <- 5
nrep <- 10

nrows <- nalpha * nbeta * (nfold + 1) * nrep

## Generate all the seeds

## Form destination folder
## and make file names for saving the seeds
## Form destination folder

# DON'T run the next for loop if you want to reproduce our results,
# because it will overwrite the preexisting seedtabs
for(model in c("linear", "nl")) {
    seedfile <- file.path("seedtabs", paste0(proj_name, "_seedtab_", model, ".csv"))

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
    write.csv(seedtab, file = seedfile, row.names = FALSE)
    print(paste0("Wrote table containing seeds to", seedfile))
}
