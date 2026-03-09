library(magrittr)
library(tibble)
library(dplyr)
RNGkind("L'Ecuyer-CMRG")

# Cross-validation settings
cv_gridsize <- 10
nfold <- 5
nrep <- 10

nrows <- cv_gridsize^2 * (nfold + 1) * nrep

## Generate all the seeds

## Form destination folder
## and make file names for saving the seeds
## Form destination folder
seeddir <- file.path("2_application", "23_PCs", "seedtabs")
if(!dir.exists(seeddir)) dir.create(seeddir)

jobs <- list(
  n_PCs = c(18, 27, 37), 
  seed = 2:5
)
job_grid <- expand.grid(jobs)

for(i in 1:nrow(job_grid)) {
  n_PCs <- job_grid[i,"n_PCs"]
  seed <- job_grid[i,"seed"]
  experiment <- paste0("nl_", n_PCs, "_70_", seed)
  seedfile <- file.path("2_application", "23_PCs", "seedtabs", paste0(experiment, ".csv"))

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

  s <- s %>% do.call(rbind, .)
  colnames(s) <- paste0("seed", 1:7)
  s <- s %>% as_tibble()
  # ifold == 0 corresponds to the refit step
  tab <- expand.grid(
    ialpha = 1:cv_gridsize, 
    ibeta = 1:cv_gridsize, 
    ifold = 0:nfold, 
    irep = 1:nrep
    ) %>% as_tibble()

    seedtab <- tab %>% bind_cols(s)

  ## Write the table to file
  write.csv(seedtab, file = seedfile, row.names = FALSE)
  print(paste0("Wrote table containing seeds to", seedfile))
}


