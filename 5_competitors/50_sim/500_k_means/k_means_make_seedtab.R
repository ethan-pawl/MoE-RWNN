library(dplyr)
RNGkind("L'Ecuyer-CMRG")

# Cross-validation settings
cv_gridsize <- 1
nfold <- 296
nrep <- 10

nrows <- cv_gridsize * nfold * nrep

## Generate all the seeds
seedfile <- file.path("5_competitors", "50_sim", "500_k_means", "seedtab.csv")

## Make sure the seeds have not already been generated.
if(file.exists(seedfile)) {
  stop("seedtab already exists.")
}

## Generate the random number states (7 integers each)
set.seed(NULL)
s <- list(.Random.seed)
for(ii in 2:nrows){
  s[[ii]] <- parallel::nextRNGStream(s[[ii-1]])
}

s <- s %>% do.call(rbind, .)
colnames(s) <- paste0("seed", 1:7)
s <- s %>% as_tibble()

# ifold == 0 corresponds to the refit step
tab <- expand.grid(
  ialpha = 1:cv_gridsize, 
  ibeta = 1:cv_gridsize, 
  ifold = 1:nfold, 
  irep = 1:nrep
) %>% as_tibble()

seedtab <- tab %>% bind_cols(s)

## Write the table to file
write.csv(seedtab, file = seedfile, row.names = FALSE)
print(paste0("Wrote table containing seeds to ", seedfile))
