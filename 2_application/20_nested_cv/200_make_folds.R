library(flowmix)
library(magrittr)

datobj <- readRDS(file.path("data", "MGL1704-hourly-paper.RDS"))
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()

# Cross-validation settings
nfold = 5
blocksize = 20

# Make outer cross-validation folds
out_sample_folds = make_cv_folds(ylist, nfold = nfold, blocksize = blocksize)

# Now folds[[i]] is the out-of-sample (test) dataset for the i-th
# step in the outer CV procedure

# Given a set of indices after an outer fold (test dataset) 
# has been removed, construct a set of k inner folds
make_inner_fold <- function(inds, nfold, blocksize) {
    nblocks <- ceiling(length(inds) / blocksize)
    blocks_per_fold <- nblocks / nfold 

    # Create a matrix whose columns index the 
    # blocks to select for each inner fold
    inner_inds <- sapply(1:nfold, function(k) {
        seq(k, by = nfold, length.out = blocks_per_fold)
    })

    # Extract the blocks and group them into inner folds
    folds_mat <- apply(inner_inds, 2, function(cur_inds) {
        sapply(cur_inds, function(i) {
            inds[(1 + blocksize * (i - 1)) : (blocksize * i)]
        }) %>% as.vector()
    })

    # Convert to list
    folds_list <- lapply(1:nfold, function(k) folds_mat[,k])

    if((nblocks / nfold) %% 1 == 0) {
        # If blocks_per_fold is a whole number, there will be out of bound 
        # indices to remove in the last inner fold
        y <- folds_list[[nfold]]
        y <- y[1:which(y == inds[length(inds)])]
        folds_list[[nfold]] <- y
    } else {
        # If blocks_per_fold is fractional, there will be some 
        # NAs to remove
        folds_list <- lapply(folds_list, function(x) x[!is.na(x)])
    }

    return(folds_list)
}

# Create all inner folds
in_sample_folds <- lapply(1:5, function(k) {
    outer_folds <- out_sample_folds[setdiff(1:5, k)] %>% 
        unlist()

    make_inner_fold(outer_folds, nfold, blocksize)
})

folds_dir <- file.path("2_application", "20_nested_cv", "folds")
if(!dir.exists(folds_dir)) dir.create(folds_dir)

save(in_sample_folds, file = file.path(folds_dir, "in_sample_folds.Rdata"))
save(out_sample_folds, file = file.path(folds_dir, "out_sample_folds.Rdata"))
