library(magrittr)
library(flowmix)

datobj <- readRDS(file.path("data", "MGL1704-hourly-paper.RDS"))
datobj %>% 
  list2env(envir = .GlobalEnv) %>% 
  invisible()

# Remove the changepoint variables
X_cont <- X[,3:39]

# Cross-validation and estimation settings
nfold <- 5
blocksize <- 20
a <- 0.5
TT <- length(ylist)
p <- ncol(X_cont)

# Create nested cross-validation folds and save to disk
ofolds <- make_cv_folds(ylist, nfold, blocksize)
ifolds_inner_inds <- list() # This object stores the indices with respect to 
# the outer folds. That is, the indices in ifolds_inner_inds[[i]] 
# index the entries in ofolds[[i]], NOT the original data points themselves.

# This object stores the indices with respect to the original dataset.
ifolds <- lapply(1:nfold, function(fold_ind) {
  train_inds <- unlist(ofolds[-fold_ind])
  ifold_mismatched_inds <- make_cv_folds(nfold = nfold, blocksize = blocksize, TT = length(train_inds))
  ifolds_inner_inds[[fold_ind]] <<- ifold_mismatched_inds

  lapply(ifold_mismatched_inds, function(cur_inds) {
    train_inds[cur_inds]
  })
})

save(
  ofolds, 
  ifolds, 
  ifolds_inner_inds, 
  file = file.path("data", "ofolds__ifolds__ifolds_inner_inds.Rdata")
)

# Create all PC- and RWNN- transformed datasets for each cross-validation split, 
# number of hidden nodes, number of PCs, etc. and save to disk

# This function ensures a common format for file names
save_X <- function(X, savedir, n_PCs = "NA", n_h = "NA", seed = "NA", ofold = "NA", ifold = "NA") {

  file_name <- paste("X_pc", n_PCs, "nh", n_h, "seed", seed, "ofold", ofold, "ifold", ifold, sep = "_")
  file_name <- paste0(file_name, ".RDS")
  file_name <- file.path(savedir, file_name)

  saveRDS(X, file_name)
}

X_dir <- file.path("data", "X_variations_v2")

if(!dir.exists(X_dir)) { # This is a bit intensive; don't repeat if already done
  
  dir.create(X_dir)

  # Test a range of number of PCs
  for(n_PCs in c(9, 18, 27, 37)) {

      # PC-transform the entire dataset (to refit the model after cross-validation)
      pr_out <- prcomp(X_cont, center = FALSE, scale = FALSE)
      X_pc <- X_cont %*% pr_out$rotation[,1:n_PCs]

      # To match the interpretation in 00_EDA.R (Figures 1, 10, and 11), 
      # we need to flip the signs of the PCs
      X_pc[,1] <- -X_pc[,1]
      # X_pc[,2] <- -X_pc[,2]
      # X_pc[,3] <- -X_pc[,3]
      # X_pc[,4] <- -X_pc[,4]

      save_X(X_pc, X_dir, n_PCs = n_PCs)

    # Test a range of number of hidden nodes
    n_h_vec <- if(n_PCs == 9) c(35, 70, 105, 140, 175) else 70
    for(n_h in n_h_vec) {
      # Resample the hidden layer weights several times
      seed_vec <- if(n_PCs == 9 & n_h == 70) 1:5 else 1
      for(seed in seed_vec) {

        # RWNN-transform the entire dataset (to refit the model after cross-validation)
        set.seed(seed)
        X_nl <- make_hidden_nodes(X_pc, n_h, a)
        save_X(X_nl, X_dir, n_PCs = n_PCs, n_h = n_h, seed = seed)

        # Single-layer cross-validation (for the simulation)
        for(ofold_ind in 1:nfold) {

          # Split the data
          cv_inds <- unlist(ofolds[-ofold_ind])
          test_inds <- ofolds[[ofold_ind]]

          X_cv <- X_cont[cv_inds,]
          X_test <- X_cont[test_inds,]

          # Calculate the PC rotation matrix using only training data
          cv_pr_out <- prcomp(X_cv, center = FALSE, scale = FALSE)
          cv_rot_mat <- cv_pr_out$rotation[,1:n_PCs]

          # PC-transform the train and test datasets, then merge into a single matrix for convenience
          # (The order is reversed here for efficient block matmul)
          X_outer_cv <- matrix(NA, TT, p)
          X_outer_cv[cv_inds,] <- X_cv 
          X_outer_cv[test_inds,] <- X_test 
          X_outer_cv_pc <- X_outer_cv %*% cv_rot_mat
          save_X(X_outer_cv_pc, X_dir, n_PCs = n_PCs, ofold = ofold_ind)

          # RWNN-transform the train and test datasets
          # (We can do this as a single operation since this transformation is fixed, not learned during model training)
          set.seed(seed)
          X_outer_cv_nl <- make_hidden_nodes(X_outer_cv_pc, n_h, a)
          save_X(X_outer_cv_nl, X_dir, n_PCs = n_PCs, n_h = n_h, seed = seed, ofold = ofold_ind)

          # Nested cross-validation (for application)
          if(n_PCs == 9 & n_h == 70) {

            for(ifold_ind in 1:nfold) {

              # Split the data
              train_inds <- unlist(ifolds_inner_inds[[ofold_ind]][-ifold_ind])
              val_inds <- ifolds_inner_inds[[ofold_ind]][[ifold_ind]]

              X_train <- X_cv[train_inds,]
              X_val <- X_cv[val_inds,]

              # Calculate the PC rotation matrix using only training data
              train_pr_out <- prcomp(X_train, center = FALSE, scale = FALSE)
              train_rot_mat <- train_pr_out$rotation[,1:n_PCs]

              # PC-transform the train and test datasets, then merge into a single matrix for convenience
              # (The order is reversed here for efficient block matmul)
              X_inner_cv <- matrix(NA, nrow(X_cv), p)
              X_inner_cv[train_inds,] <- X_train 
              X_inner_cv[val_inds,] <- X_val 
              X_inner_cv_pc <- X_inner_cv %*% train_rot_mat
              save_X(X_inner_cv_pc, X_dir, n_PCs = n_PCs, ofold = ofold_ind, ifold = ifold_ind)

              # RWNN-transform the train and test datasets
              # (We can do this as a single operation since this transformation is fixed, not learned during model training)
              set.seed(seed)
              X_inner_cv_nl <- make_hidden_nodes(X_inner_cv_pc, n_h, a)
              save_X(X_inner_cv_nl, X_dir, n_PCs = n_PCs, n_h = n_h, seed = seed, ofold = ofold_ind, ifold = ifold_ind)
            }
          }
        }
      }
    }
  }
}
