library(flowmix)
library(magrittr)
library(ggplot2)
library(reshape2)

# Collect results
l_files <- list.files(file.path("2_application", "20_nested_cv", "results"), 
                      recursive = TRUE, 
                      pattern = "nh_NA_seed_NA_ofold_[[:digit:]]\\.RDS$", 
                      full.names = TRUE)

l_res_list <- lapply(l_files, function(cur_fname) readRDS(cur_fname))

nl_files <- list.files(file.path("2_application", "20_nested_cv", "results"), 
                      recursive = TRUE, 
                      pattern = "nh_70_seed_1_ofold_[[:digit:]]\\.RDS$", 
                      full.names = TRUE)

nl_res_list <- lapply(nl_files, function(cur_fname) readRDS(cur_fname))

# Load all data
datobj <- readRDS(file.path("data", "MGL1704-hourly-paper.RDS"))
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()

# Get indices of out-of-sample data
X_dir <- file.path(
  "data", 
  "X_variations"
)

load(file.path("data", "ofolds__ifolds__ifolds_inner_inds.Rdata"))

load_X <- function(readdir, n_PCs = "NA", n_h = "NA", seed = "NA", ofold = "NA", ifold = "NA") {

  file_name <- paste("X_pc", n_PCs, "nh", n_h, "seed", seed, "ofold", ofold, "ifold", ifold, sep = "_")
  file_name <- paste0(file_name, ".RDS")
  file_name <- file.path(readdir, file_name)

  X <- readRDS(file_name)
}

############################

# Calculate out-of-sample negative log-likelihoods (OOS NLLs)
l_oos_nll <- nl_oos_nll <- numeric(5)
for(ofold_ind in 1:5) { # For each outer fold (held-out test dataset)
    X_pc <- load_X(X_dir, n_PCs = 9, ofold = ofold_ind)
    X_nl <- load_X(X_dir, n_PCs = 9, n_h = 70, seed = 1, ofold = ofold_ind)

    outsample_inds <- ofolds[[ofold_ind]]

    # Estimate linear model parameters for out-of-sample data
    l_pred <- predict(l_res_list[[ofold_ind]]$bestres, newx = X_pc[outsample_inds,])    
    # Calculate linear model OOS NLL based on estimated parameters
    l_oos_nll[ofold_ind] <- objective_newdat(l_pred, ylist[outsample_inds], 
                                     countslist[outsample_inds])

    # Estimate nonlinear model parameters for out-of-sample data
    nl_pred <- predict(nl_res_list[[ofold_ind]]$bestres, newx = X_nl[outsample_inds,])
    # Calculate nonlinear model OOS NLL based on estimated parameters
    nl_oos_nll[i] <- objective_newdat(nl_pred, ylist[outsample_inds], 
                                     countslist[outsample_inds])
}

# These are the estimates of out-of-sample predictive performance (NLPLs) 
# mentioned in the article
mean(l_oos_nll)
mean(nl_oos_nll)

res_df <- data.frame(Linear = l_oos_nll, Nonlinear = nl_oos_nll, `Outer Fold` = 1:5, check.names = FALSE)
res_df <- melt(res_df, id.vars = "Outer Fold", variable.name = "Model", value.name = "NLL")

ggplot(res_df) + 
    geom_boxplot(aes(Model, NLL))
