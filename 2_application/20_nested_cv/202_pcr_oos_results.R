library(flowmix)
library(magrittr)

# Collect results
l_files <- list.files(file.path("2_application", "20_nested_cv", "results"), 
                      recursive = TRUE, 
                      pattern = "pcr_oos_linear_[[:digit:]].RDS", 
                      full.names = TRUE)

l_res_list <- lapply(l_files, function(cur_fname) readRDS(cur_fname))

nl_files <- list.files(file.path("2_application", "20_nested_cv", "results"), 
                      recursive = TRUE, 
                      pattern = "pcr_oos_nl_[[:digit:]].RDS", 
                      full.names = TRUE)

nl_res_list <- lapply(nl_files, function(cur_fname) readRDS(cur_fname))

# Load all data
datobj <- readRDS(file.path("data", "MGL1704-hourly-paper.RDS"))
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()

# Get indices of out-of-sample data
load(file.path("2_application", 
               "20_nested_cv", 
               "folds", 
               "out_sample_folds.RData"))

load(file.path("data", "X_pc.Rdata"))
X_pc <- X

load(file.path("data", "X_nl.Rdata"))
X_nl <- X

############################

# Calculate out-of-sample negative log-likelihoods (OOS NLLs)
l_oos_nll <- nl_oos_nll <- numeric(5)
for(i in 1:5) { # For each outer fold (held-out test dataset)

    # Estimate linear model parameters for out-of-sample data
    l_pred <- predict(l_res_list[[i]]$bestres, newx = X_pc[out_sample_folds[[i]],])    
    # Calculate linear model OOS NLL based on estimated parameters
    l_oos_nll[i] <- objective_newdat(l_pred, ylist[out_sample_folds[[i]]], 
                                     countslist[out_sample_folds[[i]]])

    # Estimate nonlinear model parameters for out-of-sample data
    nl_pred <- predict(nl_res_list[[i]]$bestres, newx = X_nl[out_sample_folds[[i]],])
    # Calculate nonlinear model OOS NLL based on estimated parameters
    nl_oos_nll[i] <- objective_newdat(nl_pred, ylist[out_sample_folds[[i]]], 
                                     countslist[out_sample_folds[[i]]])
}

# These are the estimates of out-of-sample predictive performance (NLPLs) 
# mentioned in the article
mean(l_oos_nll)
mean(nl_oos_nll)
