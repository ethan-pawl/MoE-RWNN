library(ggplot2)

cv_files <- list.files("results", full.names = TRUE)

nl_res_list <- list()
for(i in 1:5) {
    nl_res_list[[i]] <- readRDS(cv_files[i])
}

l_res_list <- list()
for(i in 1:5) {
    l_res_list[[i]] <- readRDS(cv_files[i + 5])
}

library(flowmix)

# Load all data
datobj = readRDS(file.path("~", 
                           "00_Cyto", 
                           "data", 
                           "paper-data-v2", 
                           "MGL1704-hourly-paper.RDS"))

ylist <- datobj$ylist
countslist <- datobj$countslist

# Get indices of out-of-sample data
load(file.path("~", 
               "00_Cyto", 
               "scripts", 
               "2024_02_28_cv", 
               "out_sample_folds.RData"))

# Log-likelihoods
X_dir <- file.path("~", 
                   "00_Cyto", 
                   "data", 
                   "X_data")

load(file.path(X_dir,
               "X_pc.Rdata"))
X_pc <- X

load(file.path(X_dir,
               "X_nl_70.Rdata"))
X_nl <- X

l_oos_nll <- nl_oos_nll <- numeric(5)
for(i in 1:5) {
    l_pred <- predict(l_res_list[[i]]$bestres, newx = X_pc[out_sample_folds[[i]],])
    l_oos_nll[i] <- objective_newdat(l_pred, ylist[out_sample_folds[[i]]], 
                                     countslist[out_sample_folds[[i]]])

    nl_pred <- predict(nl_res_list[[i]]$bestres, newx = X_nl[out_sample_folds[[i]],])
    nl_oos_nll[i] <- objective_newdat(nl_pred, ylist[out_sample_folds[[i]]], 
                                     countslist[out_sample_folds[[i]]])
}

mean(l_oos_nll)
mean(nl_oos_nll)

nll_df <- data.frame(NLL = c(l_oos_nll, nl_oos_nll), 
                     Model = rep(c("Linear", "Nonlinear"), each = 5), 
                     Heldout = rep(1:5, 2))

ggplot(nll_df, aes(Model, NLL, fill = Model)) + 
    geom_boxplot() + 
    geom_point()

# figure out which held out point the nonlinear does so much better
# out_3 <- out_sample_folds[[3]]
# 
# plot((1:296)[-out_3], datobj$lat[-out_3], ylim = range(datobj$lat))
# points((1:296)[out_3], datobj$lat[out_3], col = "red")
