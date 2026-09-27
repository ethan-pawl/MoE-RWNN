library(flowmix)
library(abind)

datobj <- readRDS(
  file.path(
    "4_3dsim", 
    "3dsimdata", 
    "simdata_imean_2_iprob_1_iint_5_dataSeed_1.RDS"
  )
)

ybin_list <- datobj$ybin_list
countslist <- datobj$countslist

# Estimation settings
maxdev <- diff(range(datobj$mean_spec[,1])) / 2 
numclust <- 2
nrep <- 10

X_dir <- file.path("data", "X_variations")
X <- readRDS(file.path(X_dir, "X_pc_9_nh_NA_seed_NA_ofold_NA_ifold_NA.RDS"))

load(file.path("6_competitors", "4_linear", "settings.Rdata"))

res <- flowmix_once(
  ylist = ybin_list, 
  X = X, 
  countslist = countslist, 
  numclust = numclust, 
  prob_lambda = prob_lambda, 
  mean_lambda = mean_lambda, 
  verbose = TRUE, 
  maxdev = maxdev,
  seed = seed
)

# Store ground truth
mn_true <- abind(datobj$mean_spec, datobj$pico_mu, along = 3)
prob1_true <- datobj$prob_spec
cov_true <- abind(datobj$clust1_cov, datobj$clust2_cov, along = 3)

# Calculate error metrics

# Calculate RMSE (root mean l2 error)
mn_fit <- res$mn
resids <- mn_true - mn_fit 
resid_dotprods <- apply(resids, 3, function(x) {
  crossprod(as.vector(t(x)))
})

TT <- 296
rmse <- sqrt(resid_dotprods / TT)

prob <- res$prob
prob_rmse <- sqrt(mean((prob1_true - prob[,1])^2))

K <- 2
cov_fit <- res$sigma |> aperm(c(2, 3, 1))
cov_err <- sapply(1:K, function(k) { 
  sqrt(sum((cov_true[,,k] - cov_fit[,,k])^2)) 
})

results <- data.frame(
  Model = "Linear Model",
  Metric = rep(c("RMSE, Mean", "RMSE, Probability", "Frobenius Error, Covariance"), times = c(3, 1, 3)), 
  Cluster = c("1", "2", "Total", "1", "1", "2", "Total"),
  Value = c(rmse, sum(rmse), prob_rmse, cov_err, sum(cov_err))
)

write.csv(results, file.path("6_competitors", "metrics", "linear_flowmix.csv"), row.names = FALSE)
