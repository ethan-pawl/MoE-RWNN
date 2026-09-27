library(randomForestSRC)
library(abind)
library(dplyr)

# Load the data
datobj <- readRDS(
  file.path(
    "4_3dsim", 
    "3dsimdata", 
    "simdata_imean_2_iprob_1_iint_5_dataSeed_1.RDS"
  )
)

ybin_list <- datobj$ybin_list
countslist <- datobj$countslist

# Store ground truth for later
mn_true <- abind(datobj$mean_spec, datobj$pico_mu, along = 3)
prob1_true <- datobj$prob_spec
cov_true <- abind(datobj$clust1_cov, datobj$clust2_cov, along = 3)

rm(datobj)

# Load covariates
real_dataobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper.RDS"))
X <- real_dataobj$X
X_df <- as.data.frame(X[,3:ncol(X)])

rm(real_dataobj)

ybin_list <- lapply(seq_along(ybin_list), function(tt) {
  cur_y <- as.data.frame(ybin_list[[tt]])
  colnames(cur_y) <- c("y1", "y2", "y3")
  
  return(cur_y)
})

# Cluster each cytogram using sidCluster
set.seed(0)
clust_res <- lapply(seq_along(ybin_list), function(tt) {
  print(tt)
  sidClustering(
    ybin_list[[tt]], 
    k = 2, 
    case.wt = countslist[[tt]], 
    reduce = FALSE
  )$clustering
})

# Output is a length TT list of length nt vectors of cluster assignments
TT <- length(ybin_list)
K <- 2

# Match clusters across time
# Use flowMatch helper functions to match clusters
for(tt in 2:TT) {
  print(tt)
  prev_samp <- flowMatch::ClusteredSample(labels = clust_res[[tt-1]], sample = ybin_list[[tt-1]])
  cur_samp <- flowMatch::ClusteredSample(labels = clust_res[[tt]], sample = ybin_list[[tt]])
  
  D <- flowMatch::dist.matrix(prev_samp, cur_samp, dist.type = "KL")
  rownames(D) <- colnames(D) <- as.character(1:K)

  cur_samp_labels <- integer(K)
  for(k in 1:K) {
    # Greedily match clusters based on minimal KL divergence
    min_ind <- which.min(D) |> arrayInd(dim(D))
    cur_samp_labels[as.integer(colnames(D)[min_ind[,2]])] <- as.integer(rownames(D)[min_ind[,1]])
    D <- D[-min_ind[,1],-min_ind[,2], drop = FALSE]
  }

  reordering <- order(cur_samp_labels)
  clust_res[[tt]] <- clust_res[[tt]] |> recode_values(from = 1:K, to = cur_samp_labels)
}

# Data prep 

data_long <- lapply(seq_along(ybin_list), function(tt) {
  cur_y <- ybin_list[[tt]]
  nt <- nrow(cur_y)

  which_covs <- 3:ncol(X)
  xt_mat <- matrix(rep(X[tt,which_covs], each = nt), nt)
  colnames(xt_mat) <- colnames(X)[which_covs]

  Weight <- countslist[[tt]]
  Cluster <- clust_res[[tt]]
  Time <- tt
  res <- cbind(cur_y, xt_mat, Weight, Cluster, Time) |> 
    as.data.frame()
  
  return(res)
}) %>% bind_rows()

data_long_by_clust_list <- data_long |> 
  group_by(Cluster) |> 
  group_split(.keep = FALSE)

# Regress each cluster on covariates using multivariate random forests
# Multithreading is too memory-intensive
RhpcBLASctl::omp_set_num_threads(1L)
RhpcBLASctl::blas_set_num_threads(1L)

set.seed(0)
rf_res <- lapply(data_long_by_clust_list, function(cur_data_long) {
  print("cluster")
  clust_rf_res <- rfsrc.fast(
    cbind(y1, y2, y3) ~ . - Time - Weight,
    data = cur_data_long,
    case.wt = cur_data_long$Weight, 
    forest = TRUE,
    do.trace = 5
  )

  return(clust_rf_res)
})

# Analyze the results
TT <- length(ybin_list)
d <- ncol(ybin_list[[1]])
K <- length(rf_res)

# Make mean predictions
mn_fit <- array(NA, c(TT, d, K))
Time <- rep(1, 296)
Weight <- rep(1, 296)
for(k in 1:K) {
  cur_preds <- predict(
    rf_res[[k]],
    newdata = cbind(X_df, Time, Weight)
  )

  mn_fit[,,k] <- cbind(
    cur_preds$regrOutput$y1$predicted,
    cur_preds$regrOutput$y2$predicted,
    cur_preds$regrOutput$y3$predicted
  )
}

# Empirical biomass proportions
prob <- matrix(0, length(countslist), 2)
for(tt in 1:TT) {
  totals <- rowsum(countslist[[tt]], clust_res[[tt]])
  totals <- totals / sum(totals)
  prob[tt, as.integer(rownames(totals))] <- as.numeric(totals)
}

# Calculate error metrics

# Calculate RMSE (root mean l2 error)
resids <- mn_true - mn_fit
resid_dotprods <- apply(resids, 3, function(x) {
  crossprod(as.vector(t(x)))
})

rmse <- sqrt(resid_dotprods / TT)

prob_rmse <- sqrt(mean((prob1_true - prob[,1])^2))

# Covariance estimates
resid_cov <- array(NA, c(d, d, K))
for(k in 1:K) {
  clust_resid <- lapply(1:TT, function(tt) {
    t(t(ybin_list[[tt]][clust_res[[tt]] == k,]) - mn_fit[tt,,k]) # Transpose to recycle across columns, then tranpose back to maintain order
  }) %>% do.call(rbind, .)

  resid_cov[,,k] <- cov(clust_resid)  
}

cov_err <- sapply(1:K, function(k) {
  sqrt(sum((cov_true[,,k] - resid_cov[,,k])^2))
})

results <- data.frame(
  Model = "Random Forest",
  Metric = rep(c("RMSE, Mean", "RMSE, Probability", "Frobenius Error, Covariance"), times = c(3, 1, 3)), 
  Cluster = c("1", "2", "Total", "1", "1", "2", "Total"),
  Value = c(rmse, sum(rmse), prob_rmse, cov_err, sum(cov_err))
)

write.csv(results, file.path("6_competitors", "metrics", "rf.csv"), row.names = FALSE)
