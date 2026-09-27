
library(mgcv)
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

# Store ground truth params for later
mn_true <- abind(datobj$mean_spec, datobj$pico_mu, along = 3)
prob1_true <- datobj$prob_spec
cov_true <- abind(datobj$clust1_cov, datobj$clust2_cov, along = 3)

# Free up space
rm(datobj)

# Load covariates
real_dataobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper.RDS"))
X <- real_dataobj$X
X_df <- as.data.frame(X[,3:ncol(X)])

# Free up space
rm(real_dataobj)

# Load the clustering results
load(file.path("6_competitors", "0_k_means", "best_kmeans__scores-matched.Rdata"))

# Data prep
data_long <- lapply(seq_along(ybin_list), function(tt) {
  cur_y <- ybin_list[[tt]]
  colnames(cur_y) <- c("y1", "y2", "y3")
  nt <- nrow(cur_y)

  which_covs <- 3:ncol(X)
  xt_mat <- matrix(rep(X[tt,which_covs], each = nt), nt)
  colnames(xt_mat) <- colnames(X)[which_covs]

  Weight <- countslist[[tt]]
  Cluster <- best_kmeans[[tt]]$cluster
  Time <- tt
  res <- cbind(cur_y, xt_mat, Weight, Cluster, Time) |> 
    as.data.frame()
  
  return(res)
}) %>% bind_rows()

data_long_by_clust_list <- data_long |> 
  group_by(Cluster) |> 
  group_split(.keep = FALSE)

# Model setup 

# Last three columns are Weight, Cluster, and Time; we don't want these as covariates
rhs <- paste0("s(", colnames(data_long)[4:(ncol(data_long) - 3)], ", bs = \"cs\")") |> 
  paste(collapse = " + ")

y1_formula <- paste0("y1", " ~ ", rhs) |> 
  as.formula()

y2_formula <- paste0("y2", " ~ ", rhs) |> 
  as.formula()

y3_formula <- paste0("y3", " ~ ", rhs) |> 
  as.formula()

# Model fitting
y1_res <- lapply(data_long_by_clust_list, function(cur_data_long) {
  cur_y1_gam <- bam(
    y1_formula, 
    data = cur_data_long, 
    weights = Weight, 
    select = TRUE,
    gamma = 529.8317
  )

  return(cur_y1_gam)
})

y2_res <- lapply(data_long_by_clust_list, function(cur_data_long) {
  cur_y2_gam <- bam(
    y2_formula, 
    data = cur_data_long, 
    weights = Weight, 
    select = TRUE,
    gamma = 529.8317
  )

  return(cur_y2_gam)
})

y3_res <- lapply(data_long_by_clust_list, function(cur_data_long) {
  cur_y3_gam <- bam(
    y3_formula, 
    data = cur_data_long, 
    weights = Weight, 
    select = TRUE,
    gamma = 529.8317
  )

  return(cur_y3_gam)
})

per_cluster_res <- list(
  y1_res = y1_res, 
  y2_res = y2_res, 
  y3_res = y3_res
)

# Model predictions
TT <- length(ybin_list)
d <- length(per_cluster_res)
K <- length(per_cluster_res[[1]])

mn_fit <- array(NA, c(TT, d, K))

for(j in 1:d) {
  for(k in 1:K) {
    mn_fit[,j,k] <- predict(
      per_cluster_res[[j]][[k]],
      newdata = X_df, 
      type = "response"
    )
  }
}

# Calculate error metrics 

# Calculate RMSE (root mean l2 error)
resids <- mn_true - mn_fit[,,2:1] # LABELS SWITCHED
resid_dotprods <- apply(resids, 3, function(x) {
  crossprod(as.vector(t(x)))
})

rmse <- sqrt(resid_dotprods / TT)

# Empirical biomass proportions
prob <- matrix(0, length(countslist), 2)
for(tt in 1:TT) {
  totals <- rowsum(countslist[[tt]], best_kmeans[[tt]]$cluster)
  totals <- totals / sum(totals)
  prob[tt, as.integer(rownames(totals))] <- as.numeric(totals)
}

prob_rmse <- sqrt(mean((prob1_true - prob[,2])^2))

cov_fit <- sapply(1:K, function(k) {
  sapply(1:d, function(j) {
    per_cluster_res[[j]][[k]]$sig2
  }) |> diag()
}, simplify = "array")

cov_err <- sapply(1:K, function(k) {
  sqrt(sum((cov_true[,,k] - cov_fit[,,if(k == 1) 2 else 1])^2))
})

results <- data.frame(
  Model = "GAM",
  Metric = rep(c("RMSE, Mean", "RMSE, Probability", "Frobenius Error, Covariance"), times = c(3, 1, 3)), 
  Cluster = c("1", "2", "Total", "1", "1", "2", "Total"),
  Value = c(rmse, sum(rmse), prob_rmse, cov_err, sum(cov_err))
)

write.csv(results, file.path("6_competitors", "metrics", "gam.csv"), row.names = FALSE)
