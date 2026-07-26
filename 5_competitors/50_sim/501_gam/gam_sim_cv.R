library(parallel)
library(mgcv)
library(dplyr)
library(matrixStats)

nfold <- 5
cv_gridsize <- 10
blocksize <- 20
TT <- 296
d <- 3
K <- 2

min_gamma <- 0.1
max_gamma <- 10000

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
mn_true <- abind::abind(datobj$mean_spec, datobj$pico_mu, along = 3)
rm(datobj)

real_dataobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper.RDS"))
X <- real_dataobj$X
X_df <- as.data.frame(X[,3:ncol(X)])
rm(real_dataobj)

# Load the clustering results
load(file.path("5_competitors", "50_sim", "500_k_means", "best_kmeans__scores-matched.Rdata"))
rea
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
  # TODO: test that cbind with Time works
  res <- cbind(cur_y, xt_mat, Weight, Cluster, Time) |> 
    as.data.frame()
  
  return(res)
}) %>% bind_rows()
# TODO: save for future use

# head(data_long)

data_long_by_clust_list <- data_long |> 
  group_by(Cluster) |> 
  group_split(.keep = FALSE)

# Last three columns are Weight, Cluster, and Time; we don't want these as covariates
rhs <- paste0("s(", colnames(data_long)[4:(ncol(data_long) - 3)], ", bs = \"cs\")") |> 
  paste(collapse = " + ")

y1_formula <- paste0("y1", " ~ ", rhs) |> 
  as.formula()

y2_formula <- paste0("y2", " ~ ", rhs) |> 
  as.formula()

y3_formula <- paste0("y3", " ~ ", rhs) |> 
  as.formula()

# TODO: setup cross-validation grid and folds
gamma_grid <- flowmix::logspace(min_gamma, max_gamma, cv_gridsize)
folds <- flowmix::make_cv_folds(nfold = nfold, blocksize = blocksize, TT = TT)

cv_jobs <- expand.grid(
  list(
    igamma = 1:cv_gridsize, 
    ifold = 1:nfold
  )
)

prob <- matrix(0, length(countslist), 2)
for(tt in 1:TT) {
  totals <- rowsum(countslist[[tt]], best_kmeans[[tt]]$cluster)
  totals <- totals / sum(totals)
  prob[tt, as.integer(rownames(totals))] <- as.numeric(totals)
}

# Cross-validate using the weighted negative log likelihood (assuming Gaussian mixture)
cvscores <- matrix(NA, cv_gridsize, nfold)
for(ijob in 1:nrow(cv_jobs)) {
  igamma <- cv_jobs[ijob, "igamma"]
  ifold <- cv_jobs[ijob, "ifold"]

  in_sample_inds <- unlist(folds[-ifold])
  out_sample_inds <- folds[[ifold]]

  # Set up vectorized data to easily calculate the likelihood
  out_sample_data <- filter(data_long, Time %in% out_sample_inds)
  n_per_tt_out <- sapply(out_sample_inds, function(tt) {
    length(countslist[[tt]])
  })
  probs_long <- prob[rep(out_sample_inds, times = n_per_tt_out),]
  counts_out <- unlist(countslist[out_sample_inds])
  y1 <- sapply(ybin_list[out_sample_inds], function(cur_y) cur_y[,1])
  y1 <- unlist(y1)
  y2 <- sapply(ybin_list[out_sample_inds], function(cur_y) cur_y[,2])
  y2 <- unlist(y2)
  y3 <- sapply(ybin_list[out_sample_inds], function(cur_y) cur_y[,3])
  y3 <- unlist(y3)

  likelihood_contributions <- matrix(0, sum(n_per_tt_out), 2)
  for(k in 1:length(data_long_by_clust_list)) {
    cur_data_long <- data_long_by_clust_list[[k]]
    cur_data_long_in <- filter(cur_data_long, Time %in% in_sample_inds)

    cur_y1_gam <- bam(
      y1_formula, 
      data = cur_data_long_in, 
      weights = Weight, 
      select = TRUE, 
      gamma = gamma_grid[igamma]
    )

    y1_preds <- predict(cur_y1_gam, newdata = out_sample_data, type = "response")

    cur_y2_gam <- bam(
      y2_formula, 
      data = cur_data_long_in, 
      weights = Weight, 
      select = TRUE, 
      gamma = gamma_grid[igamma]
    )

    y2_preds <- predict(cur_y2_gam, newdata = out_sample_data, type = "response")
    
    cur_y3_gam <- bam(
      y3_formula, 
      data = cur_data_long_in, 
      weights = Weight, 
      select = TRUE, 
      gamma = gamma_grid[igamma]
    )

    y3_preds <- predict(cur_y3_gam, newdata = out_sample_data, type = "response")

    likelihood_contributions[,k] <- log(probs_long[,k]) + 
      dnorm(y1, y1_preds, sqrt(cur_y1_gam$sig2), log = TRUE) + 
      dnorm(y2, y2_preds, sqrt(cur_y2_gam$sig2), log = TRUE) + 
      dnorm(y3, y3_preds, sqrt(cur_y3_gam$sig2), log = TRUE)
  }

  likelihood_contributions <- rowLogSumExps(likelihood_contributions)
  WNLL <- -sum(counts_out * likelihood_contributions)

  cvscores[igamma,ifold] <- WNLL
}

saveRDS(cvscores, file = "5_competitors/50_sim/501_gam/gam_cvscores.RDS")

avg_scores <- rowMeans(cvscores)
avg_scores
gamma_grid[which.min(avg_scores)]

# TODO: refit with gamma = 1000 and see what the model is like
# TODO: do one big cross-validation from 1000 to 0.1 and see where it levels off

plot(gamma_grid, log(avg_scores - min(avg_scores) + 1))
gamma_grid
avg_scores
gamma_grid[7]
# 148.7352
# levels out, then right after this it goes up

plot(gamma_grid, log(cvscores[,1] + abs(min(cvscores)) + 1), type = "l", ylim = c(0, 16))
lines(gamma_grid, log(cvscores[,2] + abs(min(cvscores)) + 1), ylim = c(0, 16))
lines(gamma_grid, log(cvscores[,3] + abs(min(cvscores)) + 1), ylim = c(0, 16))
lines(gamma_grid, log(cvscores[,4] + abs(min(cvscores)) + 1), ylim = c(0, 16))
lines(gamma_grid, log(cvscores[,5] + abs(min(cvscores)) + 1), ylim = c(0, 16))

gamma_grid[3] # [1] 529.8317

# TODO: refit with this value of gamma and see how it looks

# TODO: try more penalization
new_cvscores <- numeric(nfold)
for(ifold in 1:nfold) {
  new_gamma <- 10000

  in_sample_inds <- unlist(folds[-ifold])
  out_sample_inds <- folds[[ifold]]

  # Set up vectorized data to easily calculate the likelihood
  out_sample_data <- filter(data_long, Time %in% out_sample_inds)
  n_per_tt_out <- sapply(out_sample_inds, function(tt) {
    length(countslist[[tt]])
  })
  probs_long <- prob[rep(out_sample_inds, times = n_per_tt_out),]
  counts_out <- unlist(countslist[out_sample_inds])
  y1 <- sapply(ybin_list[out_sample_inds], function(cur_y) cur_y[,1])
  y1 <- unlist(y1)
  y2 <- sapply(ybin_list[out_sample_inds], function(cur_y) cur_y[,2])
  y2 <- unlist(y2)
  y3 <- sapply(ybin_list[out_sample_inds], function(cur_y) cur_y[,3])
  y3 <- unlist(y3)

  likelihood_contributions <- matrix(0, sum(n_per_tt_out), 2)
  for(k in 1:length(data_long_by_clust_list)) {
    cur_data_long <- data_long_by_clust_list[[k]]
    cur_data_long_in <- filter(cur_data_long, Time %in% in_sample_inds)

    cur_y1_gam <- bam(
      y1_formula, 
      data = cur_data_long_in, 
      weights = Weight, 
      select = TRUE, 
      gamma = new_gamma
    )

    y1_preds <- predict(cur_y1_gam, newdata = out_sample_data, type = "response")

    cur_y2_gam <- bam(
      y2_formula, 
      data = cur_data_long_in, 
      weights = Weight, 
      select = TRUE, 
      gamma = new_gamma
    )

    y2_preds <- predict(cur_y2_gam, newdata = out_sample_data, type = "response")
    
    cur_y3_gam <- bam(
      y3_formula, 
      data = cur_data_long_in, 
      weights = Weight, 
      select = TRUE, 
      gamma = new_gamma
    )

    y3_preds <- predict(cur_y3_gam, newdata = out_sample_data, type = "response")

    likelihood_contributions[,k] <- log(probs_long[,k]) + 
      dnorm(y1, y1_preds, sqrt(cur_y1_gam$sig2), log = TRUE) + 
      dnorm(y2, y2_preds, sqrt(cur_y2_gam$sig2), log = TRUE) + 
      dnorm(y3, y3_preds, sqrt(cur_y3_gam$sig2), log = TRUE)
  }

  likelihood_contributions <- rowLogSumExps(likelihood_contributions)
  WNLL <- -sum(counts_out * likelihood_contributions)

  new_cvscores[ifold] <- WNLL
}

new_cvscores_mat <- rbind(new_cvscores, cvscores)

plot(c(10000, gamma_grid), log(new_cvscores_mat[,1] + abs(min(new_cvscores_mat)) + 1), type = "l", ylim = c(0, 16))
lines(c(10000, gamma_grid), log(new_cvscores_mat[,2] + abs(min(new_cvscores_mat)) + 1), ylim = c(0, 16))
lines(c(10000, gamma_grid), log(new_cvscores_mat[,3] + abs(min(new_cvscores_mat)) + 1), ylim = c(0, 16))
lines(c(10000, gamma_grid), log(new_cvscores_mat[,4] + abs(min(new_cvscores_mat)) + 1), ylim = c(0, 16))
lines(c(10000, gamma_grid), log(new_cvscores_mat[,5] + abs(min(new_cvscores_mat)) + 1), ylim = c(0, 16))