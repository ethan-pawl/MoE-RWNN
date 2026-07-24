library(parallel)
library(mgcv)
library(dplyr)

nfold <- 5
cv_gridsize <- 30
blocksize <- 20
TT <- 296
d <- 3
K <- 2

min_gamma <- 0.1
max_gamma <- 1000

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

# TODO: insert code which loads the data and creates the data frame

# TODO: might have to slurm this
cvscores <- matrix(NA, cv_gridsize, nfold)
for(ijob in 1:nrow(cv_jobs)) {
  igamma <- cv_jobs[ijob, "igamma"]
  ifold <- cv_jobs[ijob, "ifold"]

  in_sample_inds <- unlist(folds[-ifold])
  out_sample_inds <- folds[[ifold]]

  y1_res <- lapply(data_long_by_clust_list, function(cur_data_long) {
    cur_y1_gam <- bam(
      y1_formula, 
      data = filter(cur_data_long, Time %in% in_sample_inds), 
      weights = Weight, 
      select = TRUE,
      gamma = gamma_grid[igamma]
    )

    return(cur_y1_gam)
  })

  y2_res <- lapply(data_long_by_clust_list, function(cur_data_long) {
    cur_y2_gam <- bam(
      y2_formula, 
      data = filter(cur_data_long, Time %in% in_sample_inds), 
      weights = Weight, 
      select = TRUE,
      gamma = gamma_grid[igamma]
    )

    return(cur_y2_gam)
  })

  y3_res <- lapply(data_long_by_clust_list, function(cur_data_long) {
    cur_y3_gam <- bam(
      y3_formula, 
      data = filter(cur_data_long, Time %in% in_sample_inds), 
      weights = Weight, 
      select = TRUE,
      gamma = gamma_grid[igamma]
    )

    return(cur_y3_gam)
  })

  mn_fit <- array(NA, c(length(out_sample_inds), d, K))
  for(j in 1:3) {
    res <- get(paste0("y", j, "_res"))
    for(k in 1:K) {
      mn_fit[,j,k] <- predict(
        res[[k]],
        newdata = X_df[out_sample_inds,], 
        type = "response"
      )
    }
  }

  resids <- mn_true[out_sample_inds,,] - mn_fit
  resid_dotprods <- apply(resids, 3, function(x) {
    crossprod(as.vector(t(x)))
  })

  cvscores[igamma,ifold] <- sum(resid_dotprods)
  print(cvscores)
}

saveRDS(cvscores, file = "5_competitors/50_sim/501_gam/gam_cvscores.RDS")

avg_scores <- rowMeans(cvscores)
avg_scores
gamma_grid[which.min(avg_scores)]

# TODO: refit with gamma = 1000 and see what the model is like
# TODO: do one big cross-validation from 1000 to 0.1 and see where it levels off

plot(gamma_grid, log(avg_scores))
gamma_grid
avg_scores
gamma_grid[7]
# 148.7352
# levels out, then right after this it goes up

plot(gamma_grid, log(cvscores[,1]), type = "l")
lines(gamma_grid, log(cvscores[,2]))
lines(gamma_grid, log(cvscores[,3]))
lines(gamma_grid, log(cvscores[,4]))
lines(gamma_grid, log(cvscores[,5]))

gamma_grid[3] # [1] 529.8317

# TODO: refit with this value of gamma and see how it looks
