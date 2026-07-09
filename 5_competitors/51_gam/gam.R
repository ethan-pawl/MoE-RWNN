library(dplyr)

# Load the data
datobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper.RDS"))
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()

# Load the clustering results
load(file.path("5_competitors", "50_k_means", "best_kmeans__scores-matched.Rdata"))

# FIGUREOUT: what's the closest nlme option?
# - Nonlinear
# - Lasso-penalized
# - PCA?
# - nonparametric: stats::selfStart 
# - or try mgcv::gam 
#   - advantage of mgcv is that it has a `select` argument and `gamma` argument to do variable selection 
#   - another advantage is that we can use nonparametric smoothers, like splines & GPs
#   - advantage of glmnet is that we can do justice to reviewer's comment about elastic net (check this is the correct reviewer)
#     - but I'd have to do slightly more work to specify nonlinear terms (or could just do formula with selfStart -> model.matrix)
# - probably do all of these and report at least the best one
# selfStart: try all of them, but prioritize asymp, logis, gompertz/weibull for more general shape

# TODO: first do the easiest possible thing, which is mgcv::gam with s(x1) + s(x2) + ..., 
# cross-validating in the same way as flowmix. Don't need to worry about lambda alpha here.
# FIGUREOUT: do I need to do nested CV?

library(mgcv)

# Replicate x_t for every t to get a tabular dataset

data_long <- lapply(seq_along(time), function(tt) {
  cur_y <- ylist[[tt]]
  nt <- nrow(cur_y)

  which_covs <- 3:ncol(X)
  xt_mat <- matrix(rep(X[tt,which_covs], each = nt), nt)
  colnames(xt_mat) <- colnames(X)[which_covs]
  
  res <- cbind(cur_y, xt_mat) |> 
    as.data.frame()
  
  return(res)
}) %>% bind_rows()

head(data_long)

counts_long <- unlist(unname(countslist))

rhs <- paste0("s(", colnames(data_long)[4:ncol(data_long)], ")") |> 
  paste(collapse = " + ")

diam_formula <- paste0("diam_mid", " ~ ", rhs) |> 
  as.formula()

# FIXME: this crashes
# diam_gam_1 <- gam(
#   diam_formula, 
#   data = data_long, 
#   weights = counts_long, 
#   method = "BFGS"
# )


# TODO: use mvn instead for multivariate 
# responses

rhs_2 <- paste0("s(", colnames(data_long)[4:ncol(data_long)], ", bs = \"cc\")") |> 
  paste(collapse = " + ")

diam_formula_2 <- paste0("diam_mid", " ~ ", rhs_2) |> 
  as.formula()

rhs_3 <- paste0("s(", colnames(data_long)[4:ncol(data_long)], ", bs = \"cs\")") |> 
  paste(collapse = " + ")

diam_formula_3 <- paste0("diam_mid", " ~ ", rhs_3) |> 
  as.formula()

# print("TPRS Time:")
# system.time({
#   diam_gam_1 <- bam(
#     diam_formula, 
#     data = data_long, 
#     weights = counts_long
#   )
# })
# # elapsed: 116.776 

# print("Cyclic Cubic Spline Time:")
# system.time({
#   diam_gam_2 <- bam(
#     diam_formula_2, 
#     data = data_long, 
#     weights = counts_long
#   )
# })
# # elapsed: 90.648

# print("Cubic Spline Time:")
# system.time({
#   diam_gam_3 <- bam(
#     diam_formula_3, 
#     data = data_long, 
#     weights = counts_long
#   )
# })
# elapsed: NAs produced after 64.25 seconds

# Cyclic cubic splines seem to be the way 
# to go

# Move to multivariate, see if feasible
chl_formula_2 <- paste0("chl_small", " ~ ", rhs_2) |> 
  as.formula()

pe_formula_2 <- paste0("pe", " ~ ", rhs_2) |> 
  as.formula()

# # TODO: split into each cluster
# print("Multivariate Cyclic Cubic Spline Time:")
# system.time({
#   mv_gam <- bam(
#     list(
#       diam_formula_2, 
#       chl_formula_2, 
#       pe_formula_2
#     ),
#     family = mvn(d = 3),
#     data = data_long
#   )
# })

# Problem: bam can't use mvn

data_long <- lapply(seq_along(time), function(tt) {
  cur_y <- ylist[[tt]]
  nt <- nrow(cur_y)

  which_covs <- 3:ncol(X)
  xt_mat <- matrix(rep(X[tt,which_covs], each = nt), nt)
  colnames(xt_mat) <- colnames(X)[which_covs]

  Weight <- countslist[[tt]]
  Cluster <- best_kmeans[[tt]]$cluster
  res <- cbind(cur_y, xt_mat, Weight, Cluster) |> 
    as.data.frame()
  
  return(res)
}) %>% bind_rows()

data_long$Cluster <- factor(data_long$Cluster)
head(data_long)

data_long_by_clust_list <- data_long |> 
  group_by(Cluster) |> 
  group_split()

RhpcBLASctl::blas_set_num_threads(1)
system.time({
  per_cluster_res <- parallel::mclapply(data_long_by_clust_list, function(cur_data_long) {
    cur_diam_gam <- bam(
      diam_formula_3, 
      data = cur_data_long, 
      weights = Weight, 
      select = TRUE, 
      gamma = 1
    )

    return(cur_diam_gam)
  }, mc.cores = 1, mc.preschedule = FALSE)
})
# select = TRUE seems to fix this
# TODO: cross-validate over gamma

cur_data_long <- data_long_by_clust_list[[1]]
head(cur_data_long)

# TODO: for fixed gamma, make the relevant 
# plots

# TODO: also try per-cluster multivariate 
# plots