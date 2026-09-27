library(mixdistreg)
library(dplyr)
library(abind)

# Load the data
datobj <- readRDS(
  file.path(
    "4_3dsim", 
    "3dsimdata", 
    "simdata_imean_2_iprob_1_iint_5_dataSeed_1.RDS"
  )
)

# Load data
ybin_list <- datobj$ybin_list
countslist <- datobj$countslist

# Store ground truth params for later
mn_true <- abind(datobj$mean_spec, datobj$pico_mu, along = 3)
prob1_true <- datobj$prob_spec
cov_true <- abind(datobj$clust1_cov, datobj$clust2_cov, along = 3)

rm(datobj)

# Load covariates
real_dataobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper.RDS"))
X <- real_dataobj$X
X_df <- as.data.frame(X[,3:ncol(X)])
rm(real_dataobj)

# Data prep
ybin_list <- lapply(seq_along(ybin_list), function(tt) {
  cur_y <- as.data.frame(ybin_list[[tt]])
  colnames(cur_y) <- c("y1", "y2", "y3")
  
  return(cur_y)
})

data_long <- lapply(seq_along(ybin_list), function(tt) {
  cur_y <- ybin_list[[tt]]
  nt <- nrow(cur_y)

  which_covs <- 3:ncol(X)
  xt_mat <- matrix(rep(X[tt,which_covs], each = nt), nt)
  colnames(xt_mat) <- colnames(X)[which_covs]

  Weight <- countslist[[tt]]
  Time <- tt
  res <- cbind(cur_y, xt_mat, Weight, Time) |> 
    as.data.frame()
  
  return(res)
}) %>% bind_rows()

# PCA transformation

y_mat <- data_long[,1:3]
y_pca <- prcomp(y_mat, center = TRUE, scale = FALSE)
y_pcs <- y_pca$x

data_long <- cbind(y_pcs, data_long)
colnames(data_long)[1:3] <- paste0("PC", 1:3)

# Model setup

deep_model <- function(x) {
  x %>%
    layer_dense(
      units = 64,
      activation = "relu",
      kernel_regularizer = regularizer_l2(1e-3)
    ) %>%
    layer_dense(
      units = 64,
      activation = "relu",
      kernel_regularizer = regularizer_l2(1e-3)
    ) %>%
    layer_dense(
      units = 1,
      activation = "linear"
    )
}

deep_model_mix <- function(x) {
  x %>%
    layer_dense(
      units = 64,
      activation = "relu",
      kernel_regularizer = regularizer_l2(1e-5)
    ) %>%
    layer_dense(
      units = 64,
      activation = "relu",
      kernel_regularizer = regularizer_l2(1e-5)
    ) %>%
    layer_dense(
      units = 2,
      activation = "linear"
    )
}

PC1_formula_str <- paste0("~ 1 + deep_model(", paste(colnames(data_long[,7:(ncol(data_long) - 2)]), collapse = ","), ")")
PC1_formula <- as.formula(PC1_formula_str)

PC1_formula_mix_str <- paste0("~ 1 + deep_model_mix(", paste(colnames(data_long[,7:(ncol(data_long) - 2)]), collapse = ","), ")")
PC1_formula_mix <- as.formula(PC1_formula_mix_str)

# Model fitting
PC1_mod <- mixdistreg(
  data_long$PC1,
  families = "normal",
  nr_comps = 2L,
  list_of_formulas = list(
    loc_1 = PC1_formula, scale_1 = ~ 1, 
    loc_2 = PC1_formula, scale_2 = ~ 1
  ),
  trafos_each_param = list(
    list(
      function(x) x,
      function(x) tf$exp(x)
    ),
    list(
      function(x) x,
      function(x) tf$exp(x)
    )
  ),
  formula_mixture = PC1_formula_mix,
  list_of_deep_models = list(deep_model = deep_model, deep_model_mix = deep_model_mix),
  data = data_long, 
  optimizer = optimizer_adam()
)

RhpcBLASctl::omp_set_num_threads(1L)
RhpcBLASctl::blas_set_num_threads(1L)

history <- PC1_mod %>% fit(
  epochs = 500,
  early_stopping = TRUE, 
  sample_weight = data_long$Weight, 
  batch_size = 512, 
  callbacks = list(
    callback_early_stopping(
      monitor = "val_loss",
      patience = 10,
      restore_best_weights = TRUE
    )
  )
)

dist_obj <- get_distribution(PC1_mod)

# Gating probabilities
pi_hat <- as.matrix(
  tf$squeeze(
    dist_obj$submodules[[2]]$probs,
    axis = 1L
  )
)

# Component density terms
dens <- as.matrix(
  tf$squeeze(
    dist_obj$submodules[[1]]$prob(
      array(data_long$PC1, dim = c(nrow(data_long),1, 1))
    ), 
    axis = 1L
  )
)

# Multiply and normalize (Bayes' Rule)
post_probs <- pi_hat * dens
post_probs <- post_probs / rowSums(post_probs)
# rowSums(post_probs) is the posterior predictive distribution with the cluster assignments marginalized out (the mixture)

PC2_mod1 <- deepregression(
  y = data_long$PC2,
  family = "normal", 
  list_of_formulas = list(
    loc = PC1_formula,
    scale = ~ 1
  ), 
  list_of_deep_models = list(deep_model = deep_model),
  data = data_long, 
  optimizer = optimizer_adam()
)

history_21 <- PC2_mod1 %>% fit(
  epochs = 500,
  early_stopping = TRUE, 
  sample_weight = data_long$Weight * post_probs[,1], 
  batch_size = 512, 
  callbacks = list(
    callback_early_stopping(
      monitor = "val_loss",
      patience = 10,
      restore_best_weights = TRUE
    )
  )
)

PC2_mod2 <- deepregression(
  y = data_long$PC2,
  family = "normal", 
  list_of_formulas = list(
    loc = PC1_formula,
    scale = ~ 1
  ), 
  list_of_deep_models = list(deep_model = deep_model),
  data = data_long, 
  optimizer = optimizer_adam()
)

history_22 <- PC2_mod2 %>% fit(
  epochs = 500,
  early_stopping = TRUE, 
  sample_weight = data_long$Weight * post_probs[,2], 
  batch_size = 512, 
  callbacks = list(
    callback_early_stopping(
      monitor = "val_loss",
      patience = 10,
      restore_best_weights = TRUE
    )
  )
)

PC3_mod1 <- deepregression(
  y = data_long$PC3,
  family = "normal", 
  list_of_formulas = list(
    loc = PC1_formula,
    scale = ~ 1
  ), 
  list_of_deep_models = list(deep_model = deep_model),
  data = data_long, 
  optimizer = optimizer_adam()
)

history_31 <- PC3_mod1 %>% fit(
  epochs = 500,
  early_stopping = TRUE, 
  sample_weight = data_long$Weight * post_probs[,1], 
  batch_size = 512, 
  callbacks = list(
    callback_early_stopping(
      monitor = "val_loss",
      patience = 10,
      restore_best_weights = TRUE
    )
  )
)

PC3_mod2 <- deepregression(
  y = data_long$PC3,
  family = "normal", 
  list_of_formulas = list(
    loc = PC1_formula,
    scale = ~ 1
  ), 
  list_of_deep_models = list(deep_model = deep_model),
  data = data_long, 
  optimizer = optimizer_adam()
)

history_32 <- PC3_mod2 %>% fit(
  epochs = 500,
  early_stopping = TRUE, 
  sample_weight = data_long$Weight * post_probs[,2], 
  batch_size = 512, 
  callbacks = list(
    callback_early_stopping(
      monitor = "val_loss",
      patience = 10,
      restore_best_weights = TRUE
    )
  )
)

## TODO: continue here

# Analyze the results
TT <- length(ybin_list)
d <- ncol(ybin_list[[1]])
K <- 2

wide_X_df <- cbind(1, 1, 1, 1, 1, 1, X_df, 1:296, 1)
colnames(wide_X_df)[1:6] <- c(paste0("PC", 1:3), paste0("y", 1:3))
colnames(wide_X_df)[ncol(wide_X_df) - 1:0] <- c("Weight", "Time")

# Extract the means and probs over time
dist_pc1 <- get_distribution(PC1_mod, data = wide_X_df)
dist_pc21 <- get_distribution(PC2_mod1, data = wide_X_df)
dist_pc22 <- get_distribution(PC2_mod2, data = wide_X_df)
dist_pc31 <- get_distribution(PC3_mod1, data = wide_X_df)
dist_pc32 <- get_distribution(PC3_mod2, data = wide_X_df)

# Gather the mean predictions in principal space
mn_fit_pc <- abind(
  dist_pc1$submodules[[1]]$loc |> tf$squeeze(1L) |> as.matrix(), 
  cbind(
    dist_pc21$submodules[[1]]$loc |> as.matrix(), 
    dist_pc22$submodules[[1]]$loc |> as.matrix()
  ), 
  cbind(
    dist_pc31$submodules[[1]]$loc |> as.matrix(), 
    dist_pc32$submodules[[1]]$loc |> as.matrix()
  ), 
  along = 3
) |> aperm(c(1, 3, 2))

# Back-transform from principal space to original y space
mn_fit_list <- apply(mn_fit_pc, 3, function(x) {
  t(t(tcrossprod(x, y_pca$rotation)) + y_pca$center)
}, simplify = FALSE)

mn_fit <- do.call(abind, list(mn_fit_list, along = 3))

# Gather probability predictions
prob <- dist_pc1$submodules[[2]]$probs |> tf$squeeze(1) |> as.matrix()

# Gather variance estimates
pc1_vars <- as.numeric(dist_pc1$submodules[[1]]$scale[1,1,])^2
pc2_vars <- c(
  as.numeric(dist_pc21$submodules[[1]]$scale[1,])^2, 
  as.numeric(dist_pc22$submodules[[1]]$scale[1,])^2
)
pc3_vars <- c(
  as.numeric(dist_pc31$submodules[[1]]$scale[1,])^2, 
  as.numeric(dist_pc32$submodules[[1]]$scale[1,])^2
)

pc_clust1_cov <- diag(c(pc1_vars[1], pc2_vars[1], pc3_vars[1]))
pc_clust2_cov <- diag(c(pc1_vars[2], pc2_vars[2], pc3_vars[2]))
y_clust1_cov <- y_pca$rotation %*% pc_clust1_cov %*% t(y_pca$rotation)
y_clust2_cov <- y_pca$rotation %*% pc_clust2_cov %*% t(y_pca$rotation)

# Calculate error metrics

# Calculate RMSE (root mean l2 error)
resids <- mn_true - mn_fit
resid_dotprods <- apply(resids, 3, function(x) {
  crossprod(as.vector(t(x)))
})

rmse <- sqrt(resid_dotprods / TT)
prob_rmse <- sqrt(mean((prob1_true - prob[,1])^2))

# Covariance estimates would be a diagonal matrix
cov_fit <- abind(
  y_clust1_cov, y_clust2_cov, along = 3
)

cov_err <- sapply(1:K, function(k) {
  sqrt(sum((cov_true[,,k] - cov_fit[,,k])^2))
})

results <- data.frame(
  Model = "mixdistreg",
  Metric = rep(c("RMSE, Mean", "RMSE, Probability", "Frobenius Error, Covariance"), times = c(3, 1, 3)), 
  Cluster = c("1", "2", "Total", "1", "1", "2", "Total"),
  Value = c(rmse, sum(rmse), prob_rmse, cov_err, sum(cov_err))
)

write.csv(results, file.path("6_competitors", "metrics", "mixdistreg.csv"), row.names = FALSE)
