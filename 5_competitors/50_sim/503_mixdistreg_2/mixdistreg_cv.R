library(mixdistreg)
library(dplyr)
library(matrixStats)

nfold <- 5 
blocksize <- 20 

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
TT <- length(ybin_list)
rm(datobj)

folds <- flowmix::make_cv_folds(nfold = nfold, blocksize = blocksize, TT = TT)
l1_grid <- c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100)

job_grid <- expand.grid(
  list(
    ialpha = 1:length(l1_grid), 
    ibeta = 1:length(l1_grid),
    ifold = 1:nfold
  )
)

ijob <- commandArgs(trailingOnly = TRUE)[1]
ialpha <- job_grid[ijob, "ialpha"]
ibeta <- job_grid[ijob, "ibeta"]
ifold <- job_grid[ijob, "ifold"]

real_dataobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper.RDS"))
X <- real_dataobj$X
X_df <- as.data.frame(X[,3:ncol(X)])
rm(real_dataobj)

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

# Train-validation-outer validation split
load(file.path("data", "ofolds__ifolds__ifolds_inner_inds.Rdata"))

# Nested 5-fold will take forever and is probably unnecessary; instead, do an outer 5-fold 
# and an inner split which provides even temporal coverage
# -> hold a different inner fold each time so I'm early stopping on a variety of different 
#    times
out_sample_inds <- ofolds[[ifold]]
train_inds <- unlist(ifolds[[ifold]][-ifold])
val_inds <- ifolds[[ifold]][[ifold]]

y_mat_train <- subset(data_long, Time %in% train_inds)[,1:3]
y_mat_out <- subset(data_long, Time %in% out_sample_inds)[,1:3]
y_mat_val <- subset(data_long, Time %in% val_inds)[,1:3]

y_pca <- prcomp(y_mat_train, center = TRUE, scale = FALSE)
y_train_pcs <- y_pca$x
y_center <- y_pca$center 
y_val_pcs <- predict(y_pca, newdata = y_mat_val)
y_out_pcs <- predict(y_pca, newdata = y_mat_out)
# No need to transform y_mat_out because I will evaluate predictive error on the y scale anyway

data_long_train <- cbind(y_train_pcs, subset(data_long, Time %in% train_inds))
colnames(data_long_train)[1:3] <- paste0("PC", 1:3)

data_long_val <- subset(data_long, Time %in% val_inds)

deep_model <- function(x) {
  x %>%
    layer_dense(
      units = 64,
      activation = "relu",
      kernel_regularizer = regularizer_l1(l1_grid[ibeta])
    ) %>%
    layer_dense(
      units = 64,
      activation = "relu",
      kernel_regularizer = regularizer_l1(l1_grid[ibeta])
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
      kernel_regularizer = regularizer_l1(l1_grid[ialpha])
    ) %>%
    layer_dense(
      units = 64,
      activation = "relu",
      kernel_regularizer = regularizer_l1(l1_grid[ialpha])
    ) %>%
    layer_dense(
      units = 2,
      activation = "linear"
    )
}

PC1_formula_str <- paste0("~ 1 + deep_model(", paste(colnames(data_long_train[,7:(ncol(data_long_train) - 2)]), collapse = ","), ")")
PC1_formula <- as.formula(PC1_formula_str)

PC1_formula_mix_str <- paste0("~ 1 + deep_model_mix(", paste(colnames(data_long_train[,7:(ncol(data_long_train) - 2)]), collapse = ","), ")")
PC1_formula_mix <- as.formula(PC1_formula_mix_str)

PC1_mod <- mixdistreg(
  y_train_pcs[,1],
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
  data = data_long_train, 
  optimizer = optimizer_adam()
)

RhpcBLASctl::omp_set_num_threads(1L)
RhpcBLASctl::blas_set_num_threads(1L)

history <- PC1_mod %>% fit(
  epochs = 500,
  early_stopping = TRUE, 
  sample_weight = data_long_train$Weight, 
  batch_size = 512, 
  callbacks = list(
    callback_early_stopping(
      monitor = "val_loss",
      patience = 10,
      restore_best_weights = TRUE
    )
  ), 
  validation_data = list(
    data_long_val,
    y_val_pcs[,1]
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
      array(y_train_pcs[,1], dim = c(nrow(y_train_pcs),1, 1))
    ), 
    axis = 1L
  )
)

# Multiply and normalize (Bayes' Rule) to get responsibilities
post_probs <- pi_hat * dens
post_probs <- post_probs / rowSums(post_probs)

PC2_mod1 <- deepregression(
  y = y_train_pcs[,2],
  family = "normal", 
  list_of_formulas = list(
    loc = PC1_formula,
    scale = ~ 1
  ), 
  list_of_deep_models = list(deep_model = deep_model),
  data = data_long_train, 
  optimizer = optimizer_adam()
)

history_21 <- PC2_mod1 %>% fit(
  epochs = 500,
  early_stopping = TRUE, 
  sample_weight = data_long_train$Weight * post_probs[,1], 
  batch_size = 512, 
  callbacks = list(
    callback_early_stopping(
      monitor = "val_loss",
      patience = 10,
      restore_best_weights = TRUE
    )
  ), 
  validation_data = list(
    data_long_val,
    y_val_pcs[,2]
  )
)

PC2_mod2 <- deepregression(
  y = y_train_pcs[,2],
  family = "normal", 
  list_of_formulas = list(
    loc = PC1_formula,
    scale = ~ 1
  ), 
  list_of_deep_models = list(deep_model = deep_model),
  data = data_long_train, 
  optimizer = optimizer_adam()
)

history_22 <- PC2_mod2 %>% fit(
  epochs = 500,
  early_stopping = TRUE, 
  sample_weight = data_long_train$Weight * post_probs[,2], 
  batch_size = 512, 
  callbacks = list(
    callback_early_stopping(
      monitor = "val_loss",
      patience = 10,
      restore_best_weights = TRUE
    )
  ), 
  validation_data = list(
    data_long_val,
    y_val_pcs[,2]
  )
)

PC3_mod1 <- deepregression(
  y = y_train_pcs[,3],
  family = "normal", 
  list_of_formulas = list(
    loc = PC1_formula,
    scale = ~ 1
  ), 
  list_of_deep_models = list(deep_model = deep_model),
  data = data_long_train, 
  optimizer = optimizer_adam()
)

history_31 <- PC3_mod1 %>% fit(
  epochs = 500,
  early_stopping = TRUE, 
  sample_weight = data_long_train$Weight * post_probs[,1], 
  batch_size = 512, 
  callbacks = list(
    callback_early_stopping(
      monitor = "val_loss",
      patience = 10,
      restore_best_weights = TRUE
    )
  ), 
  validation_data = list(
    data_long_val,
    y_val_pcs[,3]
  )
)

PC3_mod2 <- deepregression(
  y = y_train_pcs[,3],
  family = "normal", 
  list_of_formulas = list(
    loc = PC1_formula,
    scale = ~ 1
  ), 
  list_of_deep_models = list(deep_model = deep_model),
  data = data_long_train, 
  optimizer = optimizer_adam()
)

history_32 <- PC3_mod2 %>% fit(
  epochs = 500,
  early_stopping = TRUE, 
  sample_weight = data_long_train$Weight * post_probs[,2], 
  batch_size = 512, 
  callbacks = list(
    callback_early_stopping(
      monitor = "val_loss",
      patience = 10,
      restore_best_weights = TRUE
    )
  ), 
  validation_data = list(
    data_long_val,
    y_val_pcs[,3]
  )
)

# Analyze the results
TT <- length(ybin_list)
d <- ncol(ybin_list[[1]])
K <- 2

data_long_out <- cbind(y_out_pcs, subset(data_long, Time %in% out_sample_inds))

# Extract the means and probs over time
dist_pc1 <- get_distribution(PC1_mod, data = data_long_out)
dist_pc21 <- get_distribution(PC2_mod1, data = data_long_out)
dist_pc22 <- get_distribution(PC2_mod2, data = data_long_out)
dist_pc31 <- get_distribution(PC3_mod1, data = data_long_out)
dist_pc32 <- get_distribution(PC3_mod2, data = data_long_out)

# Extract the likelihood contributions
n_out <- nrow(data_long_out)

# Gating probabilities
pi_hat <- as.matrix(
  tf$squeeze(
    dist_pc1$submodules[[2]]$probs,
    axis = 1L
  )
)

K <- ncol(pi_hat)

# Component densities
logdens_PC1 <- as.matrix(
  tf$squeeze(
    dist_pc1$submodules[[1]]$log_prob(
      array(y_out_pcs[,1], dim = c(n_out, 1, 1))
    ),
    axis = 1L
  )
)

logdens_PC2 <- cbind(
  dist_pc21$log_prob(
    array(y_out_pcs[,2], dim = c(n_out,1))
  ),
  dist_pc22$log_prob(
    array(y_out_pcs[,2], dim = c(n_out,1))
  )
)

logdens_PC3 <- cbind(
  dist_pc31$log_prob(
    array(y_out_pcs[,3], dim = c(n_out,1))
  ),
  dist_pc32$log_prob(
    array(y_out_pcs[,3], dim = c(n_out,1))
  )
)

logdens_PC2 <- as.matrix(logdens_PC2)
logdens_PC3 <- as.matrix(logdens_PC3)

# Log likelihood contributions within each component
log_component_density <- 
  log(pi_hat) + logdens_PC1 + logdens_PC2 + logdens_PC3

# Marginalize out cluster membership to get marginal likelihood contributions
log_mix_density <- rowLogSumExps(log_component_density)

# Calculate weighted negative log-likelihood (out of sample)
weights <- data_long_out$Weight
WNLL <- -sum(weights * log_mix_density) / sum(weights)

cvscore_dir <- file.path("5_competitors", "50_sim", "503_mixdistreg", "cvscores")
if(!dir.exists(cvscore_dir)) dir.create(cvscore_dir)

cvscore <- WNLL
save(ialpha, ibeta, ifold, cvscore, file = file.path(cvscore_dir, paste(ialpha, ibeta, ifold, "cvscore.Rdata", sep = "-")))
