library(mixdistreg)
library(dplyr)

job_grid <- expand.grid(
  list(
    ialpha = 1:10, 
    ibeta = 1:10
  )
)

ijob <- commandArgs(trailingOnly = TRUE)

ialpha <- job_grid[ijob, "ialpha"]
ibeta <- job_grid[ijob, "ibeta"]

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

rm(datobj)

real_dataobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper.RDS"))
X <- real_dataobj$X
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

### 

y_mat <- data_long[,1:3]
y_pca <- prcomp(y_mat, center = TRUE, scale = FALSE)
y_pcs <- y_pca$x

data_long <- cbind(y_pcs, data_long)
colnames(data_long)[1:3] <- paste0("PC", 1:3)

###

l1_grid <- c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100)

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

PC1_formula_str <- paste0("~ 1 + deep_model(", paste(colnames(data_long[,7:(ncol(data_long) - 2)]), collapse = ","), ")")
PC1_formula <- as.formula(PC1_formula_str)

PC1_formula_mix_str <- paste0("~ 1 + deep_model_mix(", paste(colnames(data_long[,7:(ncol(data_long) - 2)]), collapse = ","), ")")
PC1_formula_mix <- as.formula(PC1_formula_mix_str)

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

cvscore_dir <- file.path("5_competitors", "50_sim", "503_mixdistreg", "cvscores")
if(!dir.exists(cvscore_dir)) dir.create(cvscore_dir)

cvscore <- min(history$metrics$val_loss)
save(ialpha, ibeta, cvscore, file = file.path(cvscore_dir, paste0(ialpha, "-", ibeta, "-cvscore.Rdata")))
