library(mixdistreg)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(gifski)
library(ggrepel)
library(tidyr)
library(ellipse)

plot_animation <- TRUE
rerun_regression <- FALSE
plot_regression_animation <- TRUE
plot_mean_responses <- TRUE

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

gt_df <- lapply(seq_along(ybin_list), function(tt) {
  rbind(
    data.frame( # Cluster 1
      mean_1 = datobj$mean_spec[tt,1],
      mean_2 = datobj$mean_spec[tt,2],
      mean_3 = datobj$mean_spec[tt,3],
      prob = datobj$prob_spec[tt,],
      cluster = "Cluster 1 Truth",
      time = tt
    ), 
    data.frame( # Cluster 2
      mean_1 = datobj$pico_mu[tt,1],
      mean_2 = datobj$pico_mu[tt,2],
      mean_3 = datobj$pico_mu[tt,3],
      prob = 1 - datobj$prob_spec[tt,],
      cluster = "Cluster 2 Truth",
      time = tt
    )
  )
}) |> bind_rows()

mn_true <- abind::abind(datobj$mean_spec, datobj$pico_mu, along = 3)
prob1_true <- datobj$prob_spec
cov_true <- abind::abind(datobj$clust1_cov, datobj$clust2_cov, along = 3)

rm(datobj)

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

### 

y_mat <- data_long[,1:3]
y_pca <- prcomp(y_mat, center = TRUE, scale = FALSE)
y_pcs <- y_pca$x

data_long <- cbind(y_pcs, data_long)
colnames(data_long)[1:3] <- paste0("PC", 1:3)

###

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

# dim(pi_hat)
# head(pi_hat)
# apply(pi_hat, 2, range)
# apply(pi_hat, 2, sd)

# logits <- as.matrix(
#   tf$squeeze(
#     dist_dr$submodules[[2]]$logits_parameter(),
#     axis = 1L
#   )
# )

# apply(logits, 2, range)
# apply(logits, 2, sd)

# Get posterior probabilities of cluster membership

# mu <- dist_obj$submodules[[1]]$loc |> tf$squeeze(1L) |> as.matrix()

# dim(mu)
# head(mu)
# apply(mu, 2, range)
# apply(mu, 2, sd)

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

mn_fit_pc <- abind::abind(
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

# Invert the PCA transformation
mn_fit_list <- apply(mn_fit_pc, 3, function(x) {
  t(t(tcrossprod(x, y_pca$rotation)) + y_pca$center)
}, simplify = FALSE)

mn_fit <- do.call(abind::abind, list(mn_fit_list, along = 3))

prob <- dist_pc1$submodules[[2]]$probs |> tf$squeeze(1) |> as.matrix()

means_probs_df <- lapply(1:TT, function(tt) {
  data.frame(
    # mean_1 = best_kmeans[[tt]]$centers[,1],
    # mean_2 = best_kmeans[[tt]]$centers[,2],
    # mean_3 = best_kmeans[[tt]]$centers[,3],
    mean_1 = mn_fit[tt,1,],
    mean_2 = mn_fit[tt,2,],
    mean_3 = mn_fit[tt,3,],
    prob = prob[tt,],
    cluster = 1:K,
    time = tt
  )
}) |> bind_rows()

ylist_df <- lapply(1:TT, function(tt) {
  cur_ylist_df <- ybin_list[[tt]] |> as.data.frame()
  cur_ylist_df <- cbind(cur_ylist_df, countslist[[tt]])
  colnames(cur_ylist_df) <- c("y1", "y2", "y3", "Bin_Biomass")# , "Cluster Membership") # No ground truth here
  cur_ylist_df$time <- tt

  cur_ylist_df 
}) |> bind_rows()

# FIXME: also need the ellipse_df

set1 <- brewer.pal(9, "Set1")
plots_path <- file.path(
  "5_competitors", "50_sim", "503_mixdistreg", "plots_reg"
)
if(!dir.exists(plots_path)) dir.create(plots_path)

if(plot_mean_responses) {
  means_probs_df$cluster <- paste0("Cluster ", means_probs_df$cluster, " Estimate")
  means_probs_df$Type <- means_probs_df$cluster 
  gt_df$Type <- "Truth"

  full_df <- rbind(means_probs_df, gt_df)

  full_df_long <- full_df |> 
    rename(
      `Mean: Axis 1` = mean_1,
      `Mean: Axis 2` = mean_2,
      `Mean: Axis 3` = mean_3,
      `Relative Abundance` = prob, 
      Cluster = cluster, 
      `Time (t)` = time
    ) |> 
    pivot_longer(
      c(starts_with("Mean"), "Relative Abundance"), 
      names_to = "Data Type", 
      values_to = "Value"
    )
  
  # Plot mean responses over time, estimates and ground truth

  leg_breaks <- c(
    "Cluster 1 Estimate", 
    "Cluster 2 Estimate", 
    "Cluster 1 Truth"
  )

  leg_labels <- c(
    "Cluster 1 Estimate" = "Cluster 1 Estimate", 
    "Cluster 2 Estimate" = "Cluster 2 Estimate", 
    "Cluster 1 Truth" = "Truth"
  )

  mean_response_plot <- ggplot(full_df_long |> filter(!(`Data Type` == "Relative Abundance" & Cluster %in% paste0("Cluster 2 ", c("Truth", "Estimate"))))) + 
    geom_line(aes(`Time (t)`, Value, linetype = Cluster, color = Cluster)) + 
    facet_wrap(~ `Data Type`) + 
    theme_gray() + 
    scale_color_manual(
      name = "",
      values = c(
        "Cluster 1 Estimate" = set1[5], 
        "Cluster 2 Estimate" = set1[2], 
        "Cluster 1 Truth" = "black", 
        "Cluster 2 Truth" = "black"
      ), 
      labels = leg_labels,
      breaks = leg_breaks
    ) + 
    scale_linetype_manual(
      name = "",
      values = c(
        "Cluster 1 Estimate" = "solid", 
        "Cluster 2 Estimate" = "solid", 
        "Cluster 1 Truth" = "twodash", 
        "Cluster 2 Truth" = "twodash"
      ), 
      labels = leg_labels,
      breaks = leg_breaks
    ) + 
    theme(legend.position = "bottom")
  
  ggsave(file.path(plots_path, "mean_response_plot.pdf"), mean_response_plot, width = 15.25, height = 6.75, 
    units = "in"
  )
}

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

# save(
#   mn_fit_pc, mn_fit, prob, pc_clust1_cov, pc_clust2_cov, 
#   y_clust1_cov, y_clust2_cov, 
#   file = file.path("5_competitors", "50_sim", "503_mixdistreg", "mixdistreg_results.Rdata")
# )

# load(file.path("5_competitors", "50_sim", "503_mixdistreg", "mixdistreg_results.Rdata"))

ellipse_df_dim1 <- lapply(1:TT, function(tt) {
  rbind(
    ellipse(y_clust1_cov, centre = mn_fit[tt,1:2,1]),
    ellipse(y_clust2_cov, centre = mn_fit[tt,1:2,2])
  )
}) %>% do.call(rbind, .)

ellipse_df_dim1 <- as.data.frame(ellipse_df_dim1)
ellipse_df_dim1$cluster <- rep(1:K, each = 100) |> rep(TT)
ellipse_df_dim1$time <- rep(1:TT, each = 100 * K)

ellipse_df_dim2 <- lapply(1:TT, function(tt) {
  rbind(
    ellipse(y_clust1_cov, centre = mn_fit[tt,2:3,1], which = c(2, 3)),
    ellipse(y_clust2_cov, centre = mn_fit[tt,2:3,2], which = c(2, 3))
  )
}) %>% do.call(rbind, .)

ellipse_df_dim2 <- as.data.frame(ellipse_df_dim2)
ellipse_df_dim2$cluster <- rep(1:K, each = 100) |> rep(TT)
ellipse_df_dim2$time <- rep(1:TT, each = 100 * K)

if(plot_regression_animation) {  

  frames_path <- file.path(plots_path, "frames")

  if(!dir.exists(frames_path)) {
    dir.create(frames_path, recursive = TRUE)
  }

  y1_breaks <- seq(0, 2.5, 0.5)
  y2_breaks <- seq(0.75, 2.25, 0.5)
  y3_breaks <- seq(0.25, 2.75, 0.5)


  prob_range <- range(means_probs_df$prob, na.rm = TRUE)

  for(tt in 1:TT) {
    cur_ylist_df <- subset(ylist_df, time == tt)
    cur_mp_df <- subset(means_probs_df, time == tt)
    cur_ellipse_df_dim1 <- subset(ellipse_df_dim1, time == tt)
    cur_ellipse_df_dim2 <- subset(ellipse_df_dim2, time == tt)

    p1 <- ggplot(cur_ylist_df) +
      geom_tile(aes(y1, y2, alpha = Bin_Biomass)) + 
      scale_x_continuous(breaks = y1_breaks) + 
      scale_y_continuous(breaks = y2_breaks) + 
      coord_cartesian(xlim = range(y1_breaks), ylim = range(y2_breaks)) + 
      scale_fill_manual(values = c(set1, "black")) + 
      scale_alpha_continuous(range = c(0.2, 1)) + 
      theme_gray(base_family = "sans") + 
      theme(axis.title = element_text(size = 8), plot.title = element_text(size = 10), 
        axis.text = element_text(size = 6), legend.position = "none", 
        axis.line = element_line(linewidth = 0.25)
      ) + 
      geom_point(
        data = cur_mp_df, 
        mapping = aes(x = mean_1, y = mean_2, size = prob)
      ) + 
      scale_size_continuous(range = c(0, 4)) + 
      geom_text_repel(data = cur_mp_df, 
        mapping = aes(mean_1, mean_2, label = cluster), size = 3, box.padding = 0.25, 
        min.segment.length = 0.25, segment.size = 0.25
      ) + 
      geom_path(
        data = cur_ellipse_df_dim1, 
        mapping = aes(y1, y2, group = cluster), 
        linetype = "dashed", linewidth = 0.25
      )
  
    p2 <- ggplot(cur_ylist_df) +
      geom_tile(aes(y2, y3, alpha = Bin_Biomass)) + 
      scale_x_continuous(breaks = y2_breaks) + 
      scale_y_continuous(breaks = y3_breaks) + 
      coord_cartesian(xlim = range(y2_breaks), ylim = range(y3_breaks)) + 
      scale_fill_manual(values = c(set1, "black")) + 
      scale_alpha_continuous(range = c(0.2, 1)) + 
      theme_gray(base_family = "sans") + 
      theme(axis.title = element_text(size = 8), plot.title = element_text(size = 10), 
        axis.text = element_text(size = 6), legend.position = "none", 
        axis.line = element_line(linewidth = 0.25)
      ) + 
      geom_point(
        data = cur_mp_df, 
        mapping = aes(x = mean_2, y = mean_3, size = prob)
      ) +
      scale_size_continuous(range = c(0, 4)) + 
      geom_text_repel(data = cur_mp_df, 
        mapping = aes(mean_2, mean_3, label = cluster), size = 3, box.padding = 0.25, 
        min.segment.length = 0.25, segment.size = 0.25
      ) + 
      geom_path(
        data = cur_ellipse_df_dim2, 
        mapping = aes(y2, y3, group = cluster), 
        linetype = "dashed", linewidth = 0.25
      )

    png(
      sprintf(file.path(frames_path, "frame_%04d.png"), tt), 
      width = 1080, height = 540, res = 80, type = "cairo"
    )

    grid.arrange(
      p1, p2,
      ncol = 2,
      top = sprintf("Time %d", tt)
    )

    dev.off()
  }

  # memory cache exhausted; try to rerun
  # imgs <- image_read(list.files(file.path("5_competitors", "50_sim", "501_gam", "plots", "frames"), full.names = TRUE))
  # gif <- image_animate(imgs, fps = 10)
  # image_write(gif, "animation.gif")

  gifski(
    list.files(frames_path, full.names = TRUE),
    gif_file = file.path(plots_path, "animation.gif"),
    width = 1080,
    height = 540,
    delay = 1/10
  )
}

# Calculate RMSE (root mean l2 error)
resids <- mn_true - mn_fit
resid_dotprods <- apply(resids, 3, function(x) {
  crossprod(as.vector(t(x)))
})

rmse <- sqrt(resid_dotprods / TT)
rmse # [1] 0.2230055 0.1223208
sum(rmse) # [1] 0.3453264

prob_rmse <- sqrt(mean((prob1_true - prob[,1])^2))
prob_rmse # [1] 0.3582185

# Covariance estimates would be a diagonal matrix
cov_fit <- abind::abind(
  y_clust1_cov, y_clust2_cov, along = 3
)

cov_err <- sapply(1:K, function(k) {
  sqrt(sum((cov_true[,,k] - cov_fit[,,k])^2))
})
cov_err
# [1] 0.02400551 0.04493897

results <- data.frame(
  Model = "mixdistreg",
  Metric = rep(c("RMSE, Mean", "RMSE, Probability", "Frobenius Error, Covariance"), times = c(3, 1, 3)), 
  Cluster = c("1", "2", "Total", "1", "1", "2", "Total"),
  Value = c(rmse, sum(rmse), prob_rmse, cov_err, sum(cov_err))
)

write.csv(results, file.path("5_competitors", "50_sim", "metrics", "mixdistreg.csv"), row.names = FALSE)
