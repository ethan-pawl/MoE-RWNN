# TODO: rerun with weights
# TODO: double-check filepaths

library(randomForestSRC)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(gifski)
library(ggrepel)
library(tidyr)

rerun_clustering_and_matching <- FALSE
plot_animation <- FALSE
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

match_file <- file.path("5_competitors", "50_sim", "506_rf_wts", "sidcluster-matched.Rdata")
if(!file.exists(match_file) | rerun_clustering_and_matching) {
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

    # # Use the Hungarian algorithm to find the lowest total KL-divergence cluster match 
    # # (I don't use this because solution is pretty nonsensical)
    # solution <- RcppHungarian::HungarianSolver(D)
    # cur_samp_labels <- solution$pairs[,2][solution$pairs[,1]]

    # TODO: if time, try minimax instead of minimin
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

  save(clust_res, file = match_file)
} else {
  load(match_file)
}

# Check clustering

ylist_df <- lapply(seq_along(ybin_list), function(tt) {
  cur_ylist_df <- ybin_list[[tt]] |> as.data.frame()
  cur_ylist_df <- cbind(cur_ylist_df, countslist[[tt]])
  colnames(cur_ylist_df) <- c("y1", "y2", "y3", "Bin_Biomass")# , "Cluster Membership") # No ground truth here
  cur_ylist_df$time <- tt
  cur_ylist_df$cluster <- factor(clust_res[[tt]])

  cur_ylist_df 
}) |> bind_rows()

if(plot_animation) {
  
  y1_breaks <- seq(0, 2.5, 0.5)
  y2_breaks <- seq(0.75, 2.25, 0.5)
  y3_breaks <- seq(0.25, 2.75, 0.5)

  set1 <- brewer.pal(3, "Set1")

  for(tt in 1:TT) {
    cur_ylist_df <- subset(ylist_df, time == tt)

    p1 <- ggplot(cur_ylist_df) +
      geom_tile(aes(y1, y2, alpha = Bin_Biomass, fill = cluster)) + 
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
      scale_size_continuous(range = c(0, 4)) 
  
    p2 <- ggplot(cur_ylist_df) +
      geom_tile(aes(y2, y3, alpha = Bin_Biomass, fill = cluster)) + 
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
      scale_size_continuous(range = c(0, 4))

    plots_path <- file.path(
      "5_competitors", "50_sim", "506_rf_wts", "plots"
    )

    frames_path <- file.path(
      plots_path, "frames"
    )

    if(!dir.exists(frames_path)) dir.create(frames_path, recursive = TRUE)

    png(
      sprintf(
        file.path(frames_path, "frame_%04d.png"), 
        tt
      ), width = 1080, height = 540, res = 80, type = "cairo")

    grid.arrange(
      p1, p2,
      ncol = 2,
      top = sprintf("Time %d", tt)
    )

    graphics.off()
  }

  gifski(
    list.files(frames_path, full.names = TRUE),
    gif_file = file.path(plots_path, "animation.gif"),
    width = 1080,
    height = 540,
    delay = 1 / 10
  )
}

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

# TODO: try fast variant because this takes a long time
rf_res_file <- file.path("5_competitors", "50_sim", "506_rf_wts", "rf_res.RDS")
if(!file.exists(rf_res_file) | rerun_regression) {
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

  saveRDS(rf_res, rf_res_file)
} else {
  rf_res <- readRDS(rf_res_file)
}

# Analyze the results
TT <- length(ybin_list)
d <- ncol(ybin_list[[1]])
K <- length(rf_res)

mn_fit <- array(NA, c(TT, d, K))
Time <- rep(1, 296)
Weight <- rep(1, 296)
for(k in 1:K) {
  cur_preds <- predict(
    rf_res[[k]],
    newdata = cbind(X_df, Time, Weight)
  )

  # str(cur_preds, max.level = 1)
  # str(cur_preds$regrOutput, max.level = 1)
  # str(cur_preds$regrOutput$y1, max.level = 1)

  mn_fit[,,k] <- cbind(
    cur_preds$regrOutput$y1$predicted,
    cur_preds$regrOutput$y2$predicted,
    cur_preds$regrOutput$y3$predicted
  )
}

prob <- matrix(NA, TT, K)
for(tt in 1:TT) {
  prob[tt,] <- table(clust_res[[tt]]) / length(clust_res[[tt]])
}

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
  cur_ylist_df$cluster <- factor(clust_res[[tt]])

  cur_ylist_df 
}) |> bind_rows()

if(plot_regression_animation) {
  
  plots_path <- file.path(
    "5_competitors", "50_sim", "506_rf_wts", "plots_reg"
  )

  frames_path <- file.path(plots_path, "frames")

  if(!dir.exists(frames_path)) {
    dir.create(frames_path, recursive = TRUE)
  }

  y1_breaks <- seq(0, 2.5, 0.5)
  y2_breaks <- seq(0.75, 2.25, 0.5)
  y3_breaks <- seq(0.25, 2.75, 0.5)

  set1 <- brewer.pal(9, "Set1")

  prob_range <- range(means_probs_df$prob, na.rm = TRUE)

  for(tt in 1:TT) {
    cur_ylist_df <- subset(ylist_df, time == tt)
    cur_mp_df <- subset(means_probs_df, time == tt)

    p1 <- ggplot(cur_ylist_df) +
      geom_tile(aes(y1, y2, alpha = Bin_Biomass, fill = cluster)) + 
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
      )
  
    p2 <- ggplot(cur_ylist_df) +
      geom_tile(aes(y2, y3, alpha = Bin_Biomass, fill = cluster)) + 
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

  gifski(
    list.files(frames_path, full.names = TRUE),
    gif_file = file.path(plots_path, "animation.gif"),
    width = 1080,
    height = 540,
    delay = 1/10
  )
}

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

# Calculate RMSE (root mean l2 error)
resids <- mn_true - mn_fit
resid_dotprods <- apply(resids, 3, function(x) {
  crossprod(as.vector(t(x)))
})

rmse <- sqrt(resid_dotprods / TT)
rmse # [1] 0.2073716 0.2764358
sum(rmse) # [1] 0.4838074

prob_rmse <- sqrt(mean((prob1_true - prob[,1])^2))
prob_rmse # [1] 0.3337359

# FIXME: there are two different kinds of covariance estimates
#   1. Estimate of the error covariance structure
#   2. Covariance of the residuals
# The RF is non-probabilistic, so it does not provide (1)
# but I can always get (2)

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
cov_err # [1] 0.03367362 0.03066351

results <- data.frame(
  Model = "Random Forest",
  Metric = rep(c("RMSE, Mean", "RMSE, Probability", "Frobenius Error, Covariance"), times = c(3, 1, 3)), 
  Cluster = c("1", "2", "Total", "1", "1", "2", "Total"),
  Value = c(rmse, sum(rmse), prob_rmse, cov_err, sum(cov_err))
)

# write.csv(results, file.path("5_competitors", "metrics", "rf_wts.csv"), row.names = FALSE)
