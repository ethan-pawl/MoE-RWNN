plot_animation <- FALSE
plot_mean_responses <- FALSE

library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(ggplot2)
# library(magick)
library(ggrepel)
library(gifski)
library(tidyr)

per_cluster_res <- readRDS(file.path("5_competitors", "50_sim", "501_gam", "bam_refit.RDS"))

datobj <- readRDS(
  file.path(
    "4_3dsim", 
    "3dsimdata", 
    "simdata_imean_2_iprob_1_iint_5_dataSeed_1.RDS"
  )
)

# str(datobj, max.level = 1)

ybin_list <- datobj$ybin_list
countslist <- datobj$countslist

TT <- length(ybin_list)
gt_df <- lapply(1:TT, function(tt) {
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
# 296 x 3 x 2

rm(datobj)

real_dataobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper.RDS"))
X <- real_dataobj$X
rm(real_dataobj)

# str(per_cluster_res, max.level = 1)
# str(per_cluster_res$y1_res, max.level = 1)
# str(per_cluster_res$y1_res[[1]], max.level = 1)

library(mgcv)

TT <- length(ybin_list)
d <- length(per_cluster_res)
K <- length(per_cluster_res[[1]])

mn_fit <- array(NA, c(TT, d, K))

X_df <- as.data.frame(X[,3:ncol(X)])

for(j in 1:d) {
  for(k in 1:K) {
    mn_fit[,j,k] <- predict(
      per_cluster_res[[j]][[k]],
      newdata = X_df, 
      type = "response"
    )
  }
}

load(file.path("5_competitors", "50_sim", "500_k_means", "best_kmeans__scores-matched.Rdata"))

prob <- matrix(NA, TT, K)
for(tt in 1:TT) {
  prob[tt,] <- best_kmeans[[tt]]$size / sum(best_kmeans[[tt]]$size)
}

# TODO: plot clustering with means over time
# TODO: plot mean responses over time, comparing with 
# ground truth (use preexisting code)

threeD_sim_summary <- readRDS("/home/ethan/00_Cyto/MoE-RWNN/4_3dsim/results/3dsim_summary.RDS")

# str(threeD_sim_summary$bestres, max.level = 1)

ybin_list <- lapply(ybin_list, function(cur_y) {
  colnames(cur_y) <- paste0("y", 1:3)
  cur_y
})
 
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
  cur_ylist_df$cluster <- factor(best_kmeans[[tt]]$cluster)

  cur_ylist_df 
}) |> bind_rows()

if(plot_animation) {
  
  plots_path <- file.path(
    "5_competitors", "50_sim", "501_gam", "plots_refit"
  )

  frames_path <- file.path(plots_path, "frames")

  if(!dir.exists(frames_path)) {
    dir.create(frames_path, recursive = TRUE)
  }

  y1_breaks <- seq(0, 2.5, 0.5)
  y2_breaks <- seq(0.75, 2.25, 0.5)
  y3_breaks <- seq(0.25, 2.75, 0.5)

  set1 <- brewer.pal(9, "Set1")

  # TODO: maybe remove this for compatibility
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
resids <- mn_true - mn_fit[,,2:1] # LABELS SWITCHED (OTHER PLOTS DON'T ACCOUNT FOR THIS)
resid_dotprods <- apply(resids, 3, function(x) {
  crossprod(as.vector(t(x)))
})

rmse <- sqrt(resid_dotprods / TT)
rmse # [1] 0.2233791 0.2681352
sum(rmse) # [1] 0.4915143

prob_rmse <- sqrt(mean((prob1_true - prob[,2])^2))
prob_rmse # [1] 0.305452

# Covariance estimates would be a diagonal matrix
cov_fit <- sapply(1:K, function(k) {
  sapply(1:d, function(j) {
    per_cluster_res[[j]][[k]]$sig2
  }) |> diag()
}, simplify = "array")

cov_err <- sapply(1:K, function(k) {
  sqrt(sum((cov_true[,,k] - cov_fit[,,if(k == 1) 2 else 1])^2))
})
cov_err
# [1] 0.02366263 0.03586512

results <- data.frame(
  Model = "GAM",
  Metric = rep(c("RMSE, Mean", "RMSE, Probability", "Frobenius Error, Covariance"), times = c(3, 1, 3)), 
  Cluster = c("1", "2", "Total", "1", "1", "2", "Total"),
  Value = c(rmse, sum(rmse), prob_rmse, cov_err, sum(cov_err))
)

# write.csv(results, file.path("5_competitors", "metrics", "gam.csv"), row.names = FALSE)