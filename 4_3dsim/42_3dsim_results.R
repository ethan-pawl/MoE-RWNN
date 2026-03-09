library(flowmix)
library(flowtrend)
library(ggplot2)
library(matrixStats)
library(gridExtra)
library(reshape2)
library(RColorBrewer)

res <- readRDS(file.path("4_3dsim", "results", "3dsim_summary.RDS"))
best <- res$bestres

in_sample <- readRDS(file.path("4_3dsim", "3dsimdata", "simdata_imean_2_iprob_1_iint_5_dataSeed_1.RDS"))

mins <- matrix(NA, 296, 3)
maxs <- matrix(NA, 296, 3)
for(tt in 1:296) {
  colnames(in_sample$ybin_list[[tt]]) <- paste0("y", 1:3)
  mins[tt,] <- colMins(in_sample$ybin_list[[tt]])
  maxs[tt,] <- colMaxs(in_sample$ybin_list[[tt]])
}

mins <- colMins(mins)
maxs <- colMaxs(maxs)

mins <- floor(mins / 0.25) * 0.25
maxs <- ceiling(maxs / 0.25) * 0.25

for(tt in 1:296) {
  plist <- plot_3d(
    in_sample$ybin_list, best, tt, in_sample$countslist, return_list_of_plots = TRUE, 
    labels = c("2", "1")
  )

  for(j in 1:3) {
    dim1 <- j 
    dim2 <- (j %% 3) + 1

    plist[[j]] <- plist[[j]] + 
      scale_x_continuous(limits = c(mins[dim1], maxs[dim1]), breaks = seq(mins[dim1], maxs[dim1], by = 0.5)) + 
      scale_y_continuous(limits = c(mins[dim2], maxs[dim2]), breaks = seq(mins[dim2], maxs[dim2], by = 0.5))
  }

  png(file.path("4_3dsim", "clustering", paste0(tt, ".png")), 10.5, 3.9, units = "in", res = 150, type = "cairo")
  print(grid.arrange(grobs = plist, ncol = 3))
  graphics.off()
}

figure_tt <- 192
plist <- plot_3d(
  in_sample$ybin_list, best, figure_tt, in_sample$countslist, return_list_of_plots = TRUE, 
  labels = c("2", "1")
)

for(j in 1:3) {
  dim1 <- j 
  dim2 <- (j %% 3) + 1

  plist[[j]] <- plist[[j]] + 
    scale_x_continuous(limits = c(mins[dim1], maxs[dim1]), breaks = seq(mins[dim1], maxs[dim1], by = 0.5)) + 
    scale_y_continuous(limits = c(mins[dim2], maxs[dim2]), breaks = seq(mins[dim2], maxs[dim2], by = 0.5))
}

pdf(file.path("plots", "Figure15_01.pdf"), 7, 3.9)
print(grid.arrange(grobs = plist[1:2], ncol = 2))
graphics.off()

figure_2_tt <- 67
plist2 <- plot_3d(
  in_sample$ybin_list, best, figure_2_tt, in_sample$countslist, return_list_of_plots = TRUE, 
  labels = c("2", "1")
)

for(j in 1:3) {
  dim1 <- j 
  dim2 <- (j %% 3) + 1

  plist2[[j]] <- plist2[[j]] + 
    scale_x_continuous(limits = c(mins[dim1], maxs[dim1]), breaks = seq(mins[dim1], maxs[dim1], by = 0.5)) + 
    scale_y_continuous(limits = c(mins[dim2], maxs[dim2]), breaks = seq(mins[dim2], maxs[dim2], by = 0.5))
}


pdf(file.path("plots", "Figure15_02.pdf"), 7, 3.9)
print(grid.arrange(grobs = plist2[1:2], ncol = 2))
graphics.off()

plist_comp <- list()
plist_comp[1:2] <- plist2[1:2]
plist_comp[3:4] <- plist[1:2]

pdf(file.path("plots", "Figure15_comp.pdf"), 7, 7.8)
print(grid.arrange(grobs = plist_comp, ncol = 2))
graphics.off()

# Mean, prob, and covariance error metrics

clust1_mn_rmse <- (in_sample$mean_spec - best$mn[,,2])^2 |> 
  rowSums() |>
  mean() |> 
  sqrt()

clust2_mn_rmse <- (in_sample$pico_mu - best$mn[,,1])^2 |> 
  rowSums() |>
  mean() |> 
  sqrt()

prob_rmse <- (in_sample$prob_spec - best$prob[,2])^2 |> 
  mean() |> 
  sqrt()

Sigma1_rfe <- sqrt(norm(in_sample$clust1_cov - best$sigma[2,,], "F"))
Sigma2_rfe <- sqrt(norm(in_sample$clust2_cov - best$sigma[1,,], "F"))
Sigma1_rmfe <- sqrt(norm(in_sample$clust1_cov - best$sigma[2,,], "F") / 9)
Sigma2_rmfe <- sqrt(norm(in_sample$clust2_cov - best$sigma[1,,], "F") / 9)

clust1_mn_rmse # 0.01043708
clust2_mn_rmse # 0.01565043
prob_rmse # 0.01036922
Sigma1_rfe # 0.02466326
Sigma2_rfe # 0.02175855
Sigma1_rmfe # 0.008221088
Sigma2_rmfe # 0.00725285

# Compare with ranges of data

# Range of cluster 1 mean is 
best$maxdev * 2 # 0.4295578

# Range of probability is 
diff(range(in_sample$prob_spec)) # 0.9280987

# Rt Frobenius norms of covariance matrices is 
sqrt(norm(in_sample$clust1_cov, "F")) # 0.2835307
sqrt(norm(in_sample$clust2_cov, "F")) # 0.2907639

# Scaled by # of components
sqrt(norm(in_sample$clust1_cov, "F") / 9) # 0.09451024
sqrt(norm(in_sample$clust2_cov, "F") / 9) # 0.09692129

# Just for cluster 1
response_comp_df <- data.frame(
  Quantity = c(
    rep(rep(c("Mean: Axis 1", "Mean: Axis 2", "Mean: Axis 3", "Relative Abundance"), each = 296), 2),
    rep(rep(c("Mean: Axis 1", "Mean: Axis 2", "Mean: Axis 3"), each = 296), 2)
  ),
  Linetype = rep( # 4144
    c(
      "Truth", 
      "Estimate", 
      "Truth", 
      "Estimate"
    ), times = c(296*4, 296*4, 296*3, 296*3)
  ), 
  Cluster = rep(
    c(
      "Cluster 1 Truth",
      "Cluster 1 Estimate", 
      "Cluster 2 Truth",
      "Cluster 2 Estimate"
    ), times = c(296*4, 296*4, 296*3, 296*3)
  ),
  Color = rep(
    c(
      "Truth",
      "Cluster 1 Estimate", 
      "Truth",
      "Cluster 2 Estimate"
    ), times = c(296*4, 296*4, 296*3, 296*3)
  ),
  Value = c(# 4144
    as.vector(in_sample$mean_spec), 
    as.vector(in_sample$prob_spec), 
    as.vector(best$mn[,,2]), 
    as.vector(best$prob[,2]), 
    as.vector(in_sample$pico_mu),  
    as.vector(best$mn[,,1])
  ), 
  `Time (t)` = rep(1:296, 14),
  check.names = FALSE
)

# Colors

brew_colors <- brewer.pal(n = 8, name = "Set1")

named_colors <- c("Cluster 1 Estimate" = brew_colors[5],
                  "Cluster 2 Estimate" = brew_colors[2],
                  "Cluster 1 Truth" = "black", 
                  "Cluster 2 Truth" = "black")

named_linetypes <- c("Estimate" = "solid",
                     "Truth" = "twodash")

labels <- c("Cluster 1 Estimate" = "Cluster 1 Estimate",
            "Cluster 2 Estimate" = "Cluster 2 Estimate",
            "Cluster 1 Truth" = "Truth", 
            "Cluster 2 Truth" = "Truth")

p1 <- ggplot(response_comp_df, aes(`Time (t)`, Value, group = Cluster, color = Cluster, linetype = Cluster)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~ Quantity, scales = "fixed") + 
  scale_color_manual(
    name = "", 
    values = c(
      "Cluster 1 Estimate" = brew_colors[5], 
      "Cluster 2 Estimate" = brew_colors[2],
      "Cluster 1 Truth" = "black", 
      "Cluster 2 Truth" = "black"
    ), 
    breaks = c(
      "Cluster 1 Estimate", 
      "Cluster 2 Estimate", 
      "Cluster 1 Truth"
    ), 
    labels = c(
      "Cluster 1 Estimate", 
      "Cluster 2 Estimate", 
      "Truth"
    )
  ) + 
  scale_linetype_manual(
    name = "", 
    values = c(
      "Cluster 1 Estimate" = "solid", 
      "Cluster 2 Estimate" = "solid",
      "Cluster 1 Truth" = "twodash", 
      "Cluster 2 Truth" = "twodash"
    ), 
    breaks = c(
      "Cluster 1 Estimate", 
      "Cluster 2 Estimate", 
      "Cluster 1 Truth"
    ), 
    labels = c(
      "Cluster 1 Estimate", 
      "Cluster 2 Estimate", 
      "Truth"
    )
  ) + 
  theme(
    legend.position = "bottom", 
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 16, margin = margin(r = 50)), 
    legend.key.size = unit(2, "lines"), 
    strip.text = element_text(size = 14)
  )

pdf(file.path("plots", "Figure15_03.pdf"), 15.125, 6.86)
print(p1)
graphics.off()