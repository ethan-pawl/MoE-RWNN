library(flowmix)
library(flowtrend)
library(ggplot2)
library(matrixStats)
library(gridExtra)
library(reshape2)
library(RColorBrewer)

destin <- file.path(
  "5_competitors", 
  "50_sim", 
  "507_linear", 
  "results"
)

res <- readRDS(file.path(destin, "3dsim_summary_linear.RDS"))
best <- res$bestres

in_sample <- readRDS(file.path("4_3dsim", "3dsimdata", "simdata_imean_2_iprob_1_iint_5_dataSeed_1.RDS"))

# Store ground truth
mn_true <- abind::abind(in_sample$mean_spec, in_sample$pico_mu, along = 3)
prob1_true <- in_sample$prob_spec
cov_true <- abind::abind(in_sample$clust1_cov, in_sample$clust2_cov, along = 3)

# Plot mean responses over time
response_comp_df <- data.frame(
  Quantity = c(
    rep(rep(c("Mean: Axis 1", "Mean: Axis 2", "Mean: Axis 3", "Relative Abundance"), each = 296), 2),
    rep(rep(c("Mean: Axis 1", "Mean: Axis 2", "Mean: Axis 3"), each = 296), 2)
  ),
  Linetype = rep(
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
  Value = c(
    as.vector(in_sample$mean_spec), 
    as.vector(in_sample$prob_spec), 
    as.vector(best$mn[,,1]), 
    as.vector(best$prob[,1]), 
    as.vector(in_sample$pico_mu),  
    as.vector(best$mn[,,2])
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

pdf(file.path("5_competitors", "50_sim", "507_linear", "plots", "mean_response_plot.pdf"), 15.125, 6.86)
print(p1)
graphics.off()

# Calculate RMSE (root mean l2 error)
mn_fit <- best$mn
resids <- mn_true - mn_fit 
resid_dotprods <- apply(resids, 3, function(x) {
  crossprod(as.vector(t(x)))
})

TT <- 296
rmse <- sqrt(resid_dotprods / TT)
rmse # [1] 0.04800998 0.01513985
sum(rmse) # [1] 0.06314983

prob <- best$prob
prob_rmse <- sqrt(mean((prob1_true - prob[,1])^2))
prob_rmse # [1] 0.01373893

K <- 2
cov_fit <- best$sigma |> aperm(c(2, 3, 1))
cov_err <- sapply(1:K, function(k) { 
  sqrt(sum((cov_true[,,k] - cov_fit[,,k])^2)) 
})
cov_err # [1] 0.0008149119 0.0006253325

results <- data.frame(
  Model = "Linear Model",
  Metric = rep(c("RMSE, Mean", "RMSE, Probability", "Frobenius Error, Covariance"), times = c(3, 1, 3)), 
  Cluster = c("1", "2", "Total", "1", "1", "2", "Total"),
  Value = c(rmse, sum(rmse), prob_rmse, cov_err, sum(cov_err))
)

# write.csv(results, file.path("5_competitors", "50_sim", "metrics", "linear_flowmix.csv"), row.names = FALSE)

# Plot clustering frame-by-frame
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

  png(file.path("5_competitions", "50_sim", "507_linear", "clustering", paste0(tt, ".png")), 10.5, 3.9, units = "in", res = 150, type = "cairo")
  print(grid.arrange(grobs = plist, ncol = 3))
  graphics.off()
}
