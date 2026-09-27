library(parallel)
library(flowmix)
library(scales)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

simdata <- readRDS(
  file.path(
    "1_simulation", 
    "simdata",
    "simdata_imean_2_iprob_1_iint_5_dataSeed_1.RDS"
  )
)

# Data
ylist <- simdata$ybin_list
countslist <- simdata$countslist

# Estimation settings
maxdev <- diff(range(simdata$mean_spec)) / 2

X_dir <- file.path(
  "data", 
  "X_variations"
)

res_list <- mclapply(1:5, function(NNseed) {
  X <- readRDS(file.path(X_dir, paste0("X_pc_9_nh_70_seed_", NNseed, "_ofold_NA_ifold_NA.RDS")))
  load(file.path("3_seedsim", paste0(NNseed, "_settings.Rdata")))

  flowmix_once(
    ylist = ylist, 
    X = X, 
    countslist = countslist, 
    numclust = 2, 
    prob_lambda = prob_lambda, 
    mean_lambda = mean_lambda, 
    maxdev = maxdev, 
    seed = seed
  )
}, mc.cores = 5)

oos_data <- readRDS(
  file.path(
    "1_simulation", 
    "simdata", 
    "simdata_imean_2_iprob_1_iint_5_dataSeed_0.RDS"
  )
)

# Compute average of 5 models and standard errors

means <- array(NA, c(296, 5, 2))
means[,1,] <- res_list[[1]]$mn[,,2:1]
for(i in 2:5) {
  means[,i,] <- res_list[[i]]$mn[,1,]
}

avg_mns <- apply(means, MARGIN = c(1, 3), mean)
se_mns <- apply(means, MARGIN = c(1, 3), function(x) sd(x) / sqrt(5))

model_comp_df <- data.frame(
  `Time (t)` = rep(1:296, 2), 
  `Model Average` = as.numeric(avg_mns), 
  `Model Standard Error` = as.numeric(se_mns), 
  Cluster = rep(c("Cluster 1 Model-Averaged Mean", "Cluster 2 Model-Averaged Mean"), each = 296), 
  check.names = FALSE
)

model_comp_df <- model_comp_df |> 
  mutate(
    `Upper Uncertainty Bound` = `Model Average` + qt(1 - 0.05 / 2, 4) * `Model Standard Error`, 
    `Lower Uncertainty Bound` = `Model Average` - qt(1 - 0.05 / 2, 4) * `Model Standard Error`
  )

true_means_df <- data.frame(
  `Time (t)` = rep(1:296, 2),
  `Model Average` = c(oos_data$mean_spec, oos_data$pico_mu), 
  `Cluster` = rep(c("Cluster 1 True Mean", "Cluster 2 True Mean"), each = 296), 
  check.names = FALSE
)

brew_colors <- brewer.pal(n = 8, name = "Set1")

pdf(file.path("plots", "SuppFigure12.pdf"), 14.6, 8.5)
flowtrend::plot_1d(
  oos_data$ybin_list, oos_data$countslist, bin = TRUE, plot_band = FALSE
) + 
  geom_ribbon(
    data = model_comp_df, 
    aes(x = `Time (t)`, ymin = `Lower Uncertainty Bound`, ymax = `Upper Uncertainty Bound`, group = Cluster), 
    alpha = 0.85, 
    fill = brew_colors[4]
  ) + 
  geom_line(
    data = model_comp_df, 
    mapping = aes(x = `Time (t)`, y = `Model Average`, group = Cluster, color = Cluster, linetype = Cluster),
    lineend = "round", 
    linewidth = 1
  ) + 
  geom_line(
    data = true_means_df, 
    mapping = aes(x = `Time (t)`, y = `Model Average`, group = Cluster, color = Cluster, linetype = Cluster), 
    lineend = "round", 
    linewidth = 1
  ) + 
  scale_y_continuous(limits = c(0.6, 1.75)) + 
  scale_color_manual(
    name = "", 
    values = c(
      "Cluster 1 Model-Averaged Mean" = brew_colors[5], 
      "Cluster 2 Model-Averaged Mean" = brew_colors[2], 
      "Cluster 1 True Mean" = "black", 
      "Cluster 2 True Mean" = "black"
    ), 
    breaks = c(
      "Cluster 1 Model-Averaged Mean",
      "Cluster 2 Model-Averaged Mean",
      "Cluster 1 True Mean"
    ),
    labels = c(
      "Cluster 1 Model-Averaged Mean" = "Cluster 1 Model-Averaged Mean", 
      "Cluster 2 Model-Averaged Mean" = "Cluster 2 Model-Averaged Mean", 
      "Cluster 1 True Mean" = "True Mean"
    )
  ) + 
  scale_linetype_manual(
    name = "", 
    values = c(
      "Cluster 1 Model-Averaged Mean" = "solid", 
      "Cluster 2 Model-Averaged Mean" = "solid", 
      "Cluster 1 True Mean" = "twodash", 
      "Cluster 2 True Mean" = "twodash"
    ), 
    breaks = c(
      "Cluster 1 Model-Averaged Mean",
      "Cluster 2 Model-Averaged Mean",
      "Cluster 1 True Mean"
    ),
    labels = c(
      "Cluster 1 Model-Averaged Mean" = "Cluster 1 Model-Averaged Mean", 
      "Cluster 2 Model-Averaged Mean" = "Cluster 2 Model-Averaged Mean", 
      "Cluster 1 True Mean" = "True Mean"
    )
  ) + 
  guides(fill = "none") + 
  labs(title = "Model Averaged over Random Weight Resampling") + 
  scale_fill_gradient(
    low = alpha("grey92", 0.1), 
    high = alpha("black", 0.8)
  ) + 
  theme_gray() + 
  theme(
    legend.position = "bottom", 
    plot.title = element_text(size = 20, margin = margin(t = 5, b = 10)), 
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 18, margin = margin(r = 50)), 
    legend.key.size = unit(2, "lines"), 
    legend.margin = margin(0, 0, 0, 0)
  )
graphics.off()
