library(flowmix)
library(parallel)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)

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

res_list <- mclapply(c(35, 70, 105, 140, 175), function(nh) {
  if(nh == 70) {
    res <- readRDS(file.path("1_simulation", "results", "2-1-5", "2-1-5-70-1_fit_model.RDS"))
    return(res)
  } else {
    X <- readRDS(file.path(X_dir, paste0("X_pc_9_nh_", nh, "_seed_1_ofold_NA_ifold_NA.RDS")))
    load(file.path("5_widthsim", paste0(nh, "_settings.Rdata")))
  
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
  }
}, mc.cores = 5)

# Fix label switching in 140 & 175
for(i in 4:5) {
  res_list[[i]]$mn <- res_list[[i]]$mn[,,2:1,drop = FALSE]
  res_list[[i]]$prob <- res_list[[i]]$prob[,2:1,drop = FALSE]
  res_list[[i]]$sigma <- res_list[[i]]$sigma[2:1,,,drop = FALSE]
}

# Remake Figure 2 for each number of hidden nodes
brew_colors <- brewer.pal(n = 8, name = "Set1")

named_colors <- c("1" = brew_colors[5],
                  "2" = brew_colors[2],
                  "Truth" = "black", 
                  "True 95% CI" = "black")

named_linetypes <- c("1" = "solid", 
                     "2" = "solid", 
                    "Truth" = "twodash", 
                    "True 95% CI" = "dotted")

# Labels for the legend
labels <- c("1" = "Cluster 1 Estimated Mean",
            "2" = "Cluster 2 Estimated Mean",
            "Truth" = "True Mean", 
            "True 95% CI" = "True 95% Probability Region")

oos_data <- readRDS(
  file.path(
    "1_simulation", 
    "simdata", 
    "simdata_imean_2_iprob_1_iint_5_dataSeed_0.RDS"
  )
)

mns <- c(oos_data$mean_spec, oos_data$pico_mu)
probs <- c(oos_data$prob_spec, 1 - oos_data$prob_spec)

simdat_df <- data.frame(
  x = rep(1:296, 2), 
  y = mns, 
  prob = probs,
  clust = "Truth",
  clust_group = rep(c("1, Truth", "2, Truth"), each = 296)
)

true_1_upper <- oos_data$mean_spec[,1] + 1.96 * 0.2
true_1_lower <- oos_data$mean_spec[,1] - 1.96 * 0.2

true_2_upper <- oos_data$pico_mu + 1.96 * 0.2
true_2_lower <- oos_data$pico_mu - 1.96 * 0.2

simdat_ci_df <- data.frame(
  x = rep(1:296, 4), 
  y = c(true_1_upper, true_1_lower, true_2_upper, true_2_lower), 
  clust = rep(c("1, Truth", "2_Truth"), each = 296 * 2), 
  ci_type = rep(rep(c("upper", "lower"), each = 296), 2), 
  clust_type = "True 95% CI"
)

library(flowtrend)

plist <- vector("list", 5)
plist[[1]] <- plot_1d(oos_data$ybin_list, oos_data$countslist, res_list[[1]], bin = TRUE) + 
  geom_line(
    data = simdat_df,
    mapping = aes(x = x, y = y, group = clust_group, color = clust, linetype = clust), 
    lineend = "round", 
    linewidth = 1
  ) + 
  geom_line(
    data = simdat_ci_df, 
    mapping = aes(
      x = x, 
      y = y, 
      group = interaction(clust, ci_type), 
      linetype = clust_type, 
      color = clust_type
    ), 
    lineend = "round", 
    linewidth = 1
  ) + 
  scale_color_manual(
    name = "", 
    values = named_colors, 
    labels = labels
  ) + 
  scale_linetype_manual(
    name = "", 
    values = named_linetypes, 
    labels = labels
  ) + 
  scale_linewidth_continuous(
    name = "Cluster Probability", 
    range = c(1, 4)
  ) + 
  theme_gray() + 
  theme(
    plot.title = element_text(size = 20, margin = margin(t = 5, b = 10)), 
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 14),
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18, margin = margin(r = 50)), 
    legend.key.size = unit(2, "lines"), 
    legend.margin = margin(0, 0, 0, 0)
  ) + 
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = c(
          "solid", "solid", "dotted", "twodash"
        )
      ), 
      nrow = 2, 
      ncol = 2
    ),
    linetype = "none",
    linewidth = "none",
    fill = "none"
  ) + 
  scale_y_continuous(limits = c(0.25, 2)) + 
  labs(x = "", y = "Data", title = "Estimated Linear Model")

leg <- get_legend(plist[[1]]) %>% 
    as_ggplot()

plist[[1]] <- plist[[1]] + 
    theme(legend.position = "none")

for(i in 2:5) {
  plist[[i]] <- plot_1d(oos_data$ybin_list, oos_data$countslist, res_list[[i]], bin = TRUE) + 
    geom_line(
      data = simdat_df,
      mapping = aes(x = x, y = y, group = clust_group, color = clust, linetype = clust), 
      lineend = "round", 
      linewidth = 1
    ) + 
    geom_line(
      data = simdat_ci_df, 
      mapping = aes(
        x = x, 
        y = y, 
        group = interaction(clust, ci_type), 
        linetype = clust_type, 
        color = clust_type
      ), 
      lineend = "round", 
      linewidth = 1
    ) + 
    scale_color_manual(
      name = "", 
      values = named_colors, 
      labels = labels
    ) + 
    scale_linetype_manual(
      name = "", 
      values = named_linetypes, 
      labels = labels
    ) + 
    scale_linewidth_continuous(
      name = "Cluster Probability", 
      range = c(1, 4)
    ) + 
    theme_gray() + 
    theme(
      plot.title = element_text(size = 20, margin = margin(t = 5, b = 10)), 
      axis.title = element_text(size = 18), 
      axis.text = element_text(size = 14),
      legend.position = "none"
    ) + 
    scale_y_continuous(limits = c(0.25, 2)) + 
    labs(x = "", y = "Data", title = "Estimated Linear Model")
}

plist[[6]] <- text_grob("Time (t)", size = 18, x = 0.5, y = 0)
plist[[7]] <- leg
plist[[8]] <- nullGrob()

pdf(file.path("plots", "SuppFigure11.pdf"), 14.6, 14)
grid.arrange(grobs = plist,
             layout_matrix = matrix(c(1, 2, 3, 
                                      4, 5, 8, 
                                      6, 6, 6, 
                                      7, 7, 7), byrow = TRUE, ncol = 3), 
             widths = c(1, 1, 1), 
             heights = c(0.3, 0.3, 0.0001, 0.1))
graphics.off()
