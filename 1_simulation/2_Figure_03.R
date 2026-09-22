library(RColorBrewer)
library(ggplot2)
library(flowtrend)
library(ggpubr)
library(gridExtra)
library(grid)

simdat <- readRDS(file.path("1_simulation", "simdata", "simdata_imean_2_iprob_1_iint_5_dataSeed_1.RDS"))

ybin_list <- simdat$ybin_list 
countslist <- simdat$countslist
mns <- c(simdat$mean_spec, simdat$pico_mu)
probs <- c(simdat$prob_spec, 1 - simdat$prob_spec)
TT <- length(ybin_list)
numclust <- 2

simdat_df <- data.frame(
  x = rep(1:TT, numclust), 
  y = mns, 
  prob = probs,
  clust = "Truth",
  clust_group = rep(c("1, Truth", "2, Truth"), each = TT)
)

true_1_upper <- simdat$mean_spec[,1] + 1.96 * 0.2
true_1_lower <- simdat$mean_spec[,1] - 1.96 * 0.2

true_2_upper <- simdat$pico_mu + 1.96 * 0.2
true_2_lower <- simdat$pico_mu - 1.96 * 0.2

simdat_ci_df <- data.frame(
  x = rep(1:296, 4), 
  y = c(true_1_upper, true_1_lower, true_2_upper, true_2_lower), 
  clust = rep(c("1, Truth", "2_Truth"), each = 296 * 2), 
  ci_type = rep(rep(c("upper", "lower"), each = 296), 2), 
  clust_type = "True 95% CI"
)

brew_colors <- brewer.pal(n = 8, name = "Set1")

named_colors <- c(
  "1" = brew_colors[5], "2" = brew_colors[2],
  "Truth" = "black", "True 95% CI" = "black"
)

named_linetypes <- c(
  "1" = "solid", "2" = "solid", 
  "Truth" = "twodash", "True 95% CI" = "dotted"
)

labels <- c(
  "1" = "Cluster 1 Estimated Mean", "2" = "Cluster 2 Estimated Mean",
  "Truth" = "True Mean", "True 95% CI" = "True 95% Probability Region"
)

lin_fit <- readRDS(file.path("1_simulation", "results", "2-1-5", "2-1-5-NA-NA_fit_model.RDS"))
nl_fit <- readRDS(file.path("1_simulation", "results", "2-1-5", "2-1-5-70-1_fit_model.RDS"))

lin_plot <- plot_1d(ybin_list, countslist, lin_fit, bin = TRUE) + 
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

leg <- get_legend(lin_plot) %>% 
    as_ggplot()

lin_plot <- lin_plot + 
    theme(legend.position = "none")

# Exact same plot as above but for nonlinear model and no legend
nl_plot <- plot_1d(ybin_list, countslist, nl_fit, bin = TRUE) + 
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
    legend.position = "none",
  ) + 
  scale_y_continuous(limits = c(0.25, 2)) + 
  labs(x = "", y = "Data", title = "Estimated Nonlinear Model")

x_lab <- text_grob("Time (t)", size = 18, x = 0.53, y = 0.9)

pdf(file.path("plots", "Figure03.pdf"), 14.6, 8.5)
grid.arrange(
  lin_plot, nl_plot, x_lab, leg, nullGrob(), 
  layout_matrix = matrix(
    c(
      1, 1, 2, 2,
      3, 3, 3, 3,
      5, 4, 4, 5
    ), 
    byrow = TRUE, 
    nrow = 3
  ), 
  widths = c(0.27, 0.73, 0.95, 0.05), heights = c(1, 0.01, 0.2)
)
graphics.off()