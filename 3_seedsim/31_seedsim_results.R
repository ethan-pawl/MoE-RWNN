library(flowmix)
library(flowtrend)
library(RColorBrewer)
library(ggpubr)
library(ggtext)
library(grid)
library(gridExtra)

cv_sums <- list()
reslist <- list()
for(i in 1:5) {
  cv_sums[[i]] <- readRDS(
    file.path(
      "3_seedsim", 
      "results", 
      i, 
      paste0(i, "_summary.RDS")
    )
  )

  reslist[[i]] <- cv_sums[[i]]$bestres 
  
}

oos_data <- readRDS(
    file.path(
      "1_simulation", 
      "simdata", 
      "simdata_imean_2_iprob_1_iint_5_dataSeed_0.RDS"
    )
  )

plot_sim_means <- function(simdat, cv, plot_band = TRUE, plot_model = TRUE) {
    if(plot_model) {
        flowtrend::plot_1d(simdat$ybin_list, simdat$countslist, cv$bestres, 
                           bin = TRUE, plot_band = plot_band) + 
            geom_line(data = data.frame(x = rep(1:296, 2), 
                                        y = c(simdat$mean_spec, simdat$pico_mu), 
                                        prob = c(simdat$prob_spec, 1 - simdat$prob_spec),
                                        clust = "Truth",
                                        clust_group = rep(c("1, Truth", "2, Truth"), 
                                                    each = 296)), 
                      mapping = aes(x = x, y = y, group = clust_group, 
                                    color = clust, linetype = clust), 
                      lineend = "round", linewidth = 1) + 
            scale_color_manual(name = "", 
                               values = named_colors, 
                               labels = labels) + 
            scale_linetype_manual(name = "", 
                                  values = named_linetypes, 
                                  labels = labels) + 
            scale_linewidth_continuous(name = "Cluster Probability", 
                                       range = c(1, 4))
    } else {
        flowtrend::plot_1d(simdat$ybin_list, simdat$countslist, 
                           bin = TRUE) + 
            geom_line(data = data.frame(x = rep(1:296, 2), 
                                        y = c(simdat$mean_spec, simdat$pico_mu), 
                                        prob = c(simdat$prob_spec, 1 - simdat$prob_spec),
                                        clust = rep(c("1, Truth", "2, Truth"), 
                                                    each = 296)), 
                      mapping = aes(x = x, y = y, color = clust, linewidth = prob), lineend = "round") + 
            scale_color_manual(name = "Cluster", 
                               values = c("1, Truth" = "purple", 
                                          "2, Truth" = "#E69F00")) + 
            scale_linewidth_continuous(name = "Cluster Probability", 
                                       range = c(0.25, 4)) + 
            labs(title = get_title(cv, 1))
    } 
}

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

true_1_upper <- oos_data$mean_spec[,1] + 1.96 * 0.2
true_1_lower <- oos_data$mean_spec[,1] - 1.96 * 0.2

true_2_upper <- oos_data$pico_mu + 1.96 * 0.2
true_2_lower <- oos_data$pico_mu - 1.96 * 0.2

ci_df <- data.frame(x = rep(1:296, 4), 
                    y = c(true_1_upper, true_1_lower, true_2_upper, true_2_lower), 
                    Cluster = rep(c("1, Truth", "2_Truth"), each = 296 * 2), 
                    type = rep(rep(c("upper", "lower"), each = 296), 2), 
                    clust = "True 95% CI")

p11 <- plot_sim_means(oos_data, cv_sums[[1]], plot_model = TRUE) + 
    geom_line(data = ci_df, 
             mapping = aes(x = x, y = y,
                           group = interaction(Cluster, type), 
                           linetype = clust, 
                           color = clust), 
             lineend = "round", 
             linewidth = 1) + 
    theme_gray() + 
    theme(plot.title = element_text(size = 20, margin = margin(t = 5, b = 10)), 
          axis.title = element_text(size = 18), 
          axis.text = element_text(size = 14),
          legend.position = "right",
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18, margin = margin(r = 50)), 
          legend.key.size = unit(2, "lines"), 
          legend.margin = margin(0, 0, 0, 0)) + 
    guides(color = guide_legend(override.aes = list(linetype = c("solid",  
                                                                 "solid", 
                                                                 "dotted", 
                                                                 "twodash")), 
                                nrow = 2, ncol = 2),
           linetype = "none",
           linewidth = "none",
           fill = "none") + 
    scale_y_continuous(limits = c(0.25, 2)) + 
    labs(x = "", y = "Data", title = "35 Hidden Nodes")

leg <- get_legend(p11) |> 
    as_ggplot()

p11 <- p11 + 
    theme(legend.position = "none")

plot_list <- list(p11)

for(i in 2:5) {
  plot_list[[i]] <- plot_sim_means(oos_data, cv_sums[[i]], plot_model = TRUE) + 
      geom_line(data = ci_df, 
               mapping = aes(x = x, y = y,
                             group = interaction(Cluster, type), 
                             linetype = clust, 
                             color = clust), 
               lineend = "round", 
               linewidth = 1) + 
      theme_gray() + 
      theme(plot.title = element_text(size = 20, margin = margin(t = 5, b = 10)), 
            axis.title = element_text(size = 18), 
            axis.text = element_text(size = 14),
            legend.position = "none",
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 18, margin = margin(r = 50)), 
            legend.key.size = unit(2, "lines"), 
            legend.margin = margin(0, 0, 0, 0)) + 
      guides(color = guide_legend(override.aes = list(linetype = c("solid",  
                                                                   "solid", 
                                                                   "dotted", 
                                                                   "twodash")), 
                                  nrow = 2, ncol = 2),
             linetype = "none",
             linewidth = "none",
             fill = "none") + 
      scale_y_continuous(limits = c(0.25, 2)) + 
      labs(x = "", y = "Data", title = paste0("Seed ", i))
}

plot_list[[6]] <- text_grob("Time (t)", size = 18, x = 0.5, y = 0)
plot_list[[7]] <- leg
plot_list[[8]] <- nullGrob()

pdf(file.path("plots", "Figure02_seeds_comp.pdf"), 14.6, 14)
grid.arrange(grobs = plot_list,
             layout_matrix = matrix(c(1, 2, 3, 
                                      4, 5, 8, 
                                      6, 6, 6, 
                                      7, 7, 7), byrow = TRUE, ncol = 3), 
             widths = c(1, 1, 1), 
             heights = c(0.3, 0.3, 0.0001, 0.1))
graphics.off()

# First model has switched labels, others are fine

# Compute average of 5 models and standard errors

means <- array(NA, c(296, 5, 2))
means[,1,] <- cv_sums[[1]]$bestres$mn[,,2:1]
for(i in 2:5) {
  means[,i,] <- cv_sums[[i]]$bestres$mn[,1,]
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

pdf(file.path("plots", "Figure02_ribbon.pdf"), 14.6, 8.5)
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
    low = scales::alpha("grey92", 0.1), 
    high = scales::alpha("black", 0.8)
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
  #  + 
  # scale_color_manual(
  #   name = "", 
  #   values = named_colors, 
  #   labels = labels
  # ) + 
  # scale_linetype_manual(
  #   name = "", 
  #   values = named_linetypes, 
  #   labels = labels
  # )
named_colors <- c("1" = brew_colors[5],
                  "2" = brew_colors[2],
                  "Truth" = "black", 
                  "True 95% CI" = "black")