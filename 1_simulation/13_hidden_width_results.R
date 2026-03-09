library(flowmix)
library(flowtrend)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(grid)

cvres_folders <- paste0("2-1-5-", c(35, 70, 105, 140, 175), "-1")
cvres_files <- paste0(cvres_folders, "_summary.RDS")

oos_data <- file.path(
  "1_simulation", 
  "simdata", 
  "simdata_imean_2_iprob_1_iint_5_dataSeed_0.RDS"
) |> 
  readRDS()

# Create oracle model (the model which generates the data)
oracle <- list(
  mn = array(c(oos_data$mean_spec, oos_data$pico_mu), dim = c(296, 1, 2)), 
  prob = array(c(oos_data$prob_spec, 1 - oos_data$prob_spec), dim = c(296, 2)), 
  sigma = array(oos_data$clust_sig^2, dim = c(2, 1, 1))
)

# Get oracle NLL
oracle_oos_nll <- objective_newdat(oracle, oos_data$ybin_list, oos_data$countslist) 

oos_nlls <- numeric(5)
reslist <- list()
for(i in 1:5) {
  # Load CV results
  cv <- readRDS(file.path("1_simulation", "results", cvres_folders[i], cvres_files[i]))

  # Fix label switching between cluster 1 and cluster 2
  clust_1_pos <- cv$bestres$mn[,1,] |>
      apply(2, min) |>
      which.min()

  if(clust_1_pos == 2) {
      cv$bestres$alpha <- cv$bestres$alpha[2:1,1:10]
      cv$bestres$beta <- cv$bestres$beta[2:1]
      cv$bestres$mn <- cv$bestres$mn[,,2:1, drop = FALSE]
      cv$bestres$prob <- cv$bestres$prob[,2:1, drop = FALSE]
      cv$bestres$sigma <- cv$bestres$sigma[2:1, 1, 1, drop = FALSE]
  }

  # Calculate out-of-sample negative log-likelihoods (NLLs) in a data frame
  oos_nlls[i] <- objective_newdat(cv$bestres, oos_data$ybin_list, oos_data$countslist)
  reslist[[i]] <- cv
}

# Get model NLLs, less oracle NLL
nll_oos_oracle_diff <- oos_nlls - oracle_oos_nll

######################3

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

true_1_upper <- oos_data$mean_spec[,1] + 1.96 * 0.2
true_1_lower <- oos_data$mean_spec[,1] - 1.96 * 0.2

true_2_upper <- oos_data$pico_mu + 1.96 * 0.2
true_2_lower <- oos_data$pico_mu - 1.96 * 0.2

ci_df <- data.frame(x = rep(1:296, 4), 
                    y = c(true_1_upper, true_1_lower, true_2_upper, true_2_lower), 
                    Cluster = rep(c("1, Truth", "2_Truth"), each = 296 * 2), 
                    type = rep(rep(c("upper", "lower"), each = 296), 2), 
                    clust = "True 95% CI")

p11 <- plot_sim_means(oos_data, reslist[[1]], plot_model = TRUE) + 
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

leg <- get_legend(p11) %>% 
    as_ggplot()

p11 <- p11 + 
    theme(legend.position = "none")

plot_list <- list(p11)
nh <- seq(35, by = 35, length.out = 5)
for(i in 2:5) {
  cur_nh <- nh[i]

  plot_list[[i]] <- plot_sim_means(oos_data, reslist[[i]], plot_model = TRUE) + 
    geom_line(data = ci_df, 
             mapping = aes(x = x, y = y,
                           group = interaction(Cluster, type), 
                           linetype = clust, 
                           color = clust), 
             lineend = "round", 
             linewidth = 1)  + 
    scale_y_continuous(limits = c(0.25, 2)) + 
    theme_gray() + 
    theme(legend.position = "none", 
          plot.title = element_text(size = 20, margin = margin(t = 5, b = 10)),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14)) + 
    labs(x = "", y = "", title = paste0(cur_nh, " Hidden Nodes"))
}

plot_list[[6]] <- text_grob("Time (t)", size = 18, x = 0.5, y = 0)
plot_list[[7]] <- leg
plot_list[[8]] <- nullGrob()

pdf(file.path("plots", "Figure02_widths_comp.pdf"), 14.6, 14)
grid.arrange(grobs = plot_list,
             layout_matrix = matrix(c(1, 2, 3, 
                                      4, 5, 8, 
                                      6, 6, 6, 
                                      7, 7, 7), byrow = TRUE, ncol = 3), 
             widths = c(1, 1, 1), 
             heights = c(0.3, 0.3, 0.0001, 0.1))
graphics.off()
