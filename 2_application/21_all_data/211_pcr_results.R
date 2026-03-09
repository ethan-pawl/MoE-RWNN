library(flowmix)
library(parallel)
library(gtools)
library(dplyr)
library(tidyr)
library(tibble)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(maps)
library(RColorBrewer)

nrep <- 10

plots_dir <- "plots"
if(!dir.exists(plots_dir)) dir.create(plots_dir)

# Load all data
datobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper.RDS"))
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()

nl_cv_file <- file.path(
  "2_application", 
  "21_all_data", 
  "results", 
  "nl_9_70_4", 
  "nl_9_70_4_pcr_summary.RDS"
)

nl_cv <- readRDS(nl_cv_file)
nl_best <- nl_cv$bestres

load(file.path("2_application", "21_all_data", "linear_fit.Rdata"))
linear_best <- res 

###################################

load_X <- function(readdir, n_PCs = "NA", n_h = "NA", seed = "NA", ofold = "NA", ifold = "NA") {

  file_name <- paste("X_pc", n_PCs, "nh", n_h, "seed", seed, "ofold", ofold, "ifold", ifold, sep = "_")
  file_name <- paste0(file_name, ".RDS")
  file_name <- file.path(readdir, file_name)

  X <- readRDS(file_name)
}

X_dir <- file.path(
  "data", 
  "X_variations"
)

X_pc <- load_X(X_dir, n_PCs = 9)
X_nl <- load_X(X_dir, n_PCs = 9, n_h = 70, seed = 4)

n.h <- ncol(X_nl)

##########################

# Figure 4
plot_list <- flowtrend::plot_3d(
  ylist, linear_best, 33, countslist, return_list_of_plots = TRUE, 
  labels = c(
    "Syn", "Other4", "Bead", "Other3", "Pico1", "Other6", "Pico2", "Other2", "Other1", "Pro"
  ), 
  mn_colours = c(
    "red", "gray", "gray", "gray", "red", "gray", "red", "gray", "gray", "red"
  )
)
p1 <- plot_list[[1]]

plot_list_nl <- flowtrend::plot_3d(
  ylist, nl_best, 33, countslist, return_list_of_plots = TRUE, 
  labels = c(
    "Pico2", "Pico1", "Other5", "Bead", "Other1", "Syn", "Other6", "Other2", "Pro", "Other4"
  ), 
  mn_colours = c(
    "red", "red", "gray", "gray", "gray", "red", "gray", "gray", "red", "gray"
  )
)
p2 <- plot_list_nl[[1]]

# Store indices for each taxon and model

clust_mat <- matrix(
  c(
    10, 9, 
    1, 6, 
    5, 2, 
    7, 1
  ), 
  ncol = 2, 
  byrow = TRUE
)

rownames(clust_mat) <- c("Pro", "Syn", "Pico1", "Pico2")
colnames(clust_mat) <- c("Linear", "Nonlinear")

p1 <- p1 + 
  labs(title = "Estimated Linear Model", y = "Chlorophyll", x = "") + 
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) + 
  scale_x_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) + 
  scale_y_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) 

p2 <- p2 + 
  labs(title = "Estimated Nonlinear Model", y = "", x = "") + 
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) + 
  scale_x_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) + 
  scale_y_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) 

x_lab <- text_grob("Diameter", size = 12, x = 0.515, y = 1)
overall_title <- text_grob("Estimated Models\n2017-07-02 08:00-09:00", 
                           size = 18)

pdf(file.path(plots_dir, "Figure04.pdf"), 10.5, 6)
grid.arrange(p1, p2, x_lab, overall_title, layout_matrix = matrix(c(4, 4, 
                                                                    1, 2, 
                                                                    3, 3), 
                                                                    byrow = TRUE, 
                                                                    nrow = 3), 
             heights = c(0.175, 1, 0.05), widths = c(1, 1))
graphics.off()

#######################

# Figure 5 and Appendix Figures 12-14
lin_probs <- linear_best$prob[,clust_mat[,"Linear"]]
nl_probs <- nl_best$prob[,clust_mat[,"Nonlinear"]]

dt <- as.POSIXct(time, format = "%Y-%m-%dT%T")

probs <- data.frame(Prediction = c(lin_probs, nl_probs), 
                    Cluster = rep(c("Pro", "Syn", "Pico1", "Pico2", 
                                  "Pro", "Syn", "Pico1", "Pico2"), each = 296), 
                    Model = rep(c("Linear", "Nonlinear"), each = 296 * 4), 
                    Time = rep(dt, 8), 
                    Response = "Relative Abundance", 
                    Response_Type = "Cluster Probability")

probs$Cluster <- factor(probs$Cluster, levels = c("Pro", "Syn", "Pico1", "Pico2"))

# 296 x 3 x 4 
lin_mn <- c(linear_best$mn[,,clust_mat[,"Linear"]])
nl_mn <- c(nl_best$mn[,,clust_mat[,"Nonlinear"]])

mns <- data.frame(Prediction = c(lin_mn, nl_mn), 
                  Cluster = rep(rep(c("Pro", 
                                      "Syn", 
                                      "Pico1", 
                                      "Pico2"), each = 296 * 3), 
                                2),
                  Time = rep(dt, 24), 
                  Model = rep(c("Linear", "Nonlinear"), each = 296 * 3 * 4),
                  Response = rep(rep(c("Diameter", 
                                       "Chlorophyll", 
                                       "Phycoerythrin"), each = 296), 
                                 8), 
                  Response_Type = "Cluster Mean Components")

mns$Cluster <- factor(mns$Cluster, levels = c("Pro", "Syn", "Pico1", "Pico2"))

resp_df <- rbind(probs, mns)

cols <- brewer.pal(12, "Set3")
colors <- c(cols[4], 
            cols[5], 
            cols[6], 
            cols[10])

names(colors) <- c("Diameter", 
                   "Chlorophyll", 
                   "Phycoerythrin",
                   "Relative Abundance")

# Figure 5
pro <- ggplot(filter(resp_df, Cluster  == "Pro")) + 
  geom_line(aes(Time, Prediction, linetype = Model, color = Response), 
            linewidth = 1) + 
  geom_vline(xintercept = dt[203], linetype = "dotted") + 
  geom_vline(xintercept = dt[221], linetype = "dotted") + 
  facet_wrap(~ Response_Type, scales = "free_y") + 
  scale_linetype_manual(values = c("Linear" = "dotdash", 
                                   "Nonlinear" = "solid")) + 
  scale_color_manual(name = "Estimated Quantity", values = colors) + 
  scale_x_datetime(date_breaks = "2 days", date_labels = "%b %d") + 
  theme(text = element_text(size = 22), 
        plot.title = element_text(size = 24, hjust = 0.5, margin = margin(b = 20)),
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 35, hjust = 1), 
        legend.title = element_text(size = 18, margin = margin(b = 10), 
                                    hjust = 0.5), 
        legend.text = element_text(size = 16, margin = margin(r = 10)),
        legend.key.size = unit(2, "lines"),
        legend.position = "bottom", 
        legend.box = "horizontal", 
        legend.box.just = "top",
        legend.spacing.x = unit(1.5, "lines")) + 
  labs(title = expression(italic("Prochlorococcus") * " " * "Estimated Means and Probabilities"), x = "") + 
  guides(color = guide_legend(title.position = "top", nrow = 2, order = 2), 
         linetype = guide_legend(title.position = "top", order = 1))

pro_leg <- get_legend(pro) %>% 
  as_ggplot()

pro <- pro + 
  theme(legend.position = "none")

pdf(file.path(plots_dir, "Figure05.pdf"), 17.7, 8.6)
grid.arrange(pro, pro_leg, nullGrob(), 
             layout_matrix = matrix(c(1, 1, 1, 
                                      3, 2, 3), byrow = TRUE, nrow = 2), 
             widths = c(0.28, 1, 0.11), heights = c(1, 0.2))
graphics.off()

# Appendix Figure 12
syn <- ggplot(filter(resp_df, Cluster  == "Syn")) + 
  geom_line(aes(Time, Prediction, linetype = Model, color = Response), 
            linewidth = 1) + 
  geom_vline(xintercept = dt[203], linetype = "dotted") + 
  geom_vline(xintercept = dt[221], linetype = "dotted") + 
  facet_wrap(~ Response_Type, scales = "free_y") + 
  scale_linetype_manual(values = c("Linear" = "dotdash", 
                                   "Nonlinear" = "solid")) + 
  scale_color_manual(name = "Estimated Quantity", values = colors) + 
  scale_x_datetime(date_breaks = "2 days", date_labels = "%b %d") + 
  theme(text = element_text(size = 22), 
        plot.title = element_text(size = 24, hjust = 0.5, margin = margin(b = 20)),
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 35, hjust = 1), 
        legend.title = element_text(size = 18, margin = margin(b = 10), 
                                    hjust = 0.5), 
        legend.text = element_text(size = 16, margin = margin(r = 10)),
        legend.key.size = unit(2, "lines"),
        legend.position = "bottom", 
        legend.box = "horizontal", 
        legend.box.just = "top",
        legend.spacing.x = unit(1.5, "lines")) + 
  labs(title = expression(italic("Synechococcus") * " " * "Estimated Means and Probabilities"), x = "") + 
  guides(color = guide_legend(title.position = "top", nrow = 2, order = 2), 
         linetype = guide_legend(title.position = "top", order = 1))

syn_leg <- get_legend(syn) %>% 
  as_ggplot()

syn <- syn + 
  theme(legend.position = "none")

pdf(file.path(plots_dir, "Figure12.pdf"), 17.7, 8.6)
grid.arrange(syn, syn_leg, nullGrob(), 
             layout_matrix = matrix(c(1, 1, 1, 
                                      3, 2, 3), byrow = TRUE, nrow = 2), 
             widths = c(0.28, 1, 0.11), heights = c(1, 0.2))
graphics.off()

# Appendix Figure 13
pico1 <- ggplot(filter(resp_df, Cluster  == "Pico1")) + 
  geom_line(aes(Time, Prediction, linetype = Model, color = Response), 
            linewidth = 1) + 
  geom_vline(xintercept = dt[203], linetype = "dotted") + 
  geom_vline(xintercept = dt[221], linetype = "dotted") + 
  facet_wrap(~ Response_Type, scales = "free_y") + 
  scale_linetype_manual(values = c("Linear" = "dotdash", 
                                   "Nonlinear" = "solid")) + 
  scale_color_manual(name = "Estimated Quantity", values = colors) + 
  scale_x_datetime(date_breaks = "2 days", date_labels = "%b %d") + 
  theme(text = element_text(size = 22), 
        plot.title = element_text(size = 24, hjust = 0.5, margin = margin(b = 20)),
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 35, hjust = 1), 
        legend.title = element_text(size = 18, margin = margin(b = 10), 
                                    hjust = 0.5), 
        legend.text = element_text(size = 16, margin = margin(r = 10)),
        legend.key.size = unit(2, "lines"),
        legend.position = "bottom", 
        legend.box = "horizontal", 
        legend.box.just = "top",
        legend.spacing.x = unit(1.5, "lines")) + 
  labs(title = expression(italic("PicoEukaryote") * " " * "Cluster 1 Estimated Means and Probabilities"), x = "") + 
  guides(color = guide_legend(title.position = "top", nrow = 2, order = 2), 
         linetype = guide_legend(title.position = "top", order = 1))

pico1_leg <- get_legend(pico1) %>% 
  as_ggplot()

pico1 <- pico1 + 
  theme(legend.position = "none")

pdf(file.path(plots_dir, "Figure13.pdf"), 17.7, 8.6)
grid.arrange(pico1, pico1_leg, nullGrob(), 
             layout_matrix = matrix(c(1, 1, 1, 
                                      3, 2, 3), byrow = TRUE, nrow = 2), 
             widths = c(0.28, 1, 0.11), heights = c(1, 0.2))
graphics.off()

# Appendix Figure 14
pico2 <- ggplot(filter(resp_df, Cluster  == "Pico2")) + 
  geom_line(aes(Time, Prediction, linetype = Model, color = Response), 
            linewidth = 1) + 
  geom_vline(xintercept = dt[203], linetype = "dotted") + 
  geom_vline(xintercept = dt[221], linetype = "dotted") + 
  facet_wrap(~ Response_Type, scales = "free_y") + 
  scale_linetype_manual(values = c("Linear" = "dotdash", 
                                   "Nonlinear" = "solid")) + 
  scale_color_manual(name = "Estimated Quantity", values = colors) + 
  scale_x_datetime(date_breaks = "2 days", date_labels = "%b %d") + 
  theme(text = element_text(size = 22), 
        plot.title = element_text(size = 24, hjust = 0.5, margin = margin(b = 20)),
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 35, hjust = 1), 
        legend.title = element_text(size = 18, margin = margin(b = 10), 
                                    hjust = 0.5), 
        legend.text = element_text(size = 16, margin = margin(r = 10)),
        legend.key.size = unit(2, "lines"),
        legend.position = "bottom", 
        legend.box = "horizontal", 
        legend.box.just = "top",
        legend.spacing.x = unit(1.5, "lines")) + 
  labs(title = expression(italic("PicoEukaryote") * " " * "Cluster 2 Estimated Means and Probabilities"), x = "") + 
  guides(color = guide_legend(title.position = "top", nrow = 2, order = 2), 
         linetype = guide_legend(title.position = "top", order = 1))

pico2_leg <- get_legend(pico2) %>% 
  as_ggplot()

pico2 <- pico2 + 
  theme(legend.position = "none")

pdf(file.path(plots_dir, "Figure14.pdf"), 17.7, 8.6)
grid.arrange(pico2, pico2_leg, nullGrob(), 
             layout_matrix = matrix(c(1, 1, 1, 
                                      3, 2, 3), byrow = TRUE, nrow = 2), 
             widths = c(0.28, 1, 0.11), heights = c(1, 0.2))
graphics.off()

##################################

# Figures 6 and 7, and Appendix Figures 8 & 9

# Create fine grids of X values
PC1_seq <- seq(min(X_pc[,1]), max(X_pc[,1]), length.out = 30)
PC2_seq <- seq(min(X_pc[,2]), max(X_pc[,2]), length.out = 30)

# Start with all PCs fixed at their mean
PC1_ice_X <- tcrossprod(rep(1, 30), colMeans(X_pc))
# vary PC1 
PC1_ice_X[,1] <- PC1_seq

# Start with all PCs fixed at their mean
other_ice_X <- tcrossprod(rep(1, 30 * 3), colMeans(X_pc))

# check the ICEs of PC2, PC3, and PC4 with PC1 fixed at values 
# corresponding to low, medium, and high latitude
PC1_levels <- c(-6, -1, 4)
other_ice_X[,1] <- rep(PC1_levels, each = 30)

PC2_ice_X <- other_ice_X

PC2_ice_X[,2] <- rep(PC2_seq, times = 3)

# Convert PCs to nonlinear representation
set.seed(4)
PC1_ice_X_nl <- make_hidden_nodes(PC1_ice_X, n.h, 0.5)
set.seed(4)
PC2_ice_X_nl <- make_hidden_nodes(PC2_ice_X, n.h, 0.5)

# Make model predictions
PC1_lin_ice <- predict(linear_best, newx = PC1_ice_X, logits = FALSE)
PC1_nl_ice <- predict(nl_best, newx = PC1_ice_X_nl, logits = FALSE)
PC2_lin_ice <- predict(linear_best, newx = PC2_ice_X, logits = FALSE)
PC2_nl_ice <- predict(nl_best, newx = PC2_ice_X_nl, logits = FALSE)

# For each response (cluster probability and each mean component), 
# Re-format the data into a data frame
PC1_lin_probs <- PC1_lin_ice$prob[,clust_mat[,"Linear"]]
PC1_nl_probs <- PC1_nl_ice$prob[,clust_mat[,"Nonlinear"]]
colnames(PC1_lin_probs) <- c("Pro", "Syn", "Pico1", "Pico2")
colnames(PC1_nl_probs) <- c("Pro", "Syn", "Pico1", "Pico2")
PC1_lin_probs <- as.data.frame(PC1_lin_probs)
PC1_nl_probs <- as.data.frame(PC1_nl_probs)
PC1_lin_probs$PC1 <- PC1_nl_probs$PC1 <- PC1_seq
PC1_lin_probs$Response <- PC1_nl_probs$Response <- "Relative Abundance"
PC1_lin_probs$Model <- "Linear" 
PC1_nl_probs$Model <- "Nonlinear"

PC1_lin_diam <- PC1_lin_ice$mn[,1,clust_mat[,"Linear"]]
PC1_nl_diam <- PC1_nl_ice$mn[,1,clust_mat[,"Nonlinear"]]
colnames(PC1_lin_diam) <- c("Pro", "Syn", "Pico1", "Pico2")
colnames(PC1_nl_diam) <- c("Pro", "Syn", "Pico1", "Pico2")
PC1_lin_diam <- as.data.frame(PC1_lin_diam)
PC1_nl_diam <- as.data.frame(PC1_nl_diam)
PC1_lin_diam$PC1 <- PC1_nl_diam$PC1 <- PC1_seq
PC1_lin_diam$Response <- PC1_nl_diam$Response <- "Diameter"
PC1_lin_diam$Model <- "Linear" 
PC1_nl_diam$Model <- "Nonlinear"

PC1_lin_chl <- PC1_lin_ice$mn[,2,clust_mat[,"Linear"]]
PC1_nl_chl <- PC1_nl_ice$mn[,2,clust_mat[,"Nonlinear"]]
colnames(PC1_lin_chl) <- c("Pro", "Syn", "Pico1", "Pico2")
colnames(PC1_nl_chl) <- c("Pro", "Syn", "Pico1", "Pico2")
PC1_lin_chl <- as.data.frame(PC1_lin_chl)
PC1_nl_chl <- as.data.frame(PC1_nl_chl)
PC1_lin_chl$PC1 <- PC1_nl_chl$PC1 <- PC1_seq
PC1_lin_chl$Response <- PC1_nl_chl$Response <- "Chlorophyll"
PC1_lin_chl$Model <- "Linear" 
PC1_nl_chl$Model <- "Nonlinear"

PC1_lin_pe <- PC1_lin_ice$mn[,3,clust_mat[,"Linear"]]
PC1_nl_pe <- PC1_nl_ice$mn[,3,clust_mat[,"Nonlinear"]]
colnames(PC1_lin_pe) <- c("Pro", "Syn", "Pico1", "Pico2")
colnames(PC1_nl_pe) <- c("Pro", "Syn", "Pico1", "Pico2")
PC1_lin_pe <- as.data.frame(PC1_lin_pe)
PC1_nl_pe <- as.data.frame(PC1_nl_pe)
PC1_lin_pe$PC1 <- PC1_nl_pe$PC1 <- PC1_seq
PC1_lin_pe$Response <- PC1_nl_pe$Response <- "Phycoerythrin"
PC1_lin_pe$Model <- "Linear" 
PC1_nl_pe$Model <- "Nonlinear"

# Collect in a single data frame
PC1_pred <- rbind(PC1_lin_probs, PC1_nl_probs, 
                  PC1_lin_diam, PC1_nl_diam, 
                  PC1_lin_chl, PC1_nl_chl, 
                  PC1_lin_pe, PC1_nl_pe)

PC1_pred_melt <- melt(PC1_pred, id.vars = c("PC1", "Model", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

PC1_pred_melt$Response <- factor(PC1_pred_melt$Response, 
  levels = c("Relative Abundance", "Diameter", "Chlorophyll", "Phycoerythrin"))

# Plotting a sign flip of PC1 since it fits intuition better 
# (positive correlation with latitude)

# Appendix Figure 8
PC1_plot <- ggplot(PC1_pred_melt) + 
  geom_line(aes(x = PC1, y = Prediction, linetype = Model, color = Model), linewidth = 1) + 
  scale_linetype_manual(values = c("Linear" = "dotdash", "Nonlinear" = "solid")) + 
  facet_wrap(Population ~ Response, scales = "free") + 
  geom_vline(xintercept = X_pc[203,1], linetype = "dotted") + 
  geom_vline(xintercept = X_pc[221,1], linetype = "dotted") +
  theme(plot.title = element_text(size = 20), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 18),
        legend.key.size = unit(2, "lines"),
        axis.title = element_text(size = 16), 
        text = element_text(size = 14)) + 
  labs(title = "Expected Cluster Means and Probabilities Conditional on Principal Component 1", 
       x = "PC1 (Proxy for Latitude)")

pdf(file.path(plots_dir, "Figure08.pdf"), 13.2, 11.7)
PC1_plot
graphics.off()

# Output a list of these plots
PC1_ice_list <- lapply(c("Pro", "Syn", "Pico1", "Pico2"), function(pop) {
  lapply(c("Relative Abundance", "Diameter", "Chlorophyll", "Phycoerythrin"), function(resp) {
    ggplot(filter(PC1_pred_melt, Population == pop, Response == resp)) + 
  geom_line(aes(x = PC1, y = Prediction, linetype = Model, color = Model), linewidth = 1) + 
  scale_linetype_manual(name = "Model", labels = c("Linear", "Nonlinear"), values = c("Linear" = "dotdash", "Nonlinear" = "solid")) + 
  scale_color_discrete(name = "Model", labels = c("Linear", "Nonlinear")) +
  geom_vline(xintercept = X_pc[203,1], linetype = "dotted") + 
  geom_vline(xintercept = X_pc[221,1], linetype = "dotted") +
  theme(text = element_text(size = 22)) + 
  labs(title = paste(pop, resp))
  })
}) %>% unlist(recursive = FALSE) 

# Figure 6
xaxis_labs <- paste0(c("-8\n23.7",
                       "-4\n32.4", 
                       "0\n34.2", 
                       "4\n39.6", 
                       "8\n>42.1"), "° N")

syn_pe <- PC1_ice_list[[8]] + 
  scale_x_continuous(breaks = seq(-8, 8, by = 4), 
                     labels = xaxis_labs) + 
  labs(title = expression(italic("Synechococcus")), y = "Phycoerythrin", x = "") + 
  theme(plot.title = element_text(size = 20), 
        plot.margin = margin(r = 35, b = 5),
        legend.title = element_text(size = 22), 
        legend.key.size = unit(2, "lines"),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        axis.text = element_text(size = 16))

leg <- get_legend(syn_pe) %>% 
  as_ggplot()

syn_pe <- syn_pe + 
  theme(legend.position = "none")

pico2_prob <- PC1_ice_list[[13]] + 
  scale_x_continuous(breaks = seq(-8, 8, by = 4), 
                     labels = xaxis_labs) + 
  labs(title = expression(italic("PicoEukaryote") * " " * "Cluster 2"), y = "Relative Abundance", x = "") + 
  theme(legend.position = "none", 
        plot.title = element_text(size = 20), 
        plot.margin = margin(r = 35, b = 5),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20), 
        axis.text = element_text(size = 16))

pro_diam <- PC1_ice_list[[2]] + 
  scale_x_continuous(breaks = seq(-8, 8, by = 4), 
                     labels = xaxis_labs) + 
  labs(title = expression(italic("Prochlorococcus")), y = "Diameter", x = "") + 
  theme(legend.position = "none", 
        plot.title = element_text(size = 20), 
        plot.margin = margin(r = 35, b = 5),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20), 
        axis.text = element_text(size = 16))

x_lab <- text_grob("PC1 (Proxy for Latitude)", size = 18, 
                   hjust = 0.375, vjust = -0.5)

pdf(file.path(plots_dir, "Figure06.pdf"), 18.3, 5.6)
grid.arrange(pro_diam, syn_pe, pico2_prob, leg, x_lab, nrow = 2, 
            layout_matrix = matrix(c(1, 2, 3, 4, 
                                     NA, 5, NA, NA), 2, byrow = TRUE), 
            widths = c(1, 1, 1, 0.4), heights = c(1, 0.05))
graphics.off()

# Figure 7
plot_ice <- function(PC_lin_ice, PC_nl_ice, PC_num, PC_seq, return_plot_list = FALSE) {
  PC_lin_probs <- PC_lin_ice$prob[,clust_mat[,"Linear"]]
  PC_nl_probs <- PC_nl_ice$prob[,clust_mat[,"Nonlinear"]]
  colnames(PC_lin_probs) <- c("Pro", "Syn", "Pico1", "Pico2")
  colnames(PC_nl_probs) <- c("Pro", "Syn", "Pico1", "Pico2")
  PC_lin_probs <- as.data.frame(PC_lin_probs)
  PC_nl_probs <- as.data.frame(PC_nl_probs)
  PC_lin_probs$PC1 <- PC_nl_probs$PC1 <- factor(other_ice_X[,1], levels = c(-6, -1, 4), 
                                          labels = c("-6", "-1", "4"))
  PC_lin_probs[,PC_num] <- PC_nl_probs[,PC_num] <- PC_seq
  PC_lin_probs$Response <- PC_nl_probs$Response <- "Relative Abundance"
  PC_lin_probs$Model <- "Linear" 
  PC_nl_probs$Model <- "Nonlinear"

  PC_lin_diam <- PC_lin_ice$mn[,1,clust_mat[,"Linear"]]
  PC_nl_diam <- PC_nl_ice$mn[,1,clust_mat[,"Nonlinear"]]
  colnames(PC_lin_diam) <- c("Pro", "Syn", "Pico1", "Pico2")
  colnames(PC_nl_diam) <- c("Pro", "Syn", "Pico1", "Pico2")
  PC_lin_diam <- as.data.frame(PC_lin_diam)
  PC_nl_diam <- as.data.frame(PC_nl_diam)
  PC_lin_diam$PC1 <- PC_nl_diam$PC1 <- factor(other_ice_X[,1], levels = c(-6, -1, 4), 
                                          labels = c("-6", "-1", "4"))
  PC_lin_diam[,PC_num] <- PC_nl_diam[,PC_num] <- PC_seq
  PC_lin_diam$Response <- PC_nl_diam$Response <- "Diameter"
  PC_lin_diam$Model <- "Linear" 
  PC_nl_diam$Model <- "Nonlinear"

  PC_lin_chl <- PC_lin_ice$mn[,2,clust_mat[,"Linear"]]
  PC_nl_chl <- PC_nl_ice$mn[,2,clust_mat[,"Nonlinear"]]
  colnames(PC_lin_chl) <- c("Pro", "Syn", "Pico1", "Pico2")
  colnames(PC_nl_chl) <- c("Pro", "Syn", "Pico1", "Pico2")
  PC_lin_chl <- as.data.frame(PC_lin_chl)
  PC_nl_chl <- as.data.frame(PC_nl_chl)
  PC_lin_chl$PC1 <- PC_nl_chl$PC1 <- factor(other_ice_X[,1], levels = c(-6, -1, 4), 
                                          labels = c("-6", "-1", "4"))
  PC_lin_chl[,PC_num] <- PC_nl_chl[,PC_num] <- PC_seq
  PC_lin_chl$Response <- PC_nl_chl$Response <- "Chlorophyll"
  PC_lin_chl$Model <- "Linear" 
  PC_nl_chl$Model <- "Nonlinear"

  PC_lin_pe <- PC_lin_ice$mn[,3,clust_mat[,"Linear"]]
  PC_nl_pe <- PC_nl_ice$mn[,3,clust_mat[,"Nonlinear"]]
  colnames(PC_lin_pe) <- c("Pro", "Syn", "Pico1", "Pico2")
  colnames(PC_nl_pe) <- c("Pro", "Syn", "Pico1", "Pico2")
  PC_lin_pe <- as.data.frame(PC_lin_pe)
  PC_nl_pe <- as.data.frame(PC_nl_pe)
  PC_lin_pe$PC1 <- PC_nl_pe$PC1 <- factor(other_ice_X[,1], levels = c(-6, -1, 4), 
                                          labels = c("-6", "-1", "4"))
  PC_lin_pe[,PC_num] <- PC_nl_pe[,PC_num] <- PC_seq
  PC_lin_pe$Response <- PC_nl_pe$Response <- "Phycoerythrin"
  PC_lin_pe$Model <- "Linear" 
  PC_nl_pe$Model <- "Nonlinear"

  PC_pred <- rbind(PC_lin_probs, PC_nl_probs, 
                    PC_lin_diam, PC_nl_diam, 
                    PC_lin_chl, PC_nl_chl, 
                    PC_lin_pe, PC_nl_pe)

  PC_pred_melt <- melt(PC_pred, id.vars = c("PC1", PC_num, "Model", "Response"), 
                                    variable.name = "Population", 
                                    value.name = "Prediction")

  # Labeling PC1 values with latitudes
  PC_pred_melt$PC1 <- factor(PC_pred_melt$PC1, 
                               levels = c(-6, -1, 4), 
                               labels = c("-6 (28.5° N)", "-1 (33.2° N)", "4 (39.6° N)"))
  

  PC_pred_melt$Response <- factor(PC_pred_melt$Response, 
    levels = c("Relative Abundance", "Diameter", "Chlorophyll", "Phycoerythrin"))

  if(return_plot_list) {
    lapply(c("Pro", "Syn", "Pico1", "Pico2"), function(pop) {
      lapply(c("Relative Abundance", "Diameter", "Chlorophyll", "Phycoerythrin"), function(resp) {
        ggplot(filter(PC_pred_melt, Response == resp, Population == pop)) + 
          geom_line(aes(x = .data[[PC_num]], y = Prediction, linetype = Model, color = PC1), linewidth = 1) + 
          scale_color_discrete(name = "PC1") + 
          scale_linetype_manual(name = "Model", labels = c("Linear", "Nonlinear"), values = c("Linear" = "dotdash", "Nonlinear" = "solid")) + 
          labs(title = pop, y = resp) + 
          theme(plot.title = element_text(size = 20), 
                legend.title = element_text(size = 22), 
                legend.text = element_text(size = 18),
                legend.key.size = unit(2, "lines"),
                axis.title = element_text(size = 16), 
                text = element_text(size = 14))
      })
    }) %>% unlist(recursive = FALSE)
  } else {
    ggplot(PC_pred_melt) + 
          geom_line(aes(x = .data[[PC_num]], y = Prediction, linetype = Model, color = PC1), linewidth = 1) + 
          scale_color_discrete(name = "PC1") + 
          scale_linetype_manual(name = "Model", values = c("Linear" = "dotdash", "Nonlinear" = "solid")) + 
          facet_wrap(Population ~ Response, scales = "free") + 
          theme(plot.title = element_text(size = 20), 
                legend.title = element_text(size = 22), 
                legend.text = element_text(size = 18),
                legend.key.size = unit(2, "lines"),
                axis.title = element_text(size = 16), 
                axis.title.y = element_text(size = 18),
                axis.text = element_text(size = 16),
                text = element_text(size = 14)) +
          labs(title = "Expected Cluster Means and Probabilities Conditional on Principal Component 2")
  }
}

# Figure 7 and Appendix Figure 9
PC2_plot <- plot_ice(PC2_lin_ice, PC2_nl_ice, "PC2", PC2_seq)

# Appendix Figure 9
pdf(file.path(plots_dir, "Figure09.pdf"), 12.5, 11.4)
PC2_plot
graphics.off()

# List of plots
PC2_plist <- plot_ice(PC2_lin_ice, PC2_nl_ice, "PC2", PC2_seq, return_plot_list = TRUE)

# Figure 7
pro_prob <- PC2_plist[[1]] + 
  labs(x = "", title = expression(italic("Prochlorococcus"))) + 
  theme(legend.title = element_text(size = 22, margin = margin(b = 10), 
                                    hjust = 0), 
        legend.text = element_text(size = 18, margin = margin(l = 5, r = 10)),
        legend.key.size = unit(2, "lines"),
        legend.box = "vertical", 
        legend.box.just = "left",
        legend.spacing.x = unit(1.5, "lines")) + 
  guides(linetype = guide_legend(title.position = "top", ncol = 1), 
         color = guide_legend(title.position = "top", ncol = 1))

PC2_leg <- ggpubr::get_legend(pro_prob) %>% 
  ggpubr::as_ggplot()

pro_prob <- pro_prob + 
  theme(legend.position = "none")

syn_pe <- PC2_plist[[8]] +
  theme(legend.position = "none") + 
  labs(x = "", title = expression(italic("Synechococcus")))

pro_chl <- PC2_plist[[3]] + 
  theme(legend.position = "none") + 
  labs(x = "", title = expression(italic("Prochlorococcus")))

PC2_xlab <- ggpubr::text_grob("PC2", 
                              hjust = 0.425, vjust = 0, size = 16)

# Map (bottom panel of Figure 7)
world_map <- map_data("world")

# Define limits
lon_min <- -161
lon_max <- -100
lat_min <- 19
lat_max <- 50

world_map <- world_map %>%
  filter(long >= lon_min, long <= 0,
         lat  >= lat_min, lat  <= 70)

POIs <- c("Gradients 2 Cruise", "PC1 = 4", 
          "PC1 = -1", "PC1 = -6",
          "NPTZ Northern Limit", "NPTZ Southern Limit")

lat_df <- data.frame(
  Latitude = c(39.6, 33.2, 28.5, 34, 31.9), 
  "Points of Interest" = factor(POIs[2:6], levels = POIs), 
  check.names = FALSE
)

cruisepath_df <- data.frame(
  x = c(-158, -157), 
  y = c(33.2, 42.0), 
  xend = c(-158, -157), 
  yend = c(42.4, 22.1), 
  "Cruise Trajectory" = c("Northward", "Southward"), 
  "Points of Interest" = factor("Gradients 2 Cruise", levels = POIs),
  check.names = FALSE
)

POI_linewidths <- c(1, 0.5, 1, 0.5, 1, 0.5)
POI_colors <- c("black", "#619CFF", "#00BA38", "#F8766D", "black", "black")
POI_linetypes <- c("solid", "solid", "solid", "solid", "dotted", "dotted")

names(POI_linewidths) <- names(POI_colors) <- names(POI_linetypes) <- POIs

cruise_map <- ggplot(world_map, aes(x = long, y = lat, group = group)) + 
  geom_polygon(fill = "antiquewhite", color = "gray40") +
  coord_fixed(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
  geom_hline(data = lat_df, 
             aes(yintercept = Latitude, color = `Points of Interest`, 
                 linetype = `Points of Interest`), linewidth = 1) + 
  geom_segment(data = cruisepath_df, 
               aes(x = x, y = y, xend = xend, yend = yend, group = `Cruise Trajectory`, 
                   color = `Points of Interest`), linewidth = 1, 
               arrow = arrow(length = unit(0.25, "cm"))) + 
  geom_point(aes(x = -158, y = 33.2), size = 3) + 
  geom_point(aes(x = -157, y = 42.15), size = 3) + 
  scale_color_manual(values = POI_colors, drop = FALSE) + 
  scale_linetype_manual(values = POI_linetypes, drop = FALSE) + 
  labs(title = "Map of Gradients 2 Cruise", 
       x = "Longitude", 
       y = "Latitude") + 
  theme(plot.title = element_text(size = 20), 
        legend.title = element_text(size = 22, margin = margin(b = 10)), 
        legend.text = element_text(size = 18, margin = margin(l = 5, r = 10)),
        legend.key.size = unit(2, "lines"),
        axis.title = element_text(size = 16), 
        text = element_text(size = 14), 
        legend.margin = margin(l = 55))

cruise_map_leg <- get_legend(cruise_map) %>% 
  as_ggplot()

cruise_map <- cruise_map + 
  theme(legend.position = "none")

pdf(file.path(plots_dir, "Figure07.pdf"), 17.5, 10.6)
grid.arrange(pro_prob, syn_pe, pro_chl, PC2_leg, PC2_xlab, 
             nullGrob(), cruise_map, cruise_map_leg,
             layout_matrix = matrix(c(1, 2, 3, 4,
                                      5, 5, 5, 6,
                                      7, 7, 8, 6), nrow = 3, byrow = TRUE), 
             widths = c(1, 1, 1, 0.5), heights = c(1, 0.05, 1))
graphics.off()
