library(flowmix)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(parallel)
library(tibble)
library(reshape2)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)

# Collect results
linear_cv_file <- file.path("2_application", 
                            "21_all_data", 
                            "results", 
                            "linear", 
                            "pcr_summary_linear.RDS")

nl_cv_file <- file.path("2_application", 
                        "21_all_data", 
                        "results", 
                        "nl", 
                        "pcr_summary_nl.RDS")

linear_cv <- readRDS(linear_cv_file)
nl_cv <- readRDS(nl_cv_file)

linear_best <- linear_cv$bestres
nl_best <- nl_cv$bestres

# Load all data
datobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper.RDS"))
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()

###################################

load(file.path("data", "X_pc.Rdata"))
X_pc <- X

load(file.path("data", "X_nl.Rdata"))
X_nl <- X

##########################

# Figure 4
plot_list <- flowtrend::plot_3d(ylist, linear_best, 33, countslist, 
                                return_list_of_plots = TRUE, labels = c("Syn", 
                                                                        "Other1", 
                                                                        "Bead", 
                                                                        "Other2", 
                                                                        "Pico1", 
                                                                        "Other3", 
                                                                        "Pico2", 
                                                                        "Other4", 
                                                                        "Other5", 
                                                                        "Pro"))
p1 <- plot_list[[1]]

plot_list_nl <- flowtrend::plot_3d(ylist, nl_best, 33, countslist, 
                                return_list_of_plots = TRUE, labels = c("Other6", 
                                                                        "Other3", 
                                                                        "Syn", 
                                                                        "Other4", 
                                                                        "Pico2", 
                                                                        "Other1", 
                                                                        "Pico1", 
                                                                        "Pro", 
                                                                        "Other2", 
                                                                        "Other5"))
p2 <- plot_list_nl[[1]]

p1 <- p1 + 
  labs(title = "Linear Model Fit", y = "Chlorophyll", x = "") + 
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) + 
  scale_x_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) + 
  scale_y_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) 

p2 <- p2 + 
  labs(title = "Nonlinear Model Fit", y = "", x = "") + 
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) + 
  scale_x_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) + 
  scale_y_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) 

x_lab <- text_grob("Log Diameter", size = 12, x = 0.515, y = 1)
overall_title <- text_grob("Model Fit on July 2nd, 2017, 8 AM", size = 16)

pdf(file.path("plots", "clusters_33_onerow.pdf"), 10.5, 6)
grid.arrange(p1, p2, x_lab, overall_title, layout_matrix = matrix(c(4, 4, 
                                                                    1, 2, 
                                                                    3, 3), 
                                                                    byrow = TRUE, 
                                                                    nrow = 3), 
             heights = c(0.1, 1, 0.05), widths = c(1, 1, 0.05))
graphics.off()

# TODO: Check what edits I made to make the paper plot
trace(flowtrend::plot_3d, edit = TRUE)

# Figure 5

lin_probs <- linear_best$prob[,c(10, 1, 5, 7)]
nl_probs <- nl_best$prob[,c(8, 3, 7, 5)]

time_dt <- as.POSIXct(time, format = "%Y-%m-%d-")
dt <- as.POSIXct(time, format = "%Y-%m-%dT%T")

probs <- data.frame(Relative_Abundance = c(lin_probs, nl_probs), 
                    Cluster = rep(c("Pro", "Syn", "Pico1", "Pico2", 
                                  "Pro", "Syn", "Pico1", "Pico2"), each = 296), 
                    Model = rep(c("Linear", "Nonlinear"), each = 296 * 4), 
                    Time = rep(dt, 8))
probs$Cluster <- factor(probs$Cluster, levels = c("Pro", "Syn", "Pico1", "Pico2"))

pdf(file.path("plots", "prob_time.pdf"), 12.8, 5.3)
ggplot(probs) + 
  geom_line(aes(Time, Relative_Abundance, color = Model, linetype = Model), linewidth = 0.75) + 
  facet_wrap(~ Cluster, scales = "free") + 
  scale_linetype_manual(values = c("Linear" = "longdash", "Nonlinear" = "solid")) + 
  # scale_y_continuous(limits = c(0, 1)) + 
  scale_x_datetime(date_breaks = "2 days", date_labels = "%b %d") + 
  labs(y = "Relative Abundance") + 
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 30, hjust = 1))
graphics.off()

###### PCA Line Plots

# PC1

# Figure 10
lat_scaled <- (lat - min(lat)) / (max(lat) - min(lat)) * (10) - 1.5

cov_pc_df <- cbind(X, X_pc)
cov_pc_df <- as.data.frame(cov_pc_df)
cov_pc_df$time <- lubridate::as_datetime(rownames(cov_pc_df))
cov_pc_df$lat <- lat_scaled

colnames(cov_pc_df)[c(26, 36, 37)] <- c("sdns", "nitrate", "phosphate")

cov_pc_df_long <- melt(cov_pc_df, id.vars = "time", variable.name = "Variable", value.name = "Value")
head(cov_pc_df_long)

pc1_covs <- c("sss", "sst", "PP", "Si", "NO3", "CHL", "PHYC", "PO4", "O2", "PC1")
pc1_covs_lat <- c(pc1_covs, "lat")

pdf(file.path("plots", "PC1_line.pdf"), 14, 7.9)
ggplot(cov_pc_df_long[cov_pc_df_long$Variable %in% pc1_covs_lat,]) + 
  geom_path(data = cov_pc_df_long[cov_pc_df_long$Variable %in% pc1_covs,], 
          aes(time, Value, group = Variable, color = Variable, 
              linewidth = Variable)) + 
  scale_x_datetime(date_breaks = "2 days") + 
  guides(x = guide_axis(angle = 15)) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(9, "Set1"), "black")) + 
  scale_linewidth_manual(values = c(rep(1, 9), 2)) + 
  geom_path(data = cov_pc_df_long[cov_pc_df_long$Variable == "lat",], 
            aes(time, Value), color = "black", linetype = "dashed", linewidth = 1) + 
  annotate("text", x = cov_pc_df$time[130], y = 8, label = "latitude (for comparison;\nnot used as a covariate)", 
            fontface = "bold", size = 5) + 
  geom_vline(xintercept = cov_pc_df$time[203], linetype = "dotdash", linewidth = 1) + 
  geom_vline(xintercept = cov_pc_df$time[221], linetype = "dotdash", linewidth = 1) + 
  annotate("text", x = cov_pc_df$time[265], y = -3.5, 
          label = "latitudes 34°N-31.9°N;\nconditions change\nabruptly", 
          fontface = "bold", size = 5) + 
  annotate("text", x = cov_pc_df$time[140], y = -5.5, 
            label = "1st Principal Component:\n\"negative latitude,\"\nlearned from\n environmental conditions", 
            fontface = "bold", size = 5) + 
  theme(text = element_text(size = 22), legend.key.size = unit(1, "cm"))
graphics.off()

###### PC2

# Figure 11
pc2_covs <- c("Fe", "sla", "wind_stress", "wind_speed", "vgos", "vgosa", 
              "sdns", "northward_wind", "nitrate", "phosphate", "PC2")

pdf(file.path("plots", "PC2_line.pdf"), 14, 7.9)
ggplot(cov_pc_df_long[cov_pc_df_long$Variable %in% pc2_covs,]) + 
  geom_path(data = cov_pc_df_long[cov_pc_df_long$Variable %in% pc2_covs,], 
          aes(time, Value, group = Variable, color = Variable, 
              linewidth = Variable)) + 
  scale_x_datetime(date_breaks = "2 days") + 
  guides(x = guide_axis(angle = 15)) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(10, "Set3"), "black")) + 
  scale_linewidth_manual(values = c(rep(1, 10), 2)) + 
  geom_vline(xintercept = cov_pc_df$time[203], linetype = "dotdash", linewidth = 1) + 
  geom_vline(xintercept = cov_pc_df$time[221], linetype = "dotdash", linewidth = 1) + 
  theme(text = element_text(size = 22), legend.key.size = unit(1, "cm"))
graphics.off()

###### PC3

# Figure 12
pc3_covs <- c("p1", "p2", "ugos", "ugosa", "PC3")

pdf(file.path("plots", "PC3_line.pdf"), 14, 7.9)
ggplot(cov_pc_df_long[cov_pc_df_long$Variable %in% pc3_covs,]) + 
  geom_path(data = cov_pc_df_long[cov_pc_df_long$Variable %in% pc3_covs,], 
          aes(time, Value, group = Variable, color = Variable, 
              linewidth = Variable)) + 
  scale_x_datetime(date_breaks = "2 days") + 
  guides(x = guide_axis(angle = 15)) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(4, "Set3"), "black")) + 
  scale_linewidth_manual(values = c(rep(1, 4), 2)) + 
  geom_vline(xintercept = cov_pc_df$time[203], linetype = "dotdash", linewidth = 1) + 
  geom_vline(xintercept = cov_pc_df$time[221], linetype = "dotdash", linewidth = 1) + 
  theme(text = element_text(size = 22), legend.key.size = unit(1, "cm"))
graphics.off()

###### PC4

# Figure 13
pc4_covs <- c("p3", "p4", "par", "PC4")

pdf(file.path("plots", "PC4_line.pdf"), 14, 7.9)
ggplot(cov_pc_df_long[cov_pc_df_long$Variable %in% pc4_covs,]) + 
  geom_path(data = cov_pc_df_long[cov_pc_df_long$Variable %in% pc4_covs,], 
          aes(time, Value, group = Variable, color = Variable, 
              linewidth = Variable)) + 
  scale_x_datetime(date_breaks = "2 days") + 
  guides(x = guide_axis(angle = 15)) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(3, "Set3"), "black")) + 
  scale_linewidth_manual(values = c(rep(1, 3), 2)) + 
  geom_vline(xintercept = cov_pc_df$time[203], linetype = "dotdash", linewidth = 1) + 
  geom_vline(xintercept = cov_pc_df$time[221], linetype = "dotdash", linewidth = 1) + 
  theme(text = element_text(size = 22), legend.key.size = unit(1, "cm"))
graphics.off()

################################

# Individual Conditional Expectations 

# Create fine grids of X values
PC1_seq <- seq(min(X_pc[,1]), max(X_pc[,1]), length.out = 30)
PC2_seq <- seq(min(X_pc[,2]), max(X_pc[,2]), length.out = 30)
PC3_seq <- seq(min(X_pc[,3]), max(X_pc[,3]), length.out = 30)
PC4_seq <- seq(min(X_pc[,4]), max(X_pc[,4]), length.out = 30)

# Start with all PCs fixed at their mean
PC1_ice_X <- tcrossprod(rep(1, 30), colMeans(X_pc))
# vary PC1 
PC1_ice_X[,1] <- PC1_seq

# Start with all PCs fixed at their mean
other_ice_X <- tcrossprod(rep(1, 30 * 3), colMeans(X_pc))

# check the ICEs of PC2, PC3, and PC4 with PC1 fixed at values corresponding 
# to low, medium, and high latitude
PC1_levels <- c(-4, 1, 6)
other_ice_X[,1] <- rep(PC1_levels, each = 30)

PC2_ice_X <- PC3_ice_X <- PC4_ice_X <- other_ice_X

PC2_ice_X[,2] <- rep(PC2_seq, times = 3)
PC3_ice_X[,3] <- rep(PC3_seq, times = 3)
PC4_ice_X[,4] <- rep(PC4_seq, times = 3)

set.seed(0)
PC1_ice_X_nl <- X_hidden(PC1_ice_X, 9, n.h, 0.5)
set.seed(0)
PC2_ice_X_nl <- X_hidden(PC2_ice_X, 9, n.h, 0.5)
set.seed(0)
PC3_ice_X_nl <- X_hidden(PC3_ice_X, 9, n.h, 0.5)
set.seed(0)
PC4_ice_X_nl <- X_hidden(PC4_ice_X, 9, n.h, 0.5)

PC1_lin_ice <- predict(linear_best, newx = PC1_ice_X, logits = FALSE)
PC1_nl_ice <- predict(nl_best, newx = PC1_ice_X_nl, logits = FALSE)
PC2_lin_ice <- predict(linear_best, newx = PC2_ice_X, logits = FALSE)
PC2_nl_ice <- predict(nl_best, newx = PC2_ice_X_nl, logits = FALSE)
PC3_lin_ice <- predict(linear_best, newx = PC3_ice_X, logits = FALSE)
PC3_nl_ice <- predict(nl_best, newx = PC3_ice_X_nl, logits = FALSE)
PC4_lin_ice <- predict(linear_best, newx = PC4_ice_X, logits = FALSE)
PC4_nl_ice <- predict(nl_best, newx = PC4_ice_X_nl, logits = FALSE)

PC1_lin_probs <- PC1_lin_ice$prob[,c(10, 1, 5, 7)]
PC1_nl_probs <- PC1_nl_ice$prob[,c(8, 3, 7, 5)]
colnames(PC1_lin_probs) <- c("pro", "syn", "pico1", "pico2")
colnames(PC1_nl_probs) <- c("pro", "syn", "pico1", "pico2")
PC1_lin_probs <- as.data.frame(PC1_lin_probs)
PC1_nl_probs <- as.data.frame(PC1_nl_probs)
PC1_lin_probs$PC1 <- PC1_nl_probs$PC1 <- PC1_seq
PC1_lin_probs$Response <- PC1_nl_probs$Response <- "Relative_Abundance"
PC1_lin_probs$model <- "linear" 
PC1_nl_probs$model <- "nonlinear"

PC1_lin_diam <- PC1_lin_ice$mn[,1,c(10, 1, 5, 7)]
PC1_nl_diam <- PC1_nl_ice$mn[,1,c(8, 3, 7, 5)]
colnames(PC1_lin_diam) <- c("pro", "syn", "pico1", "pico2")
colnames(PC1_nl_diam) <- c("pro", "syn", "pico1", "pico2")
PC1_lin_diam <- as.data.frame(PC1_lin_diam)
PC1_nl_diam <- as.data.frame(PC1_nl_diam)
PC1_lin_diam$PC1 <- PC1_nl_diam$PC1 <- PC1_seq
PC1_lin_diam$Response <- PC1_nl_diam$Response <- "Log_Diameter"
PC1_lin_diam$model <- "linear" 
PC1_nl_diam$model <- "nonlinear"

PC1_lin_chl <- PC1_lin_ice$mn[,2,c(10, 1, 5, 7)]
PC1_nl_chl <- PC1_nl_ice$mn[,2,c(8, 3, 7, 5)]
colnames(PC1_lin_chl) <- c("pro", "syn", "pico1", "pico2")
colnames(PC1_nl_chl) <- c("pro", "syn", "pico1", "pico2")
PC1_lin_chl <- as.data.frame(PC1_lin_chl)
PC1_nl_chl <- as.data.frame(PC1_nl_chl)
PC1_lin_chl$PC1 <- PC1_nl_chl$PC1 <- PC1_seq
PC1_lin_chl$Response <- PC1_nl_chl$Response <- "Log_chl"
PC1_lin_chl$model <- "linear" 
PC1_nl_chl$model <- "nonlinear"

PC1_lin_pe <- PC1_lin_ice$mn[,3,c(10, 1, 5, 7)]
PC1_nl_pe <- PC1_nl_ice$mn[,3,c(8, 3, 7, 5)]
colnames(PC1_lin_pe) <- c("pro", "syn", "pico1", "pico2")
colnames(PC1_nl_pe) <- c("pro", "syn", "pico1", "pico2")
PC1_lin_pe <- as.data.frame(PC1_lin_pe)
PC1_nl_pe <- as.data.frame(PC1_nl_pe)
PC1_lin_pe$PC1 <- PC1_nl_pe$PC1 <- PC1_seq
PC1_lin_pe$Response <- PC1_nl_pe$Response <- "Log_pe"
PC1_lin_pe$model <- "linear" 
PC1_nl_pe$model <- "nonlinear"

PC1_pred <- rbind(PC1_lin_probs, PC1_nl_probs, 
                  PC1_lin_diam, PC1_nl_diam, 
                  PC1_lin_chl, PC1_nl_chl, 
                  PC1_lin_pe, PC1_nl_pe)

PC1_pred_melt <- melt(PC1_pred, id.vars = c("PC1", "model", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

PC1_pred_melt$Response <- factor(PC1_pred_melt$Response, 
  levels = c("Relative_Abundance", "Log_Diameter", "Log_chl", "Log_pe"))

PC1_plot <- ggplot(PC1_pred_melt) + 
  geom_line(aes(x = -PC1, y = Prediction, linetype = model, color = model), linewidth = 1) + 
  scale_linetype_manual(values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
  facet_wrap(Population ~ Response, scales = "free") + 
  geom_vline(xintercept = -X_pc[203,1], linetype = "dotted") + 
  geom_vline(xintercept = -X_pc[221,1], linetype = "dotted") +
  theme(text = element_text(size = 10))

# Figure 8
pdf(file.path("plots", "PC1_all.pdf"), 10, 8.5)
PC1_plot
graphics.off()

# code to output a list of these plots
PC1_ice_list <- lapply(c("pro", "syn", "pico1", "pico2"), function(pop) {
  lapply(c("Relative_Abundance", "Log_Diameter", "Log_chl", "Log_pe"), function(resp) {
    ggplot(filter(PC1_pred_melt, Population == pop, Response == resp)) + 
  geom_line(aes(x = -PC1, y = Prediction, linetype = model, color = model), linewidth = 1) + 
  scale_linetype_manual(name = "Model", labels = c("Linear", "Nonlinear"), values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
  scale_color_discrete(name = "Model", labels = c("Linear", "Nonlinear")) +
  geom_vline(xintercept = -X_pc[203,1], linetype = "dotted") + 
  geom_vline(xintercept = -X_pc[221,1], linetype = "dotted") +
  theme(text = element_text(size = 22)) + 
  labs(title = paste(pop, resp))
  })
}) %>% unlist(recursive = FALSE)

# PC1

# Figure 6
pro2 <- PC1_ice_list[[4]] + 
  labs(title  = "Pro Phycoerythrin", x = "", y = "") +
  labs(y = "Prediction")
  
leg <- ggpubr::get_legend(pro2)
leg <- ggpubr::as_ggplot(leg)

pro2 <- pro2 + 
  theme(title = element_text(size = 20), 
        axis.title.y = element_text(size = 18), 
        legend.text = element_text(size = 18), legend.title = element_text(size = 20), 
        legend.position = "none") + 
  scale_x_continuous(breaks = seq(-8, 8, by = 4)) + 
  scale_y_continuous(breaks = seq(0.31, 0.39, by = 0.02), limits = c(0.31, 0.39))

pico1 <- PC1_ice_list[[10]] + 
  labs(title  = "Pico1 Log Diameter", x = "") + 
  theme(legend.position = "none", title = element_text(size = 20), 
        axis.title.y = element_text(size = 18), 
        legend.text = element_text(size = 18), legend.title = element_text(size = 20)) + 
  scale_x_continuous(breaks = seq(-8, 8, by = 4)) + 
  scale_y_continuous(breaks = seq(5.4, 5.8, 0.1), limits = c(5.4, 5.8))

pico2 <- PC1_ice_list[[13]] + 
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.04), limits = c(0, 0.2)) + 
  scale_x_continuous(breaks = seq(-8, 8, by = 4)) + 
  labs(title = "Pico2 Relative Abundance", y = "", x = "") + 
  theme(legend.position = "none", title = element_text(size = 20), 
        axis.title.y = element_text(size = 18), 
        legend.text = element_text(size = 18), legend.title = element_text(size = 20))

pro_x_lab <- ggpubr::text_grob("-PC1 (Proxy for Latitude)", size = 18, hjust = 0.375, vjust = -0.5)
pdf(file.path("plots", "PC1_paper_plot.pdf"), 20.18, 6.65)
grid.arrange(pro2, pico1, pico2, leg, pro_x_lab, nrow = 2, 
            layout_matrix = matrix(c(1, 2, 3, 4, 
                                     NA, 5, NA, NA), 2, byrow = TRUE), 
            widths = c(1, 1, 1, 0.25), heights = c(1, 0.05))
graphics.off()

# other PCs

plot_ice <- function(PC_lin_ice, PC_nl_ice, PC_num, PC_seq, return_plot_list = FALSE) {
  PC_lin_probs <- PC_lin_ice$prob[,c(10, 1, 5, 7)]
  PC_nl_probs <- PC_nl_ice$prob[,c(8, 3, 7, 5)]
  colnames(PC_lin_probs) <- c("pro", "syn", "pico1", "pico2")
  colnames(PC_nl_probs) <- c("pro", "syn", "pico1", "pico2")
  PC_lin_probs <- as.data.frame(PC_lin_probs)
  PC_nl_probs <- as.data.frame(PC_nl_probs)
  PC_lin_probs$PC1 <- PC_nl_probs$PC1 <- factor(other_ice_X[,1], levels = c(-4, 1, 6), 
                                                labels = c("Subpolar", 
                                                           "Transition Zone", 
                                                           "Subtropical"))
  PC_lin_probs[,PC_num] <- PC_nl_probs[,PC_num] <- PC_seq
  PC_lin_probs$Response <- PC_nl_probs$Response <- "Relative_Abundance"
  PC_lin_probs$model <- "linear" 
  PC_nl_probs$model <- "nonlinear"

  PC_lin_diam <- PC_lin_ice$mn[,1,c(10, 1, 5, 7)]
  PC_nl_diam <- PC_nl_ice$mn[,1,c(8, 3, 7, 5)]
  colnames(PC_lin_diam) <- c("pro", "syn", "pico1", "pico2")
  colnames(PC_nl_diam) <- c("pro", "syn", "pico1", "pico2")
  PC_lin_diam <- as.data.frame(PC_lin_diam)
  PC_nl_diam <- as.data.frame(PC_nl_diam)
  PC_lin_diam$PC1 <- PC_nl_diam$PC1 <- factor(other_ice_X[,1], levels = c(-4, 1, 6), labels = c("Subpolar", 
                                                                              "Transition Zone", 
                                                                              "Subtropical"))
  PC_lin_diam[,PC_num] <- PC_nl_diam[,PC_num] <- PC_seq
  PC_lin_diam$Response <- PC_nl_diam$Response <- "Log_Diameter"
  PC_lin_diam$model <- "linear" 
  PC_nl_diam$model <- "nonlinear"

  PC_lin_chl <- PC_lin_ice$mn[,2,c(10, 1, 5, 7)]
  PC_nl_chl <- PC_nl_ice$mn[,2,c(8, 3, 7, 5)]
  colnames(PC_lin_chl) <- c("pro", "syn", "pico1", "pico2")
  colnames(PC_nl_chl) <- c("pro", "syn", "pico1", "pico2")
  PC_lin_chl <- as.data.frame(PC_lin_chl)
  PC_nl_chl <- as.data.frame(PC_nl_chl)
  PC_lin_chl$PC1 <- PC_nl_chl$PC1 <- factor(other_ice_X[,1], levels = c(-4, 1, 6), labels = c("Subpolar", 
                                                                              "Transition Zone", 
                                                                              "Subtropical"))
  PC_lin_chl[,PC_num] <- PC_nl_chl[,PC_num] <- PC_seq
  PC_lin_chl$Response <- PC_nl_chl$Response <- "Log_chl"
  PC_lin_chl$model <- "linear" 
  PC_nl_chl$model <- "nonlinear"

  PC_lin_pe <- PC_lin_ice$mn[,3,c(10, 1, 5, 7)]
  PC_nl_pe <- PC_nl_ice$mn[,3,c(8, 3, 7, 5)]
  colnames(PC_lin_pe) <- c("pro", "syn", "pico1", "pico2")
  colnames(PC_nl_pe) <- c("pro", "syn", "pico1", "pico2")
  PC_lin_pe <- as.data.frame(PC_lin_pe)
  PC_nl_pe <- as.data.frame(PC_nl_pe)
  PC_lin_pe$PC1 <- PC_nl_pe$PC1 <- factor(other_ice_X[,1], levels = c(-4, 1, 6), labels = c("Subpolar", 
                                                                              "Transition Zone", 
                                                                              "Subtropical"))
  PC_lin_pe[,PC_num] <- PC_nl_pe[,PC_num] <- PC_seq
  PC_lin_pe$Response <- PC_nl_pe$Response <- "Log_pe"
  PC_lin_pe$model <- "linear" 
  PC_nl_pe$model <- "nonlinear"

  PC_pred <- rbind(PC_lin_probs, PC_nl_probs, 
                    PC_lin_diam, PC_nl_diam, 
                    PC_lin_chl, PC_nl_chl, 
                    PC_lin_pe, PC_nl_pe)

  PC_pred_melt <- melt(PC_pred, id.vars = c("PC1", PC_num, "model", "Response"), 
                                    variable.name = "Population", 
                                    value.name = "Prediction")

  PC_pred_melt$Response <- factor(PC_pred_melt$Response, 
    levels = c("Relative_Abundance", "Log_Diameter", "Log_chl", "Log_pe"))

  subsampled_pts <- PC_pred_melt %>% 
      group_split(Population, Response, model, PC1) %>% 
      lapply(function(cur_df) {
        cur_seq <- seq(1, nrow(cur_df), length.out = 6) %>% round()
        cur_seq <- cur_seq[2:5]
        cur_df[cur_seq,]
      }) %>% bind_rows()

  if(return_plot_list) {
    lapply(c("pro", "syn", "pico1", "pico2"), function(pop) {
      lapply(c("Relative_Abundance", "Log_Diameter", "Log_chl", "Log_pe"), function(resp) {
        ggplot(filter(PC_pred_melt, Response == resp, Population == pop)) + 
          geom_line(aes(x = .data[[PC_num]], y = Prediction, linetype = model, color = PC1), linewidth = 1) + 
          geom_point(data = filter(subsampled_pts, Response == resp, Population == pop), aes(x = .data[[PC_num]], y = Prediction, group = model, color = PC1, shape = PC1), size = 4) + 
          scale_linetype_manual(name = "Model", labels = c("Linear", "Nonlinear"), values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
          # facet_wrap(Population ~ Response, scales = "free")
          labs(title = pop, y = resp) + 
          theme(text = element_text(size = 20))
      })
    }) %>% unlist(recursive = FALSE)
  } else { # put four points on lines to distinguish if plots are printed in greyscale
    ggplot(PC_pred_melt) + 
          geom_line(aes(x = .data[[PC_num]], y = Prediction, linetype = model, color = PC1), linewidth = 1) + 
          geom_point(data = subsampled_pts, aes(x = .data[[PC_num]], y = Prediction, group = model, color = PC1, shape = PC1), size = 2.5) + 
          scale_linetype_manual(name = "Model", labels = c("Linear", "Nonlinear"), values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
          facet_wrap(Population ~ Response, scales = "free") + 
          theme(text = element_text(size = 10))
  }
}

# PC2

PC2_plot <- plot_ice(PC2_lin_ice, PC2_nl_ice, "PC2", PC2_seq)
PC2_plist <- plot_ice(PC2_lin_ice, PC2_nl_ice, "PC2", PC2_seq, return_plot_list = TRUE)

# Figure 9
pdf(file.path("plots", "PC2_all.pdf"), 10, 8.5)
PC2_plot
graphics.off()

# pro abundance
PC2_1 <- PC2_plist[[1]] + 
  labs(title = "Pro Relative Abundance", y = "Prediction") + 
  theme(legend.position = "none", title = element_text(size = 20), 
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 16), 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18)) + 
  labs(x = "", y = "Prediction")

# syn pe
PC2_2 <- PC2_plist[[8]] + 
  labs(title = "Syn Phycoerythrin", y = "Prediction") + 
  theme(legend.position = "none", title = element_text(size = 20), 
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 16), 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18)) + 
  labs(x = "", y = "")

# pico2 abundance
PC2_3 <- PC2_plist[[13]] + 
  labs(title = "Pico2 Relative Abundance", y = "Prediction") + 
  labs(x = "", y = "") + 
  theme(legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18), 
        legend.key.size = unit(1, "cm"))

PC2_leg <- ggpubr::get_legend(PC2_3) %>% 
  ggpubr::as_ggplot()

PC2_3 <- PC2_3 + 
  theme(legend.position = "none", title = element_text(size = 20), 
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 16))

PC2_xlab <- ggpubr::text_grob("PC2", hjust = 2.5, vjust = -0.5, size = 18)

# Figure 7
pdf(file.path("plots", "PC2_selected.pdf"), 21.3, 6.7)
grid.arrange(PC2_1, PC2_2, PC2_3, PC2_leg, PC2_xlab, 
             layout_matrix = matrix(c(1, 2, 3, 4, 
                                       5, 5, 5, 5), nrow = 2, byrow = TRUE), 
             widths = c(1, 1, 1, 0.4), heights = c(1, 0.05))
graphics.off()

# Figures 14-17

pop_mat <- matrix(c(10, 8, 
                    1, 3, 
                    5, 7,
                    7, 5), ncol = 2, byrow = TRUE)

models <- c("Linear", "Nonlinear")
pops <- c("Pro", "Syn", "Pico1",
          "Pico2")

colnames(pop_mat) <- models
rownames(pop_mat) <- pops

# response vs time
plot_1d_resp <- function(pop_mat, pop, resp) {
  resp_dat <- switch(resp, 
                      prob = "$prob[,", 
                      diam = "$mn[,1,", 
                      chl = "$mn[,2,", 
                      pe = "$mn[,3,")
  
  clusts <- pop_mat[pop,]

  lin_resp <- linear_best %>% 
                substitute() %>% 
                deparse() %>%
                paste0(resp_dat, clusts["Linear"], "]") %>% 
                str2lang() %>%
                eval()
  half_resp <- half_best %>% 
                substitute() %>% 
                deparse() %>%
                paste0(resp_dat, clusts["Nonlinear"], "]") %>% 
                str2lang() %>%
                eval()
  
  lin_df <- data.frame(time = time, lin_resp, model = "Linear")
  half_df <- data.frame(time = time, half_resp, model = "Nonlinear")
  colnames(lin_df)[2] <- resp 
  colnames(half_df)[2] <- resp 

  resp_df <- rbind(lin_df, half_df)

  plot_title <- paste0(pop, ": Linear-", clusts["Linear"], 
                      ", Nonlinear-", clusts["Nonlinear"], ": ", resp)
  ggplot(resp_df, aes(x = time, y = .data[[resp]], group = model, color = model)) + 
    geom_line() + 
    ggtitle(plot_title) + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title = element_text(size = 16), 
          legend.key.height = unit(1, 'cm'), 
          legend.key.width = unit(1.5, "cm"), 
          legend.text = element_text(size=16), 
          legend.title = element_text(size=16), 
          title = element_text(size = 16))
}



for(pop in pops) {
  fname <- paste0("response_time_", pop, ".png")
  fpath <- file.path("plots", fname)
  p1 <- plot_1d_resp(pop_mat, pop, "prob")
  p2 <- plot_1d_resp(pop_mat, pop, "diam")
  p3 <- plot_1d_resp(pop_mat, pop, "chl")
  p4 <- plot_1d_resp(pop_mat, pop, "pe")

  pdf(fpath, 17, 9.6)
  print(grid.arrange(p1, p2, p3, p4, nrow = 2))
  graphics.off()
}