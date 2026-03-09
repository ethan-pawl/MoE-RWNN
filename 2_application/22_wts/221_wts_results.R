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

plots_dir <- "plots"
if(!dir.exists(plots_dir)) dir.create(plots_dir)

# Load all data
datobj1 <- readRDS(file = file.path("data", "MGL1704-hourly-paper-counts-and-biomass.RDS"))

# Collect results
nl_cv_file <- file.path(
  "2_application", 
  "21_all_data", 
  "results", 
  "nl_9_70_4", 
  "nl_9_70_4_pcr_summary.RDS"
)

nl_no_wts_file <- file.path(
  "2_application", 
  "22_wts", 
  "results", 
  "nl_no_wts_9_70_4", 
  "nl_no_wts_9_70_4_wts_summary.RDS"
)

nl_cv <- readRDS(nl_cv_file)
nl_cv_no_wts <- readRDS(nl_no_wts_file)

nl_best <- nl_cv$bestres
nl_no_wts_best <- nl_cv_no_wts$bestres

##########################

plot_list_nl <- flowtrend::plot_3d(
  datobj1$ybin_list, nl_best, 33, datobj1$biomass_list, return_list_of_plots = TRUE, 
  labels = c(
    "Pico2", "Pico1", "Other5", "Bead", "Other1", "Syn", "Other6", "Other2", "Pro", "Other4"
  ), 
  mn_colours = c(
    "red", "red", "gray", "gray", "gray", "red", "gray", "gray", "red", "gray"
  )
)
p1 <- plot_list_nl[[1]]

plot_list_nl_no_wts <- flowtrend::plot_3d( # Still plot the biomass even though the model doesn't see it
  datobj1$ybin_list, nl_no_wts_best, 33, datobj1$biomass_list, return_list_of_plots = TRUE, 
  labels = c(
    "Syn", "Other8", "Pico", "Pro1", "Other7", "Other1", "Pro6", "Pro5", "Pro2", "Pro3"
  ), 
  mn_colours = c(
    "red", "gray", "red", "red", "gray", "gray", "red", "red", "red", "red"
  )
)
p2 <- plot_list_nl_no_wts[[1]]

p1 <- p1 + 
  labs(title = "Estimated Nonlinear Pseudolikelihood Model", y = "Chlorophyll", x = "") + 
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) + 
  scale_x_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) + 
  scale_y_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) 

p2 <- p2 + 
  labs(title = "Estimated Nonlinear Likelihood Model", y = "", x = "") + 
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) + 
  scale_x_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) + 
  scale_y_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) 

x_lab <- text_grob("Diameter", size = 12, x = 0.515, y = 1)
overall_title <- text_grob("Estimated Models\n2017-07-02 08:00-09:00", 
                           size = 18)

pdf(file.path(plots_dir, "Figure-wts-comp.pdf"), 10.5, 6)
print(grid.arrange(p1, p2, x_lab, overall_title, layout_matrix = matrix(c(4, 4, 
                                                                    1, 2, 
                                                                    3, 3), 
                                                                    byrow = TRUE, 
                                                                    nrow = 3), 
             heights = c(0.175, 1, 0.05), widths = c(1, 1, 0.05)))
graphics.off()
