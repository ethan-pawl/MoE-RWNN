library(flowmix)

# Load the data
datobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper-counts-and-biomass.RDS"))

ylist <- datobj$ybin_list 
countslist <- datobj$countslist

X_dir <- file.path(
  "data", 
  "X_variations"
)

X <- readRDS(file.path(X_dir, "X_pc_9_nh_70_seed_4_ofold_NA_ifold_NA.RDS"))

load(file.path("2_application", "22_wts", "settings.Rdata"))

nl_no_wts_best <- flowmix_once(
  ylist = ylist, 
  countslist = countslist, 
  X = X, 
  prob_lambda = prob_lambda, 
  mean_lambda = mean_lambda, 
  maxdev = 0.5, 
  numclust = 10, 
  seed = seed, 
  verbose = TRUE
)

load(file.path("2_application", "1_all_data", "nl_fit.Rdata"))
nl_best <- res

library(flowtrend)
library(ggplot2)
library(ggpubr)
library(gridExtra)

plot_list_nl <- flowtrend::plot_3d(
  datobj$ybin_list, nl_best, 33, datobj$biomass_list, return_list_of_plots = TRUE, 
  labels = c(
    "Pico2", "Pico1", "Other5", "Bead", "Other1", "Syn", "Other6", "Other2", "Pro", "Other4"
  ), 
  mn_colours = c(
    "red", "red", "gray", "gray", "gray", "red", "gray", "gray", "red", "gray"
  )
)
p1 <- plot_list_nl[[1]]

plot_list_nl_no_wts <- flowtrend::plot_3d( # Still plot the biomass even though the model doesn't see it
  datobj$ybin_list, nl_no_wts_best, 33, datobj$biomass_list, return_list_of_plots = TRUE, 
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

pdf(file.path("plots", "SuppFigure14.pdf"), 10.5, 6)
print(grid.arrange(p1, p2, x_lab, overall_title, layout_matrix = matrix(c(4, 4, 
                                                                    1, 2, 
                                                                    3, 3), 
                                                                    byrow = TRUE, 
                                                                    nrow = 3), 
             heights = c(0.175, 1, 0.05), widths = c(1, 1, 0.05)))
graphics.off()