library(ggpubr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

datobj <- readRDS(file.path("data", "MGL1704-hourly-paper.RDS"))
datobj %>% 
  list2env(envir = .GlobalEnv) %>% 
  invisible()

# Remove the changepoint variables
X_cont <- X[,3:39]

# Some exploratory analysis on the entire dataset

# PCA; data is already centered and scaled
pr.out <- prcomp(X_cont, center = FALSE, scale = FALSE)
pr.var <- pr.out$sdev**2
pve <- pr.var / sum(pr.var)
cpve <- cumsum(pve)

# Calculate the q-dimensional projection
# of X into the principal components space
q_hat <- which(cpve >= 0.95)[1] # 9
save(q_hat, file = file.path("data", "q_hat.Rdata"))
X_pc <- X_cont %*% pr.out$rotation[,1:q_hat]

# Bind latitude to get correlations b/w latitude, 
# principal components, and covariates
lat_X_cont <- cbind(lat, X_cont)
lat_X_pc <- cbind(lat, X_pc)
X_pc_cor <- cor(lat_X_cont, lat_X_pc)

# Reformat the data frame for use with 
# ggplot2
cor_df <- melt(X_pc_cor)

# Correlation heatmap (top panel) in Figure 1
# Order the covariates and give them cleaner names
ordered_covs <- c(
  "lat",
  "silicate_WOA_clim",
  "phosphate_WOA_clim",
  "nitrate_WOA_clim",
  "oxygen_WOA_clim",
  "o2sat_WOA_clim",
  "density_WOA_clim",
  "O2",
  "PO4",
  "PHYC",
  "CHL",
  "NO3",
  "Si",
  "PP",
  "surface_downward_eastward_stress",
  "eastward_wind",
  "conductivity_WOA_clim",
  "salinity_WOA_clim",
  "AOU_WOA_clim",
  "disp_bw_sla",
  "ftle_bw_sla",
  "vgos",
  "vgosa",
  "Fe",
  "sst",
  "sss",
  "northward_wind",
  "surface_downward_northward_stress",
  "wind_speed",
  "wind_stress",
  "sla",
  "ugos",
  "ugosa",
  "par",
  "p1",
  "p2",
  "p3",
  "p4"
)

ordered_covs_labs <- c(
  "Latitude",
  "Silicate",
  "Phosphate",
  "Nitrate",
  "Oxygen",
  "O2Sat",
  "Sensity",
  "O2",
  "PO4",
  "PHYC",
  "CHL",
  "NO3",
  "Si",
  "PP",
  "SDES",
  "Eastward Wind",
  "Conductivity",
  "Salinity",
  "AOU",
  "Displacement",
  "FTLE",
  "VGOS",
  "VGOSA",
  "Fe",
  "SST",
  "SSS",
  "Northward Wind",
  "SDNS",
  "Wind Speed",
  "Wind Stress",
  "SLA",
  "UGOS",
  "UGOSA",
  "PAR",
  "p1",
  "p2",
  "p3",
  "p4"
)

cor_df$Var1 <- factor(
  cor_df$Var1, 
  levels = ordered_covs, 
  labels = ordered_covs_labs
)

cor_df$Var2 <- factor(
  cor_df$Var2, 
  levels = c(paste0("PC", 9:1), "lat")
)

# Initial correlation heatmap
corplot <- ggplot(cor_df) + 
  geom_tile(aes(Var1, Var2, fill = value)) + 
  scale_fill_viridis_c(
    name = "Correlation", 
    option = "cividis", 
    limits = c(-1, 1)
  ) + 
  theme(
    plot.title = element_text(size = 20),
    text = element_text(size = 18), 
    axis.text.x = element_text(angle = 35, hjust = 1, size = 16), 
    legend.title = element_text(size = 20, margin = margin(b = 15)), 
    legend.key.size = unit(1, "cm")
  ) + 
  labs(
    x = "", 
    y = "", 
    title = "Correlations Between Covariates, Principal Components, and Latitude"
  )

# Strongest interpretations for PCs 1-4
# More interpretable if we flip the signs

X_pc[,1] <- -X_pc[,1]
X_pc[,2] <- -X_pc[,2]
X_pc[,3] <- -X_pc[,3]
X_pc[,4] <- -X_pc[,4]

# Bind latitude to get correlations b/w latitude, 
# principal components, and covariates
lat_X_cont <- cbind(lat, X_cont)
lat_X_pc <- cbind(lat, X_pc)
X_pc_cor <- cor(lat_X_cont, lat_X_pc)

# Reformat the data frame for use with 
# ggplot2
cor_df <- melt(X_pc_cor)

# Filter, rename, and order some of the covariates
# in the plot
reduced_cor_df <- cor_df %>%
  filter(Var2 %in% c("lat", "PC1", "PC2", "PC3", "PC4"))

reduced_cor_df$Var1 <- factor(
  reduced_cor_df$Var1, 
  levels = ordered_covs, 
  labels = ordered_covs_labs
)

reduced_cor_df$Var2 <- factor(
  reduced_cor_df$Var2, 
  levels = c(paste0("PC", 4:1), "lat")
)

# Correlation heatmap (top panel of Figure 1)
corplot <- ggplot(reduced_cor_df) + 
  geom_tile(aes(Var1, Var2, fill = value)) + 
  scale_fill_viridis_c(
    name = "Correlation", 
    option = "cividis", 
    limits = c(-1, 1)
  ) + 
  theme(
    plot.title = element_text(size = 20),
    text = element_text(size = 18), 
    axis.text.x = element_text(angle = 35, hjust = 1, size = 16), 
    legend.title = element_text(size = 20, margin = margin(b = 15)), 
    legend.key.size = unit(1, "cm")
  ) + 
  labs(
    x = "", 
    y = "", 
    title = "Correlations Between Covariates, Principal Components, and Latitude"
  )

# Line plot of PC1 and important covariates (bottom-left panel in Figure 1)

# Convert timestamps to datetime objects
dt <- as.POSIXct(time, format = "%FT%T", tz = "UTC")

# Gather data in a data frame
pc_df <- data.frame(
  Time = rep(dt, 9), 
  Data = as.vector(X_pc), 
  Variable = rep(paste0("PC", 1:9), each = 296), 
  Focus = rep(c(TRUE, FALSE), times = c(296 * 4, 296 * 5)), 
  Modeled = TRUE,
  isPC = TRUE
)

# Give cleaner names
short_names <- colnames(X_cont)
names(short_names) <- short_names 
short_names["AOU_WOA_clim"] <- "AOU"
short_names["phosphate_WOA_clim"] <- "Phosphate"
short_names["conductivity_WOA_clim"] <- "Conductivity"
short_names["density_WOA_clim"] <- "Density"
short_names["salinity_WOA_clim"] <- "Salinity"
short_names["nitrate_WOA_clim"] <- "Nitrate"
short_names["silicate_WOA_clim"] <- "Silicate"
short_names["o2sat_WOA_clim"] <- "O2Sat"
short_names["oxygen_WOA_clim"] <- "Oxygen"
short_names["surface_downward_eastward_stress"] <- "SDES"
short_names["surface_downward_northward_stress"] <- "SDNS"
short_names["sst"] <- "SST"
short_names["sss"] <- "SSS"
short_names["par"] <- "PAR"
short_names["sla"] <- "SLA"
short_names["wind_speed"] <- "Wind Speed"

# Label which PCs each covariates "belongs to"
which_pc <- c(
  "PC3", "PC3", "PC4", "PC4", "PC1", "PC1", "PC2", "PC1", "PC1", "PC1", 
  "PC1", "PC1", "PC1", "PC1", "PC2", "PC2", "PC2", "PC3", "PC3", "PC2", 
  "PC1", "PC1", "PC2", "PC2", "PC2", "PC1", "PC1", "PC1", "PC1", "PC1", 
  "PC1", "PC1", "PC1", "PC1", "PC1", "PC1", "PC4"
)
names(which_pc) <- short_names

# Tell which covariates we want to plot
focus <- c(
  TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, 
  FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE
)
names(focus) <- short_names

# Gather covariates into a data frame
cov_df <- data.frame(
  Time = rep(dt, ncol(X_cont)), 
  Data = as.vector(X_cont), 
  Variable = rep(short_names, each = 296), 
  Modeled = TRUE,
  isPC = FALSE
)

# Transform latitude to be roughly on the same scale as PC1
trans_lat <- (lat - min(lat)) / (max(lat) - min(lat)) * 15 -8.5

# Gather latitude into a data frame
lat_df <- data.frame(
  Time = dt, 
  Data = trans_lat, 
  Variable = "Latitude", 
  PC = "PC1", 
  Focus = TRUE, 
  Modeled = FALSE, 
  isPC = FALSE
)

# Make the data frames have the same columns
cov_df$Focus <- focus[cov_df$Variable]
pc_df$PC <- pc_df$Variable
cov_df$PC <- which_pc[cov_df$Variable]

# Concatenate into one big data frame for plotting
big_df <- rbind(pc_df, cov_df, lat_df)

# Convert certain columns into factors
big_df$PC <- factor(big_df$PC, levels = paste0("PC", 1:9))
big_df$Variable <- factor(
  big_df$Variable, 
  levels = c(paste0("PC", 1:9), short_names[c(37, 1:36)], "Latitude")
)

# Assign colors, line widths, and colors to variables, PCs, and latitude
brewer_cols <- brewer.pal(12, "Set3")

named_colors <- c(
  "black", 
  brewer_cols[4], 
  brewer_cols[5], 
  brewer_cols[10], 
  brewer_cols[6], 
  "black"
)

named_widths <- c(
  1.5, 
  1, 
  1, 
  1, 
  1, 
  1
)

named_types <- c(rep("solid", 5), "dashed")

names(named_colors) <- names(named_widths) <- names(named_types) <- c(
  "PC1", 
  "SSS", 
  "SST", 
  "NO3", 
  "PO4", 
  "Latitude"
)

PC1_plot <- ggplot(filter(big_df, PC == "PC1", Focus)) + 
  geom_line(
    aes(
      Time, 
      Data, 
      color = Variable, 
      linetype = Variable, 
      linewidth = Variable
    )
  ) + 
  scale_linetype_manual(values = named_types) + 
  scale_color_manual(values = named_colors) + 
  scale_linewidth_manual(values = named_widths) + 
  scale_x_datetime(date_breaks = "2 days", date_labels = "%b %d") + 
  theme(
    plot.title = element_text(size = 20),
    legend.position = "bottom",
    legend.title = element_text(
      size = 18, 
      hjust = 0.5, 
      margin = margin(l = 10, b = 10)
    ),
    legend.text = element_text(size = 16, margin = margin(l = 5, r = 10)),
    legend.key.size = unit(2, "lines"), 
    text = element_text(size = 18), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14), 
    axis.text.x = element_text(angle = 35, hjust = 1), 
    legend.box.margin = margin(r = 10)) + 
  labs(title = "Principal Component 1 and Covariates over Time", x = "") + 
  guides(linetype = guide_legend(title.position = "top")) + 
  geom_vline(xintercept = dt[221], linetype = "dotted") + 
  geom_vline(xintercept = dt[204], linetype = "dotted")

# PC3 & PC4
# Repeat similar process as with PC1_plot
colors_34 <- c(
  "black", 
  brewer_cols[1],
  brewer_cols[4], 
  brewer_cols[5], 
  brewer_cols[6], 
  brewer_cols[7],
  brewer_cols[10]
)

widths_34 <- c(
  1.5, 
  1.5, 
  0.5, 
  0.5, 
  0.5, 
  0.5, 
  0.5
)

alphas_34 <- c(
  1, 
  1, 
  0.7, 
  0.65, 
  0.65, 
  0.65, 
  0.65
)

names(colors_34) <- names(widths_34) <- names(alphas_34) <- c(
  "PC3", 
  "PC4", 
  "PAR", 
  "p1", 
  "p2", 
  "p3", 
  "p4"
)

PC34_plot <- ggplot(filter(big_df, PC %in% c("PC3", "PC4"), Focus)) + 
  geom_line(
    aes(
      Time, 
      Data, 
      color = Variable, 
      linewidth = Variable, 
      alpha = Variable
    )
  ) + 
  scale_color_manual(values = colors_34) + 
  scale_linewidth_manual(values = widths_34) + 
  scale_alpha_manual(values = alphas_34) + 
  scale_x_datetime(date_breaks = "2 days", date_labels = "%b %d") + 
  theme(
    plot.title = element_text(size = 20),
    legend.position = "bottom",
    legend.title = element_text(
      size = 18, 
      hjust = 0.5, 
      margin = margin(l = 5, b = 10)
    ),
    legend.text = element_text(size = 16, margin = margin(l = 5, r = 10)), 
    legend.key.size = unit(2, "lines"), 
    text = element_text(size = 18), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14), 
    axis.text.x = element_text(angle = 35, hjust = 1), 
    legend.box.margin = margin(r = 10)) + 
  labs(
    title = "Principal Components 3 & 4 and Covariates over Time", 
    x = "", 
    y = ""
  ) + 
  guides(
    color = guide_legend(
      title.position = "top", 
      override.aes = list(
        linewidth = c(
          1.5, 
          1.5, 
          1, 
          1, 
          1, 
          1, 
          1
        )
      )
    )
  )

# Extract legends
leg1 <- get_legend(PC1_plot) %>% 
  as_ggplot()

PC1_plot <- PC1_plot + 
  theme(legend.position = "none")

leg34 <- get_legend(PC34_plot) %>% 
  as_ggplot()
PC34_plot <- PC34_plot + 
  theme(legend.position = "none")

if(!dir.exists("plots")) dir.create("plots")

# Figure 1
pdf(file.path("plots", "Figure01.pdf"), 18.4, 10)
gridExtra::grid.arrange(corplot, PC1_plot, PC34_plot, leg1, leg34, 
                        layout_matrix = matrix(c(1, 1, 
                                                 2, 3,
                                                 5, 6), byrow = TRUE, 
                                                 nrow = 3), 
                        heights = c(1, 1.025, 0.4))
graphics.off()

# Appendix Figure 10
colors_2 <- c(
  "black", 
  brewer_cols[4], 
  brewer_cols[5], 
  brewer_cols[10]
)

widths_2 <- c(
  1.5, 
  1, 
  1, 
  1
)

names(colors_2) <- names(widths_2) <- c(
  "PC2", 
  "Fe", 
  "SLA", 
  "Wind Speed"
)

PC2_plot <- ggplot(
  filter(
    big_df, 
    PC == "PC2", 
    Variable %in% c("PC2", "Fe", "SLA", "Wind Speed")
  )
) + 
  geom_line(aes(Time, Data, color = Variable, linewidth = Variable)) + 
  scale_color_manual(values = colors_2) + 
  scale_linewidth_manual(values = widths_2) + 
  scale_x_datetime(date_breaks = "2 days", date_labels = "%b %d") + 
  theme(
    plot.title = element_text(size = 20),
    legend.position = "bottom",
    legend.title = element_text(size = 18, hjust = 0.5, margin = margin(b = 10)),
    legend.text = element_text(size = 16, margin = margin(r = 10)),
    legend.key.size = unit(2, "lines"), 
    text = element_text(size = 18), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14), 
    axis.text.x = element_text(angle = 35, hjust = 1)) + 
  labs(title = "Principal Component 2 and Covariates over Time", x = "") + 
  guides(color = guide_legend(title.position = "top"))

pdf(file.path("plots", "Figure10.pdf"), 9.2, 5.88)
PC2_plot 
graphics.off()

# Appendix Figure 11
time_dictionary <- data.frame(date = dt, time = 1:296)

all_pc_plot <-
  X_pc %>% as_tibble() %>% mutate(time = 1:nrow(X_pc)) %>%
  left_join(time_dictionary, by = "time") %>%
  select(-time) %>%
  pivot_longer(-date) %>%
  ggplot() + 
  geom_line(aes(x = date, y = value, group = name), linewidth = 1.1) +
  facet_wrap(~name) +
  ylab("") +
  xlab("") +
  scale_x_datetime(date_breaks = "2 day", date_labels = "%B %d")  +
  theme(
    plot.title = element_text(size = 20),
    text = element_text(size = 18), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14), 
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  labs(title = "All Principal Components")

pdf(file.path("plots", "Figure11.pdf"), width = 12.5, height = 11.4 * .8)
all_pc_plot
graphics.off()
