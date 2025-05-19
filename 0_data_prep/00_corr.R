library(magrittr)
library(ggplot2)
library(reshape2)
library(dplyr)

datobj <- readRDS(file.path("data", "MGL1704-hourly-paper.RDS"))
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()

# Remove the artificial changepoint variables
X_cont <- X[,3:39]

# PCA; data is already centered and scaled
pr.out <- prcomp(X_cont, center = FALSE, scale = FALSE)
pr.var <- pr.out$sdev**2
pve <- pr.var / sum(pr.var)
cpve <- cumsum(pve)

# Calculate the q-dimensional approximation of the projection
# of X into the principal components space
q_hat <- which(cpve >= 0.95)[1]
X_pc <- X_cont %*% pr.out$rotation[,1:q_hat]

# Save the principal components representation 
# to be used during modeling
X <- X_pc 
save(X, file = file.path("data", "X_pc.Rdata"))

# Create the random-weight neural network-transformed 
# representation of the covariates (the hidden layer 
# outputs)
n.h <- 70
a <- 0.5
X_nl <- flowmix::make_hidden_nodes(X, n.h, a)
X <- X_nl 
save(X, file = file.path("data", "X_nl.Rdata"))

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

ordered_covs <- c("lat",
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
                  "p4")

ordered_covs_labs <- c("lat",
                       "silicate",
                       "phosphate",
                       "nitrate",
                       "oxygen",
                       "o2sat",
                       "density",
                       "O2",
                       "PO4",
                       "PHYC",
                       "CHL",
                       "NO3",
                       "Si",
                       "PP",
                       "sdes",
                       "eastward_wind",
                       "conductivity",
                       "salinity",
                       "AOU",
                       "disp_bw_sla",
                       "ftle_bw_sla",
                       "vgos",
                       "vgosa",
                       "Fe",
                       "sst",
                       "sss",
                       "northward_wind",
                       "sdns",
                       "wind_speed",
                       "wind_stress",
                       "sla",
                       "ugos",
                       "ugosa",
                       "par",
                       "p1",
                       "p2",
                       "p3",
                       "p4")

reduced_cor_df$Var1 <- factor(reduced_cor_df$Var1, 
                                    levels = ordered_covs, 
                                    labels = ordered_covs_labs)
reduced_cor_df$Var2 <- factor(reduced_cor_df$Var2, 
                            levels = c(paste0("PC", 4:1), "lat"))

if(!dir.exists("plots")) dir.create("plots")

# Figure 3
pdf(file.path("plots", "pc_cov_cor_reduced.pdf"), 18.4, 4)
ggplot(reduced_cor_df) + 
    geom_tile(aes(Var1, Var2, fill = value)) + 
    scale_fill_viridis_c(name = "Correlation", 
                        option = "cividis", limits = c(-1, 1)) + 
    theme(text = element_text(size = 18), 
          axis.text.x = element_text(angle = 35, hjust = 1, size = 16), 
          legend.title = element_text(size = 20, margin = margin(b = 15)), 
          legend.key.size = unit(1, "cm")) + 
    labs(x = "", y = "")
graphics.off()
