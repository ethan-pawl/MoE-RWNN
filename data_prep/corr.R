library(magrittr)
library(corrplot)
library(ggplot2)
library(reshape2)
library(dplyr)

datobj <- readRDS("/home/ethan/00_Cyto/data/paper-data-v2/MGL1704-hourly-paper.RDS")
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()

cormat <- cor(X)

png("corrplot0.png", width = 11.6, height = 8.5, units = "in", res = 300)
corrplot(cormat, method = "color", order = "hclust", tl.pos = "n")
graphics.off()

X_cont <- X[,3:39]
pr.out <- prcomp(X_cont, center = FALSE, scale = FALSE)
pr.var <- pr.out$sdev**2
pve <- pr.var / sum(pr.var)
cpve <- cumsum(pve)

plot(seq_along(cpve), cpve)

lat_X_cont <- cbind(lat, X_cont)

X_pc <- X_cont %*% pr.out$rotation

X_low_pc_1 <- X_pc[,1:9]

X_pc_cor <- cor(lat_X_cont, X_pc)
corrplot(X_pc_cor, method = "color")

X_low_pc <- X_cont %*% pr.out$rotation[,1:9]

X_low_pc_cor <- cor(lat_X_cont, cbind(lat, X_low_pc))
cor_df <- melt(X_low_pc_cor)
svg("cov_pc_cor.svg", width = 6, height = 10)
ggplot(cor_df) + 
    geom_tile(aes(Var2, Var1, fill = value)) + 
    scale_fill_viridis_c(option = "cividis", limits = c(-1, 1))
graphics.off()

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

reduced_cor_df <- cor_df %>%
    filter(Var2 %in% c("lat", "PC1", "PC2", "PC3", "PC4"))

reduced_cor_df$Var1 <- factor(reduced_cor_df$Var1, 
                                    levels = ordered_covs, 
                                    labels = ordered_covs_labs)
reduced_cor_df$Var2 <- factor(reduced_cor_df$Var2, 
                            levels = c(paste0("PC", 4:1), "lat"))
pdf("pc_cov_cor_reduced.pdf", 18.4, 4)
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
# TODO: improve the above plot
# if extra time, try sorting by strength of association, but 
# don't worry too much because this will probably go in the appendix

pc1_lat_df <- cbind(X_low_pc[,1], lat)

plot(1:296, X_low_pc[,1])
plot(1:296, lat)

X_low_pc[c(21, 22, 120, 121),1]
lat[X_low_pc[,1] < 0.00000000001]
X_low_pc[,1]

lat[which.min(abs(X_low_pc[,1]))]# 34.104
# also just do a heatmap of the loadings directly
loadings_df <- melt(pr.out$rotation[,1:9])

svg("loadings.svg", width = 6, height = 10)
ggplot(loadings_df) + 
    geom_tile(aes(Var2, Var1, fill = value)) + 
    scale_fill_viridis_c(option = "cividis")
graphics.off()

scaled_loadings <- apply(pr.out$rotation[,1:9], 2, function(pc) {
    pc /   
})

# High absolute correlations. Could this lead to "label switching" 
# in the PDPs?

# "b1"----------------------------------X
# "b2"----------------------------------X
# "p1"
# "p2"
# "p3"
# "p4"
# "sss"---------------------------------X
# "sst"
# "Fe"
# "PP"----------------------------------X
# "Si"----------------------------------X
# "NO3"
# "CHL"---------------------------------X
# "PHYC"
# "PO4"---------------------------------X 
# "O2"
# "vgosa"-------------------------------X
# "vgos"
# "sla"
# "ugosa"-------------------------------X
# "ugos"
# "wind_stress"-------------------------X
# "eastward_wind"
# "surface_downward_eastward_stress"----X
# "wind_speed"
# "surface_downward_northward_stress"---X
# "northward_wind"
# "ftle_bw_sla"
# "disp_bw_sla"
# "AOU_WOA_clim"
# "density_WOA_clim"
# "o2sat_WOA_clim"----------------------X
# "oxygen_WOA_clim"---------------------X
# "salinity_WOA_clim"
# "conductivity_WOA_clim"---------------X
# "nitrate_WOA_clim"
# "phosphate_WOA_clim"------------------X
# "silicate_WOA_clim"
# "par"

# try taking 
reduced <- c("p1",
            "p2",
            "p3",
            "p4",
            "sst",
            "vgos",
            "disp_bw_sla",
            "AOU_WOA_clim",
            "Fe",
            "ftle_bw_sla",
            "northward_wind",
            "nitrate_WOA_clim",
            "sla",
            "wind_speed",
            "PHYC",
            "silicate_WOA_clim",
            "NO3",
            "eastward_wind",
            "density_WOA_clim",
            "O2",
            "par",
            "ugos")

X_reduced <- X[,reduced]
corrplot(cor(X_reduced), method = "color", order = "hclust")

reduced_2 <- c("p1",
            "p2",
            "p3",
            "p4",
            "sst",
            "vgos",
            "disp_bw_sla",
            "AOU_WOA_clim",
            "Fe",
            "ftle_bw_sla",
            "northward_wind",
            "nitrate_WOA_clim",
            "sla",
            "wind_speed",
            # "PHYC",
            # "silicate_WOA_clim",
            # "NO3",
            # "eastward_wind",
            # "density_WOA_clim",
            "O2",
            "par",
            "ugos")


X_reduced_2 <- X[,reduced_2]
corrplot(cor(X_reduced_2), method = "color", order = "hclust")

reduced_3 <- c("p1",
            "p2",
            "p3",
            "p4",
            "sst",
            "vgos",
            # "disp_bw_sla",
            # "AOU_WOA_clim",
            "Fe",
            "ftle_bw_sla",
            "northward_wind",
            # "nitrate_WOA_clim",
            "sla",
            "wind_speed",
            # "PHYC",
            # "silicate_WOA_clim",
            # "NO3",
            # "eastward_wind",
            # "density_WOA_clim",
            # "O2",
            "par",
            "ugos")


X_reduced_3 <- X[,reduced_3]
corrplot(cor(X_reduced_3), method = "color", order = "hclust")

reduced_4 <- c("p1",
            "p2",
            "p3",
            "p4",
            "sst",
            # "vgos",
            # "disp_bw_sla",
            # "AOU_WOA_clim",
            "Fe",
            # "ftle_bw_sla",
            "northward_wind",
            # "nitrate_WOA_clim",
            "sla",
            "wind_speed",
            # "PHYC",
            # "silicate_WOA_clim",
            # "NO3",
            # "eastward_wind",
            # "density_WOA_clim",
            # "O2",
            "par",
            "ugos")


X_reduced_4 <- X[,reduced_4]
corrplot(cor(X_reduced_4), method = "color", order = "hclust")

setdiff(colnames(X), reduced_4)

which(cor(X_reduced_4) > 0.6, arr.ind = TRUE)
# only diags

##### X - PC correlation heatmap redo

which(colnames(lat_X_cont) == "silicate_WOA_clim")
colnames(lat_X_cont)[22:37] <- c("e_wind", "sdes", "wind_speed", "sdns", "n_wind", "ftle", "disp", "AOU", "density", "o2sat", "oxygen", 
                                "salinity", "conductivity", "nitrate", "phosphate", "silicate")
X_low_pc_cor <- cor(lat_X_cont, cbind(lat, X_low_pc))
cor_df <- melt(X_low_pc_cor)
svg("cov_pc_cor.svg", width = 6, height = 10)
ggplot(cor_df) + 
    geom_tile(aes(Var2, Var1, fill = value)) + 
    scale_fill_viridis_c(option = "cividis", limits = c(-1, 1)) + 
    xlab("") + ylab("") + labs(fill = "Corr") + 
    guides(x = guide_axis(n.dodge = 2))
graphics.off()
