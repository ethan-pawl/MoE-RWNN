library(flowmix)
library(parallel)

n_PCs <- as.character(c(9, 18, 27, 37))
NNseed <- c(4, 1, 2, 1)
names(NNseed) <- n_PCs

datobj <- file.path("data", "MGL1704-hourly-paper.RDS") |> 
  readRDS()

ylist <- datobj$ylist 
countslist <- datobj$countslist

X_dir <- file.path("data", "X_variations")

# Don't refit 9 PCs; use what we already fit
reslist <- mclapply(n_PCs, function(n_PC) {
  if(n_PC == "9") {
    load(file.path("2_application", "1_all_data", "nl_fit.Rdata"))
    return(res)
  }

  X <- readRDS(file.path(X_dir, paste0("X_pc_", n_PC, "_nh_70_seed_", NNseed[n_PC], "_ofold_NA_ifold_NA.RDS")))
  load(file.path("2_application", "3_PCs", paste0(n_PC, "_settings.Rdata")))

  flowmix_once(
    ylist = ylist, 
    countslist = countslist, 
    X = X, 
    mean_lambda = mean_lambda, 
    prob_lambda = prob_lambda, 
    seed = seed, 
    maxdev = 0.5, 
    numclust = 10
  )

}, mc.cores = 4)

###

clust_inds <- list(
  "9" = c("Pro" = 9, "Syn" = 6, "Pico1" = 2, "Pico2" = 1), 
  "18" = c("Pro" = 9, "Syn" = 1, "Pico1" = 7, "Pico2" = 6), 
  "27" = c("Pro" = 10, "Syn" = 3, "Pico1" = 8, "Pico2" = 2), 
  "37" = c("Pro" = 9, "Syn" = 4, "Pico1" = 10, "Pico2" = 3)
)

plot_resps <- c("Diameter", "Phycoerythrin", "Relative Abundance")
plot_pops <- c("Prochlorococcus", "Synechococcus", "PicoEukaryote Cluster 2")

library(dplyr)

prc_df <- lapply(1:length(n_PCs), function(i) {
  n_PC <- n_PCs[i]
  cur_seed <- NNseed[n_PC]
  cur_res <- reslist[[i]]
  X <- readRDS(file.path(X_dir, paste0("X_pc_", n_PC, "_nh_NA_seed_NA_ofold_NA_ifold_NA.RDS")))
  
  PC1_seq <- seq(min(X[,1]), max(X[,1]), length.out = 30)
  PC1_ice_X <- tcrossprod(rep(1, 30), colMeans(X))
  PC1_ice_X[,1] <- PC1_seq 

  set.seed(cur_seed)
  PC1_ice_X_nl <- make_hidden_nodes(PC1_ice_X, 70, 0.5)

  PC1_preds <- predict(cur_res, newx = PC1_ice_X_nl, logits = FALSE)

  cur_prc_df <- data.frame(
    `PC1 (Proxy for Latitude)` = rep(PC1_seq, 3), 
    Value = c(
      PC1_preds$mn[,1,clust_inds[[i]]["Pro"]], 
      PC1_preds$mn[,3,clust_inds[[i]]["Syn"]], 
      PC1_preds$prob[,clust_inds[[i]]["Pico2"]]
    ), 
    Response = rep(plot_resps, each = 30), 
    Population = rep(plot_pops, each = 30),
    `Number of PCs` = rep(n_PC, 3 * 30),
    check.names = FALSE
  )

  return(cur_prc_df)
}) %>% bind_rows()

prc_df$`Number of PCs` <- factor(prc_df$`Number of PCs`, levels = n_PCs)

library(ggplot2)
library(gridExtra)
library(ggpubr)

plist <- lapply(1:length(plot_pops), function(j) {
  pop <- plot_pops[j]
  resp <- plot_resps[j]

  p <- ggplot(filter(prc_df, Population == pop, Response == resp)) + 
    geom_line(aes(`PC1 (Proxy for Latitude)`, Value, linetype = `Number of PCs`, color = `Number of PCs`)) + 
    scale_linetype_manual(
      values = c(
        "9" = "solid", 
        "18" = "dashed", 
        "27" = "dotdash", 
        "37" = "longdash"
      )
    ) + 
    labs(title = bquote(italic(.(pop))), y = resp, x = "") + 
    theme(
      plot.margin = margin(r = 35, b = 5),
      plot.title = element_text(size = 20), 
      axis.title = element_text(size = 18), 
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 18), 
      legend.title = element_text(size = 20),
      legend.key.size = unit(2.5, "lines"), 
    )
  
  if(j == 1) {
    leg <<- get_legend(p) |> 
      as_ggplot()
  }

  p <- p + theme(legend.position = "none")
})

plist[[length(plot_pops) + 1]] <- leg 
plist[[length(plot_pops) + 2]] <- text_grob(
  "PC1 (Proxy for Latitude)", size = 18, hjust = 0.375, vjust = -0.5
)

pdf(file.path("plots", "SuppFigure13.pdf"), 18.3, 5.6)
print({
  grid.arrange(
    grobs = plist, 
    layout_matrix = matrix(
      c(
        1, 2, 3, 4, 
        NA, 5, NA, NA  
      ), 
      nrow = 2, 
      byrow = TRUE
    ), 
    widths = c(1, 1, 1, 0.4), 
    heights = c(1, 0.05)
  )
})
graphics.off()
