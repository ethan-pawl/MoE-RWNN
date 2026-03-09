library(flowtrend)
library(flowmix)
library(ggpubr)
library(gridExtra)

n_PCs <- c(9, 18, 27, 37)
n_h <- 70

datobj <- file.path("data", "MGL1704-hourly-paper.RDS") |> 
  readRDS()

ylist <- datobj$ylist 
countslist <- datobj$countslist

# overlay PRCs from Figure 6

load_X <- function(readdir, n_PCs = "NA", n_h = "NA", seed = "NA", ofold = "NA", ifold = "NA") {

  file_name <- paste("X_pc", n_PCs, "nh", n_h, "seed", seed, "ofold", ofold, "ifold", ifold, sep = "_")
  file_name <- paste0(file_name, ".RDS")
  file_name <- file.path(readdir, file_name)

  X <- readRDS(file_name)
  return(X)
}

X_dir <- file.path("data", "X_variations")

plots_dir <- "plots"
if(!dir.exists(plots_dir)) dir.create(plots_dir)

cvlist <- list()
reslist <- list()
cvscores <- matrix(NA, length(n_PCs), length(seeds))
for(i in 1:length(n_PCs)) {
  cvlist[[i]] <- list()
  reslist[[i]] <- list()
  for(cur_seed in 1:5) {

    cur_q <- n_PCs[i]

    if(!(cur_q == 27 & cur_seed == 3)) {
      cur_cv <- readRDS(
        file.path(
          "2_application", 
          "23_PCs", 
          "results", 
          paste0("nl_", cur_q, "_70_", cur_seed), 
          paste0("nl_", cur_q, "_70_", cur_seed, "_pcr_summary.RDS")
        )
      )

      cvlist[[i]][[cur_seed]] <- cur_cv
      cvscores[i,cur_seed] <- cur_cv$cvscore.mat[cur_cv$min.inds]

      cur_res <- cur_cv$bestres 
      reslist[[i]][[cur_seed]] <- cur_res

      cur_plots_dir <- file.path(plots_dir, paste0(cur_q, "_", cur_seed))
      if(!dir.exists(cur_plots_dir)) dir.create(cur_plots_dir)

      pdf(file.path(cur_plots_dir, paste0(cur_q, "_", cur_seed, "_clusters.pdf")), 15.75, 6)
      print(flowtrend::plot_3d(ylist, cur_res, 33, countslist, bin = TRUE))
      graphics.off()

    } else {
      cvlist[[i]][[cur_seed]] <- NA 
      reslist[[i]][[cur_seed]] <- NA 
    }

  }
}

#########

# Model matching

model_seeds <- expand.grid(
  `9` = c(4, 5),
  `18` = c(1, 2), 
  `27` = 2, 
  `37` = c(1, 3)
)

# Label clusters
clust_mats <- array(
  NA, c(length(n_PCs), length(seeds), 4), 
  dimnames = list(
    NULL,
    NULL,
    c("Pro", "Syn", "Pico1", "Pico2")
  )
)

clust_mats[1,4,] <- c(9, 6, 2, 1)
clust_mats[1,5,] <- c(5, 10, 1, 3)
clust_mats[2,1,] <- c(9, 1, 7, 6)
clust_mats[2,2,] <- c(8, 3, 7, 2)
clust_mats[3,2,] <- c(10, 3, 8, 2)
clust_mats[4,1,] <- c(9, 4, 10, 3)
clust_mats[4,3,] <- c(8, 5, 1, 7)

plot_resps <- c("Diameter", "Phycoerythrin", "Relative Abundance")
plot_pops <- c("Prochlorococcus", "Synechococcus", "PicoEukaryote Cluster 2")

df_list <- replicate(nrow(model_seeds), vector("list", ncol(model_seeds)), simplify = FALSE)
for(j in 1:ncol(model_seeds)) {
  cur_q <- n_PCs[j]
  cur_X <- load_X(X_dir, n_PCs = cur_q)
  PC1_seq <- seq(min(cur_X[,1]), max(cur_X[,1]), length.out = 30)
  PC1_ice_X <- tcrossprod(rep(1, 30), colMeans(cur_X))
  PC1_ice_X[,1] <- PC1_seq 

  for(i in 1:nrow(model_seeds)) {
    cur_seed <- model_seeds[i,j]
    set.seed(cur_seed)
    PC1_ice_X_nl <- make_hidden_nodes(PC1_ice_X, n_h, 0.5)

    cur_cv <- readRDS(
      file.path(
        "2_application", 
        "23_PCs", 
        "results", 
        paste0("nl_", cur_q, "_70_", cur_seed), 
        paste0("nl_", cur_q, "_70_", cur_seed, "_pcr_summary.RDS")
      )
    )

    cur_res <- cur_cv$bestres 
    PC1_preds <- predict(cur_res, newx = PC1_ice_X_nl, logits = FALSE)

    cur_prc_df <- data.frame(
      `PC1 (Proxy for Latitude)` = rep(PC1_seq, 3), 
      Value = c(
        PC1_preds$mn[,1,clust_mats[j,cur_seed,"Pro"]], 
        PC1_preds$mn[,3,clust_mats[j,cur_seed,"Syn"]], 
        PC1_preds$prob[,clust_mats[j,cur_seed,"Pico2"]]
      ), 
      Response = rep(plot_resps, each = 30), 
      Population = rep(plot_pops, each = 30),
      `Number of PCs` = rep(cur_q, 3 * 30),
      check.names = FALSE
    )

    df_list[[i]][[j]] <- cur_prc_df
  }
}

comb_plots_dir <- file.path(plots_dir, "0_seed_combs")
if(!dir.exists(comb_plots_dir)) dir.create(comb_plots_dir)

for(i in 1:nrow(model_seeds)) {
  prc_df <- bind_rows(df_list[[i]])
  prc_df$`Number of PCs` <- factor(prc_df$`Number of PCs`)

  fname <- paste(model_seeds[i,], collapse = "_") |> 
    paste0(".pdf")


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

  pdf(file.path(comb_plots_dir, fname), 18.3, 5.6)
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

}

final_cvscores <- c(cvscores[1,4], cvscores[2,1], cvscores[3,2], cvscores[4,1])
mean(final_cvscores)
sd(final_cvscores) / 2
