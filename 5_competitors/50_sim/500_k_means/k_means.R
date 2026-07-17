rerun_clustering <- FALSE
rerun_matching <- FALSE 
plot_animation <- TRUE

library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(gridExtra)
library(gganimate)
library(patchwork)

RNGkind("L'Ecuyer-CMRG")

nrep <- 10
K <- 2

# Load the data
# Load the data
datobj <- readRDS(
  file.path(
    "4_3dsim", 
    "3dsimdata", 
    "simdata_imean_2_iprob_1_iint_5_dataSeed_1.RDS"
  )
)

ylist <- datobj$ybin_list
countslist <- datobj$countslist

TT <- length(ylist)

seedfile <- file.path("5_competitors", "50_sim", "500_k_means", "seedtab.csv")
seedtab <- read.csv(seedfile)

kmeans_file <- file.path("5_competitors", "50_sim", "500_k_means", "best_kmeans__scores.Rdata")
if(!file.exists(kmeans_file) | rerun_clustering) {
  best_kmeans <- vector("list", TT)
  scores <- numeric(TT)

  # For each cytogram, use k-means clustering
  for(tt in 1:TT) {
    print(tt)
    # Select a single cytogram
    cur_y <- ylist[[tt]]

    # Filter the seedtab
    seedtab_tt <- seedtab %>% 
      filter(ifold == tt)
    
    prev_score <- Inf
    best_model <- NULL
    # scores <- numeric(nrep)
    # models <- vector("list", nrep)
    for(cur_irep in 1:nrep) {
      # Initialize the random seed to get several different (yet reproducible) mean initializations
      .Random.seed <- seedtab_tt %>% 
        filter(irep == cur_irep) %>% 
        select(starts_with("seed")) %>% 
        as.integer()

      kmeans_out <- kmeans(cur_y, K, iter.max = 15)

      # To score the different initializations, use within-cluster sum of squares / total sum of squares (lower is better)
      cur_score <- kmeans_out$tot.withinss / kmeans_out$totss
      # scores[cur_irep] <- cur_score 
      # models[[cur_irep]] <- kmeans_out
      
      # Choose the initialization which minimizes the score
      if(cur_score < prev_score) {
        prev_score <- cur_score 
        best_model <- kmeans_out
      }
    }

    best_kmeans[[tt]] <- best_model 
    scores[tt] <- prev_score
  }
  save(best_kmeans, scores, file = kmeans_file)
} else{
  load(kmeans_file)
}

match_file <- file.path("5_competitors", "50_sim", "500_k_means", "best_kmeans__scores-matched.Rdata")
if(!file.exists(match_file) | rerun_matching) {
  # Use flowMatch helper functions to match clusters
  for(tt in 2:TT) {
    print(tt)
    prev_samp <- flowMatch::ClusteredSample(labels = best_kmeans[[tt-1]]$cluster, sample = ylist[[tt-1]])
    cur_samp <- flowMatch::ClusteredSample(labels = best_kmeans[[tt]]$cluster, sample = ylist[[tt]])
    
    D <- flowMatch::dist.matrix(prev_samp, cur_samp, dist.type = "KL")
    rownames(D) <- colnames(D) <- as.character(1:K)

    # # Use the Hungarian algorithm to find the lowest total KL-divergence cluster match 
    # # (I don't use this because solution is pretty nonsensical)
    # solution <- RcppHungarian::HungarianSolver(D)
    # cur_samp_labels <- solution$pairs[,2][solution$pairs[,1]]

    # TODO: if time, try minimax instead of minimin
    cur_samp_labels <- integer(K)
    for(k in 1:K) {
      # Greedily match clusters based on minimal KL divergence
      min_ind <- which.min(D) |> arrayInd(dim(D))
      cur_samp_labels[as.integer(colnames(D)[min_ind[,2]])] <- as.integer(rownames(D)[min_ind[,1]])
      D <- D[-min_ind[,1],-min_ind[,2], drop = FALSE]
    }

    reordering <- order(cur_samp_labels)
    best_kmeans[[tt]]$size <- best_kmeans[[tt]]$size[reordering]
    best_kmeans[[tt]]$withinss <- best_kmeans[[tt]]$withinss[reordering]
    best_kmeans[[tt]]$centers <- best_kmeans[[tt]]$centers[reordering,]
    best_kmeans[[tt]]$cluster <- best_kmeans[[tt]]$cluster |> recode_values(from = 1:K, to = cur_samp_labels)
  }

  save(best_kmeans, scores, file = match_file)
} else {
  load(match_file)
}

means_probs_df <- lapply(1:TT, function(tt) {
  data.frame(
    mean_1 = best_kmeans[[tt]]$centers[,1],
    mean_2 = best_kmeans[[tt]]$centers[,2],
    mean_3 = best_kmeans[[tt]]$centers[,3],
    prob = best_kmeans[[tt]]$size / sum(best_kmeans[[tt]]$size),
    cluster = 1:K,
    time = tt
  )
}) |> bind_rows()

ylist_df <- lapply(1:TT, function(tt) {
  cur_ylist_df <- ylist[[tt]] |> as.data.frame()
  cur_ylist_df <- cbind(cur_ylist_df, countslist[[tt]])
  colnames(cur_ylist_df) <- c("y1", "y2", "y3", "Bin_Biomass")# , "Cluster Membership") # No ground truth here
  cur_ylist_df$time <- tt
  cur_ylist_df$cluster <- factor(best_kmeans[[tt]]$cluster)

  cur_ylist_df 
}) |> bind_rows()

if(plot_animation) {
  
  y1_breaks <- seq(0, 2.5, 0.5)
  y2_breaks <- seq(0.75, 2.25, 0.5)
  y3_breaks <- seq(0.25, 2.75, 0.5)

  set1 <- brewer.pal(3, "Set1")

  # TODO: maybe remove this for compatibility
  prob_range <- range(means_probs_df$prob, na.rm = TRUE)

  for(tt in 1:TT) {
    cur_ylist_df <- subset(ylist_df, time == tt)
    cur_mp_df <- subset(means_probs_df, time == tt)

    p1 <- ggplot(cur_ylist_df) +
      geom_tile(aes(y1, y2, alpha = Bin_Biomass, fill = cluster)) + 
      scale_x_continuous(breaks = y1_breaks) + 
      scale_y_continuous(breaks = y2_breaks) + 
      coord_cartesian(xlim = range(y1_breaks), ylim = range(y2_breaks)) + 
      scale_fill_manual(values = c(set1, "black")) + 
      scale_alpha_continuous(range = c(0.2, 1)) + 
      theme_gray(base_family = "sans") + 
      theme(axis.title = element_text(size = 8), plot.title = element_text(size = 10), 
        axis.text = element_text(size = 6), legend.position = "none", 
        axis.line = element_line(linewidth = 0.25)
      ) + 
      geom_point(
        data = cur_mp_df, 
        mapping = aes(x = mean_1, y = mean_2, size = prob)
      ) + 
      scale_size_continuous(range = c(0, 4)) + 
      geom_text_repel(data = cur_mp_df, 
        mapping = aes(mean_1, mean_2, label = cluster), size = 3, box.padding = 0.25, 
        min.segment.length = 0.25, segment.size = 0.25
      )
  
    p2 <- ggplot(cur_ylist_df) +
      geom_tile(aes(y2, y3, alpha = Bin_Biomass, fill = cluster)) + 
      scale_x_continuous(breaks = y2_breaks) + 
      scale_y_continuous(breaks = y3_breaks) + 
      coord_cartesian(xlim = range(y2_breaks), ylim = range(y3_breaks)) + 
      scale_fill_manual(values = c(set1, "black")) + 
      scale_alpha_continuous(range = c(0.2, 1)) + 
      theme_gray(base_family = "sans") + 
      theme(axis.title = element_text(size = 8), plot.title = element_text(size = 10), 
        axis.text = element_text(size = 6), legend.position = "none", 
        axis.line = element_line(linewidth = 0.25)
      ) + 
      geom_point(
        data = cur_mp_df, 
        mapping = aes(x = mean_2, y = mean_3, size = prob)
      ) +
      scale_size_continuous(range = c(0, 4)) + 
      geom_text_repel(data = cur_mp_df, 
        mapping = aes(mean_2, mean_3, label = cluster), size = 3, box.padding = 0.25, 
        min.segment.length = 0.25, segment.size = 0.25
      )

    png(
      sprintf(
        file.path(
          "5_competitors", "50_sim", "500_k_means", "plots", "frames", "frame_%04d.png"
        ), 
        tt
      ), width = 1080, height = 540, res = 80, type = "cairo")

    grid.arrange(
      p1, p2,
      ncol = 2,
      top = sprintf("Time %d", tt)
    )

    graphics.off()
  }

  # memory cache exhausted; try to rerun
  imgs <- magick::image_read(list.files(file.path("5_competitors", "50_sim", "500_k_means", "plots", "frames"), full.names = TRUE))
  gif <- magick::image_animate(imgs, fps = 10)
  magick::image_write(gif, "animation.gif")

}
