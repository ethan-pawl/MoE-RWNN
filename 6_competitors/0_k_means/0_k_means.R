library(dplyr)

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

seedfile <- file.path("6_competitors", "0_k_means", "seedtab.csv")
seedtab <- read.csv(seedfile)

kmeans_file <- file.path("6_competitors", "0_k_means", "best_kmeans__scores.Rdata")

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
    
    # Choose the initialization which minimizes the score
    if(cur_score < prev_score) {
      prev_score <- cur_score 
      best_model <- kmeans_out
    }
  }

  best_kmeans[[tt]] <- best_model 
  scores[tt] <- prev_score
}

# Use flowMatch helper functions to match clusters
for(tt in 2:TT) {
  print(tt)
  prev_samp <- flowMatch::ClusteredSample(labels = best_kmeans[[tt-1]]$cluster, sample = ylist[[tt-1]])
  cur_samp <- flowMatch::ClusteredSample(labels = best_kmeans[[tt]]$cluster, sample = ylist[[tt]])
  
  D <- flowMatch::dist.matrix(prev_samp, cur_samp, dist.type = "KL")
  rownames(D) <- colnames(D) <- as.character(1:K)

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

save(best_kmeans, scores, file = file.path("6_competitors", "0_k_means", "best_kmeans__scores-matched.Rdata"))