library(dplyr)

# Load the data
datobj <- readRDS(
  file.path(
    "4_3dsim", 
    "3dsimdata", 
    "simdata_imean_2_iprob_1_iint_5_dataSeed_1.RDS"
  )
)

ybin_list <- datobj$ybin_list
countslist <- datobj$countslist
rm(datobj)

real_dataobj <- readRDS(file = file.path("data", "MGL1704-hourly-paper.RDS"))
X <- real_dataobj$X
rm(real_dataobj)

# Load the clustering results
load(file.path("5_competitors", "50_sim", "500_k_means", "best_kmeans__scores-matched.Rdata"))

data_long <- lapply(seq_along(ybin_list), function(tt) {
  cur_y <- ybin_list[[tt]]
  colnames(cur_y) <- c("y1", "y2", "y3")
  nt <- nrow(cur_y)

  which_covs <- 3:ncol(X)
  xt_mat <- matrix(rep(X[tt,which_covs], each = nt), nt)
  colnames(xt_mat) <- colnames(X)[which_covs]

  Weight <- countslist[[tt]]
  Cluster <- best_kmeans[[tt]]$cluster
  res <- cbind(cur_y, xt_mat, Weight, Cluster) |> 
    as.data.frame()
  
  return(res)
}) %>% bind_rows()

# head(data_long)

data_long_by_clust_list <- data_long |> 
  group_by(Cluster) |> 
  group_split(.keep = FALSE)

# Last two columns are Weight and Cluster; we don't want these as covariates
rhs <- paste0("s(", colnames(data_long)[4:(ncol(data_long) - 2)], ", bs = \"cs\")") |> 
  paste(collapse = " + ")

y1_formula <- paste0("y1", " ~ ", rhs) |> 
  as.formula()

y2_formula <- paste0("y2", " ~ ", rhs) |> 
  as.formula()

y3_formula <- paste0("y3", " ~ ", rhs) |> 
  as.formula()

library(mgcv)

system.time({
  per_cluster_res <- lapply(data_long_by_clust_list, function(cur_data_long) {
    cur_gam <- gam(
    list(
      y1_formula, 
      y2_formula, 
      y3_formula
    ),
    data = cur_data_long,
    family = mvn(d = 3),
    select = TRUE, 
    gamma = 1
  )

    return(cur_gam)
  })
})

# TODO: skip this; always crashes
saveRDS(per_cluster_res, file.path("5_competitors", "50_sim", "501_gam", "mvgam_res.RDS"))

# TODO: run analysis script
# TODO: cross-validate over gamma