library(magrittr)
library(flowmix)
library(parallelly)
library(mvtnorm)

####################

# Settings

replic <- 0:1 # 0 is held-out dataset

nt <- 1000 # number of particles per cytogram (time index)
clust_sig <- 0.2 # cluster standard deviation
gridsize <- 50 # number of bins in each dimension
mc.cores <- max(availableCores() - 1, 1)

clust1_corr <- matrix(c(1, 0.5, 0, 
                        0.5, 1, 0, 
                        0, 0, 1), 3)

clust2_corr <- matrix(c(1, -0.7, 0, 
                        -0.7, 1, 0.3, 
                        0, 0.3, 1), 3)

clust1_sds <- clust2_sds <- c(0.2, 0.15, 0.25)

clust1_cov <- diag(clust1_sds) %*% clust1_corr %*% diag(clust1_sds)
clust2_cov <- diag(clust2_sds) %*% clust2_corr %*% diag(clust2_sds)

#################################

# Load 1D simulation data

load(file.path("1_simulation", "aux_simdata", "2-1-5_means_and_pis.Rdata"))

inter_pro_mu <- cbind(inter_pro_mu, inter_pico_int, inter_pico_int)
inter_pico_mu <- cbind(rep(inter_pico_int, length(lin_pro_pi)), inter_pico_int, inter_pico_int)

################################

# Drawing particles

draw_particles <- function(nt, pi_vec, mu1_mat, mu2_mat, cov1, cov2) {
    TT <- length(pi_vec)
    U <- runif(nt * TT)
    
    # Draw cluster membership indicators
    Z <- U < rep(pi_vec, each = nt)
    Zmat <- matrix(Z, nt, TT)
    clust1_counts <- colSums(Zmat)

    # Draw from Gaussian based on its cluster membership.
    # Data are assumed independent so particles from the 
    # same cluster are drawn as one group.
    ylist <- lapply(1:TT, function(tt) {
        rbind(
            rmvnorm(clust1_counts[tt], mu1_mat[tt,], cov1), 
            rmvnorm(nt - clust1_counts[tt], mu2_mat[tt,], cov2)
        )
    }) 

    list(ylist = ylist, clust1_counts = clust1_counts)
}

# Name the files
make_data_fname <- function(imean, iprob, iint, dataSeed) {
    paste0(
        "simdata", 
        "_imean_", 
        imean, 
        "_iprob_", 
        iprob, 
        "_iint_", 
        iint,
        "_dataSeed_", 
        dataSeed, 
        ".RDS"
    )
}

simdata_dir <- file.path("4_3dsim", "3dsimdata")
if(!dir.exists(simdata_dir)) dir.create(simdata_dir)

# training data
set.seed(1)
sim_data1 <- draw_particles(
    nt, 
    lin_pro_pi, 
    inter_pro_mu, 
    inter_pico_mu, 
    clust1_cov, 
    clust2_cov
)

manual_grid_3d <- make_grid(sim_data1$ylist, gridsize)
save(manual_grid_3d, file = file.path(simdata_dir, "manual_grid_3d.Rdata"))

binobj1 <- bin_many_cytograms(
    sim_data1$ylist, 
    manual_grid_3d, 
    mc.cores = mc.cores
)

simdata1 <- list(
    "ylist" = sim_data1$ylist, 
    "ybin_list" = binobj1$ybin_list, 
    "countslist" = binobj1$counts_list,
    "clust1_counts" = sim_data1$clust1_counts, 
    "gridsize" = gridsize,
    "mean_spec" = inter_pro_mu, 
    "prob_spec" = lin_pro_pi,
    "nt" = nt, 
    "pico_mu" = inter_pico_mu, 
    "clust1_cov" = clust1_cov,
    "clust2_cov" = clust2_cov, 
    "seed" = 1, 
    "fname" = make_data_fname(
        imean = "2", 
        iprob = "1", 
        iint = "5", 
        dataSeed = "1"
    )
)

# Save the dataset
saveRDS(
    simdata1, 
    file = file.path(
        simdata_dir,  
        simdata1$fname
    )
)

# Test set
set.seed(0)
sim_data0 <- draw_particles(
    nt, 
    lin_pro_pi, 
    inter_pro_mu, 
    inter_pico_mu, 
    clust1_cov, 
    clust2_cov
)

binobj0 <- bin_many_cytograms(
    sim_data0$ylist, 
    manual_grid_3d, 
    mc.cores = mc.cores
)

simdata0 <- list(
    "ylist" = sim_data0$ylist, 
    "ybin_list" = binobj0$ybin_list, 
    "countslist" = binobj0$counts_list,
    "clust1_counts" = sim_data0$clust1_counts, 
    "gridsize" = gridsize,
    "mean_spec" = inter_pro_mu, 
    "prob_spec" = lin_pro_pi,
    "nt" = nt, 
    "pico_mu" = inter_pico_mu, 
    "clust1_cov" = clust1_cov,
    "clust2_cov" = clust2_cov, 
    "seed" = 0, 
    "fname" = make_data_fname(
        imean = "2", 
        iprob = "1", 
        iint = "5", 
        dataSeed = "0"
    )
)

# Save the dataset
saveRDS(
    simdata0, 
    file = file.path(
        simdata_dir,  
        simdata0$fname
    )
)