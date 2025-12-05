library(magrittr)
library(flowmix)
library(parallelly)

####################

# Settings

replic <- 0:1 # 0 is held-out dataset

nt <- 1000 # number of particles per cytogram (time index)
clust_sig <- 0.2 # cluster standard deviation
gridsize <- 50 # number of bins in each dimension
mc.cores <- max(availableCores() - 1, 1)

############## 

# Importing covariates

X_pc <- readRDS(file.path("data", "X_variations", "X_pc_9_nh_NA_seed_NA_ofold_NA_ifold_NA.RDS"))
# load(file.path("data", "X_pc.Rdata"))
# X_pc <- X

# Picking PC1 & PC2 to generate the means and PC1 & PC4 to generate the probs

X_mn <- X_pc[,c(1, 4)]
X_prob <- X_pc[,c(1, 2)]

X_mn_int <- cbind(1, X_mn)
X_prob_int <- cbind(1, X_prob)

###############################

# Specify regression coefficients, then 
# generate cluster means (mu) and probabilities (pi)

# The regression coefficients corresponding to PC1 are a 
# sign flip of the ones in the paper. We use PC1 to generate 
# the data and for modeling, but we report results 
# based on -PC1 for more intuitive interpretation.

# Cluster 1 (Prochlorococcus proxy)

# Linear in both
lin_pro_bet <- c(1, -0.015, 0.035)
lin_pro_mu <- X_mn_int %*% lin_pro_bet

softmax <- function(x) {
    (1 + exp(-x))**(-1)
}

lin_pro_alph <- c(0, 0.4, 0.1)
lin_pro_pi <- softmax(X_prob_int %*% lin_pro_alph)

# Interaction 
X_mn_inter <- cbind(1, X_mn, X_mn[,1] * X_mn[,2])
X_prob_inter <- cbind(1, X_prob, X_prob[,1] * X_prob[,2])

inter_pro_bet <- c(1, -0.015, 0.02, 0.003)
inter_pro_mu <- X_mn_inter %*% inter_pro_bet

inter_pro_alph <- c(0, 0.4, 0, 0.02)
inter_pro_pi <- softmax(X_prob_inter %*% inter_pro_alph)

# Quadratic in PC2
X_mn_quad <- cbind(1, X_mn[,1], abs(X_mn[,2]) * X_mn[,2])
X_prob_quad <- cbind(1, X_prob[,1], X_prob[,2]**2)

quad_pro_bet <- c(1, -0.015, 0.01)
quad_pro_mu <- X_mn_quad %*% quad_pro_bet

quad_pro_alph <- c(0, 0.4, -0.05)
quad_pro_pi <- softmax(X_prob_quad %*% quad_pro_alph)

# Logistic function of PC1
X_mn_sig <- cbind(1, softmax(X_mn[,1]), X_mn[,2])
X_prob_sig <- cbind(1, softmax(X_prob[,1]), X_prob[,2])

sig_pro_bet <- c(1.25, -0.15, 0.035)
sig_pro_mu <- X_mn_sig %*% sig_pro_bet
 
sig_pro_alph <- c(-2, 4.5, 0.1)
sig_pro_pi <- softmax(X_prob_sig %*% sig_pro_alph)

#################################

# I want the clusters' highest points to be the same up to 3 decimal 
# places
pro_mu_list <- list(lin_pro_mu, inter_pro_mu, quad_pro_mu, sig_pro_mu)
pro_maxes <- sapply(pro_mu_list, max)

diffs <- sapply(pro_maxes, function(x) {
    x - pro_maxes[1] 
})

# Re-generate the means over time to have the same maxes

# Interaction 
inter_pro_bet[1] <- 1 - diffs[2]
inter_pro_mu <- X_mn_inter %*% inter_pro_bet

# Quadratic in PC2
quad_pro_bet[1] <- 1 - diffs[3]
quad_pro_mu <- X_mn_quad %*% quad_pro_bet

# Logistic function of PC1
sig_pro_bet[1] <- 1.25 - diffs[4]
sig_pro_mu <- X_mn_sig %*% sig_pro_bet

pro_mu_list <- list(lin_pro_mu, inter_pro_mu, quad_pro_mu, sig_pro_mu)
pro_maxes <- sapply(pro_mu_list, max)

# Now generate the cluster 2 (picoeukaryote proxy) means
# based on 10 different signal sizes
max_pico_mu <- pro_maxes + 4 * clust_sig
pro_mu_means <- sapply(pro_mu_list, mean)
min_pico_mu <- pro_mu_means

clust2_intercepts <- lapply(1:length(pro_mu_list), function(i) {
    seq(min_pico_mu[i], max_pico_mu[i], length.out = 10)
})

aux_simdata_dir <- file.path("1_simulation", "aux_simdata")
if(!dir.exists(aux_simdata_dir)) dir.create(aux_simdata_dir)

# Save for later (analysis of model results)
names(clust2_intercepts) <- c("linear", "inter_mu", "quad_mu", "sig_mu")
save(clust2_intercepts, file = file.path(aux_simdata_dir, "clust2_intercepts.Rdata"))

# Measure signal sizes and save for later
signals <- lapply(1:length(pro_mu_list), function(i) {
    clust2_intercepts[[i]] - pro_mu_means[i]
})

names(signals) <- c("linear", "inter_mu", "quad_mu", "sig_mu")
save(signals, file = file.path(aux_simdata_dir, "signals.Rdata"))

# Gather cluster 2 means for all scenarios
pico_mu_list <- lapply(clust2_intercepts, function(config) {
    lapply(config, function(int) {
        rep(int, nrow(X_pc))
    })
})

#################################

# Drawing particles

# RANDOM: must set.seed before using
make_zlist <- function(nt, pi_vec) {
    lapply(pi_vec, function(pi_kt) {
        sample(2, nt, replace = TRUE, prob = c(pi_kt, 1 - pi_kt))
    })
}

# RANDOM: must set.seed before using
make_ylist_from_zlist <- function(nt, zlist, mu1_vec, mu2_vec, clust_sig) {
    lapply(seq_along(mu1_vec), function(tt) {
        ifelse(zlist[[tt]] == 1, 
            rnorm(nt, mu1_vec[tt], clust_sig), 
            rnorm(nt, mu2_vec[tt], clust_sig)) %>% matrix()
    })
}

# RANDOM: must set.seed before using
make_ylist_and_zlist <- function(nt, pi_vec, mu1_vec, mu2_vec, clust_sig, ylist_names) {
    zlist <- make_zlist(nt, pi_vec)
    ylist <- make_ylist_from_zlist(nt, zlist, mu1_vec, mu2_vec, clust_sig)
    names(ylist) <- names(zlist) <- ylist_names
    list("ylist" = ylist, "zlist" = zlist)
}

################# 

# Data generation and binning

# for each of the 10 signal sizes, 
# generate a train and a test (held-out) dataset
# from each of the 7 data data configurations

pro_mu_list <- list(lin_pro_mu, inter_pro_mu, quad_pro_mu, sig_pro_mu)
pi_list <- list(lin_pro_pi, inter_pro_pi, quad_pro_pi, sig_pro_pi)

# Matrix of data configurations
mu_pi_mat <- matrix(c(1, 1, 
                      2, 1, 
                      3, 1, 
                      4, 1, 
                      1, 2, 
                      1, 3, 
                      1, 4), ncol = 2, byrow = TRUE)
colnames(mu_pi_mat) <- c("mu", "pi")

X_mn_list <- list(X_mn_int, X_mn_inter, X_mn_quad, X_mn_sig)
X_prob_list <- list(X_prob_int, X_prob_inter, X_prob_quad, X_prob_quad, X_prob_sig)
bet_list <- list(lin_pro_bet, inter_pro_bet, quad_pro_bet, sig_pro_bet)
alph_list <- list(lin_pro_alph, inter_pro_alph, quad_pro_alph, sig_pro_alph)

make_data_fname <- function(imean, iprob, iint, irep) {
    paste0("simdata-", imean, "-", iprob, "-", iint, "-", irep, ".Rdata")
}

simdata_dir <- file.path("1_simulation", "simdata")
if(!dir.exists(simdata_dir)) dir.create(simdata_dir)

times <- rownames(X_pc)

# Create linear training data with the largest signal size
# and define a grid for binning based on that dataset

iint <- 10 # indexes the signal sizes
replic_seed <- 1 # training data
isettings <- 1 # indexes the data configuration

imean <- mu_pi_mat[1,"mu"]
iprob <- mu_pi_mat[1,"pi"]

set.seed(replic_seed)
sim_data <- make_ylist_and_zlist(nt, pi_list[[iprob]], 
                            pro_mu_list[[imean]], pico_mu_list[[imean]][[iint]], 
                            clust_sig, times)

manual.grid <- make_grid(sim_data$ylist, gridsize)
save(manual.grid, file = file.path(aux_simdata_dir, "manual_grid.Rdata"))

binobj <- bin_many_cytograms(sim_data$ylist, manual.grid, 
                                         mc.cores = mc.cores)

simdata1 <- list("ylist" = sim_data[[1]], "ybin_list" = binobj$ybin_list, 
                "countslist" = binobj$counts_list,
                "zlist" = sim_data[[2]], 
                "X_mu" = X_mn_list[[imean]], 
                "X_pi" = X_prob_list[[iprob]], 
                "bet" = bet_list[[imean]], 
                "alph" = alph_list[[iprob]],
                "gridsize" = gridsize,
                "mean_spec" = pro_mu_list[[imean]], 
                "prob_spec" = pi_list[[iprob]], "nt" = nt, 
                "pico_mu" = pico_mu_list[[imean]][[iint]], "clust_sig" = clust_sig, 
                "seed" = replic_seed, 
                "fname" = make_data_fname(imean, iprob, iint, replic_seed))

# Save the dataset
simdata <- simdata1
save(simdata, file = file.path(simdata_dir, make_data_fname(imean, iprob, iint, replic_seed)))

# Create the rest of the datasets, save them, and gather them into a list
dat_list <- lapply(1:10, function(iint) {
    lapply(replic, function(replic_seed) {
        lapply(1:nrow(mu_pi_mat), function(isettings) {
            if(iint == 10 & replic_seed == 1 & isettings == 1) {
                return(simdata1)
            } else {
                imean <- mu_pi_mat[isettings,"mu"]
                iprob <- mu_pi_mat[isettings,"pi"]

                fname <- make_data_fname(imean, iprob, iint, replic_seed)
                if(!file.exists(file.path(simdata_dir, fname))) {
                    print(paste0(fname, " doesn't exist, making now."))
                    set.seed(replic_seed)
                    sim_data <- make_ylist_and_zlist(nt, pi_list[[iprob]], 
                                    pro_mu_list[[imean]], pico_mu_list[[imean]][[iint]], 
                                    clust_sig, times)

                    binobj <- bin_many_cytograms(sim_data$ylist, manual.grid, 
                                                 mc.cores = mc.cores)

                    simdata <- list("ylist" = sim_data[[1]], "ybin_list" = binobj$ybin_list, 
                                    "countslist" = binobj$counts_list,
                                    "zlist" = sim_data[[2]], 
                                    "X_mu" = X_mn_list[[imean]], 
                                    "X_pi" = X_prob_list[[iprob]], 
                                    "bet" = bet_list[[imean]], 
                                    "alph" = alph_list[[iprob]],
                                    "gridsize" = gridsize,
                                    "mean_spec" = pro_mu_list[[imean]], 
                                    "prob_spec" = pi_list[[iprob]], "nt" = nt, 
                                    "pico_mu" = pico_mu_list[[imean]][[iint]], "clust_sig" = clust_sig, 
                                    "seed" = replic_seed, 
                                    "fname" = fname)

                    save(simdata, file = file.path(simdata_dir, fname))
                } else {
                    print(paste0(fname, " done."))
                    load(file.path(simdata_dir, fname))
                }
                return(simdata)
            }
        })
    }) %>% unlist(recursive = FALSE)
}) %>% unlist(recursive = FALSE)
