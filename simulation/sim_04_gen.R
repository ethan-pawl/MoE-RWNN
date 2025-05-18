library(magrittr)
library(flowmix)
library(parallelly)
library(ggplot2)
library(dplyr)
library(tibble)

####################

# Settings

replic <- 0:1 # 0 is held-out

nt <- 1000
clust_sig <- 0.2
gridsize <- 50
mc.cores <- availableCores() - 1

############## 

# Importing covariates

load(file.path("~", 
               "00_Cyto", 
               "data", 
               "X_data", 
               "X_pc.Rdata"))

X_pc <- X

load(file.path("~", 
               "00_Cyto", 
               "data", 
               "X_data", 
               "X_nl_70.Rdata"))

X_nl <- X

# picking PC1 & PC2 to generate the means and PC1 & PC4 to generate the probs

X_mn <- X_pc[,c(1, 4)]
X_prob <- X_pc[,c(1, 2)]

X_mn_int <- cbind(1, X_mn)
X_prob_int <- cbind(1, X_prob)

###############################

# Pro 

# Linear in both

plots_dir <- "plots"
if(!dir.exists(plots_dir)) dir.create(plots_dir)

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

# Quadratic in X2

# X_mn_quad <- cbind(1, X_mn[,1], X_mn[,2]**2)
X_mn_quad <- cbind(1, X_mn[,1], abs(X_mn[,2]) * X_mn[,2])
X_prob_quad <- cbind(1, X_prob[,1], X_prob[,2]**2)

# quad_pro_bet <- c(1, -0.015, 0.02)
# TODO: tune this and make new data
quad_pro_bet <- c(1, -0.015, 0.01)
quad_pro_mu <- X_mn_quad %*% quad_pro_bet

# Quad resulting in a larger jump
quad_pro_alph <- c(0, 0.4, -0.05)
quad_pro_pi <- softmax(X_prob_quad %*% quad_pro_alph)

# Sigmoid in X1
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

# Quadratic in X2

quad_pro_bet[1] <- 1 - diffs[3]
quad_pro_mu <- X_mn_quad %*% quad_pro_bet

# Sigmoid in X1

sig_pro_bet[1] <- 1.25 - diffs[4]
sig_pro_mu <- X_mn_sig %*% sig_pro_bet

# Check

pro_mu_list <- list(lin_pro_mu, inter_pro_mu, quad_pro_mu, sig_pro_mu)
pro_maxes <- sapply(pro_mu_list, max)
diffs <- sapply(pro_maxes, function(x) {
    x - pro_maxes[1]
})
# diffs are now 0

#############################

# Plot the mean and probability curves

pdf(file.path(plots_dir, "00_linear_mean.pdf"), 9, 8.2)
plot(lin_pro_mu, type = "l", ylim = c(0, 3))
graphics.off()

pdf(file.path(plots_dir, "00_linear_prob.pdf"), 9, 8.2)
plot(lin_pro_pi, type = "l", ylim = c(0, 1))
graphics.off()

pdf(file.path(plots_dir, "00_inter_mean.pdf"), 9, 8.2)
plot(inter_pro_mu, type = "l", ylim = c(0, 3))
graphics.off()

pdf(file.path(plots_dir, "00_inter_prob.pdf"), 9, 8.2)
plot(inter_pro_pi, type = "l", ylim = c(0, 1))
graphics.off()

pdf(file.path(plots_dir, "00_quad_mean.pdf"), 9, 8.2)
plot(quad_pro_mu, type = "l", ylim = c(0, 3))
graphics.off()

pdf(file.path(plots_dir, "00_quad_prob.pdf"), 9, 8.2)
plot(quad_pro_pi, type = "l", ylim = c(0, 1))
graphics.off()

pdf(file.path(plots_dir, "00_logistic_mean.pdf"), 9, 8.2)
plot(sig_pro_mu, type = "l", ylim = c(0, 3))
graphics.off()

pdf(file.path(plots_dir, "00_logistic_prob.pdf"), 9, 8.2)
plot(sig_pro_pi, type = "l", ylim = c(0, 1))
graphics.off()

###

pdf(file.path(plots_dir, "00_data.pdf"), 8.5, 11)
par(mfrow = c(4, 2))
plot(lin_pro_mu, type = "l", ylim = c(0, 3), main = "Linear Mean", ylab = "Prochlorococcus Mean")
plot(lin_pro_pi, type = "l", ylim = c(0, 1), main = "Linear Logit", ylab = "Prochlorococcus Probability")
plot(inter_pro_mu, type = "l", ylim = c(0, 3), main = "Interaction Mean", ylab = "Prochlorococcus Mean")
plot(inter_pro_pi, type = "l", ylim = c(0, 1), main = "Interaction Logit", ylab = "Prochlorococcus Probability")
plot(quad_pro_mu, type = "l", ylim = c(0, 3), main = "Quadratic Mean", ylab = "Prochlorococcus Mean")
plot(quad_pro_pi, type = "l", ylim = c(0, 1), main = "Quadratic Logit", ylab = "Prochlorococcus Probability")
plot(sig_pro_mu, type = "l", ylim = c(0, 3), main = "Logistic Mean", ylab = "Prochlorococcus Mean")
plot(sig_pro_pi, type = "l", ylim = c(0, 1), main = "Logistic Logit", ylab = "Prochlorococcus Probability")
graphics.off()

#########################

max_pico_mu <- pro_maxes + 4 * clust_sig
pro_mu_means <- sapply(pro_mu_list, mean)
min_pico_mu <- pro_mu_means

ints <- lapply(1:length(pro_mu_list), function(i) {
    seq(min_pico_mu[i], max_pico_mu[i], length.out = 10)
})

names(ints) <- c("linear", "inter_mu", "quad_mu", "sig_mu")
save(ints, file = "ints.Rdata")

signals <- lapply(1:length(pro_mu_list), function(i) {
    ints[[i]] - pro_mu_means[i]
})

names(signals) <- c("linear", "inter_mu", "quad_mu", "sig_mu")
save(signals, file = "signals.Rdata")


pico_mu_list <- lapply(ints, function(config) {
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

# for each of the 4 separations and 10 replications, 
#   generate data from each of the 7 configs

pro_mu_list <- list(lin_pro_mu, inter_pro_mu, quad_pro_mu, sig_pro_mu)
pi_list <- list(lin_pro_pi, inter_pro_pi, quad_pro_pi, sig_pro_pi)

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

dname <- file.path("~", 
                   "00_Cyto", 
                   "data", 
                   "simdata",
                   "simdata_04_01")

if(!dir.exists(dname)) {
    dir.create(dname)
}

# on a single computer
times <- rownames(X_pc)

# make manual grid
iint <- 10
replic_seed <- 1
isettings <- 1

imean <- mu_pi_mat[1,"mu"]
iprob <- mu_pi_mat[1,"pi"]

set.seed(replic_seed)
sim_data <- make_ylist_and_zlist(nt, pi_list[[iprob]], 
                            pro_mu_list[[imean]], pico_mu_list[[imean]][[iint]], 
                            clust_sig, times)

manual.grid <- make_grid(sim_data$ylist, gridsize)
save(manual.grid, file = "manual_grid.Rdata")

binobj <- bin_many_cytograms(sim_data$ylist, manual.grid, 
                                         mc.cores = mc.cores)

res1 <- list("ylist" = sim_data[[1]], "ybin_list" = binobj$ybin_list, 
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

res <- res1 
save(res, file = file.path(dname, make_data_fname(imean, iprob, iint, replic_seed)))

dat_list <- lapply(1:10, function(iint) {
    lapply(replic, function(replic_seed) {
        lapply(1:nrow(mu_pi_mat), function(isettings) {
            if(iint == 10 & replic_seed == 1 & isettings == 1) {
                return(res1)
            } else {
                imean <- mu_pi_mat[isettings,"mu"]
                iprob <- mu_pi_mat[isettings,"pi"]

                fname <- make_data_fname(imean, iprob, iint, replic_seed)
                if(!file.exists(file.path(dname, fname))) {
                    print(paste0(fname, " doesn't exist, making now."))
                    set.seed(replic_seed)
                    sim_data <- make_ylist_and_zlist(nt, pi_list[[iprob]], 
                                    pro_mu_list[[imean]], pico_mu_list[[imean]][[iint]], 
                                    clust_sig, times)

                    binobj <- bin_many_cytograms(sim_data$ylist, manual.grid, 
                                                 mc.cores = mc.cores)

                    res <- list("ylist" = sim_data[[1]], "ybin_list" = binobj$ybin_list, 
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

                    save(res, file = file.path(dname, fname))
                } else {
                    print(paste0(fname, " done."))
                    load(file.path(dname, fname))
                }
                return(res)
            }
        })
    }) %>% unlist(recursive = FALSE)
}) %>% unlist(recursive = FALSE)

# plot the largest degree of separation for each scenario 
# for illustrative/reference purposes
plot_list <- lapply(dat_list, function(dat) {
        flowtrend::plot_1d(ylist = dat$ybin_list, countslist = dat$countslist, bin = TRUE) + 
            geom_line(aes(time, mn, group = cluster, linewidth = prob, col = cluster), 
                        data = data.frame(time = rep(1:296, 2), cluster = factor(rep(1:2, each = 296)), 
                                            mn = c(dat$mean_spec, dat$pico_mu), 
                                            prob = c(dat$prob_spec, 1 - dat$prob_spec)), 
                        lineend = "round")
})

# plot 1-7, 

plot_list_no_means <- lapply(dat_list, function(dat) {
        flowtrend::plot_1d(ylist = dat$ybin_list, countslist = dat$countslist, bin = TRUE)
})

fnames <- c("1-1-10-1", "2-1-10-1", "3-1-10-1", "4-1-10-1", "1-2-10-1", 
            "1-3-10-1", "1-4-10-1")

fnames_2 <- paste0(fnames, "_no_means")
fnames <- file.path(plots_dir, paste0(fnames, ".pdf"))
fnames_2 <- file.path(plots_dir, paste0(fnames_2, ".pdf"))

for(i in 1:7) {
    pdf(fnames[i], 10, 7.4)
    print(plot_list[[133 + i]])
    graphics.off()

    pdf(fnames_2[i], 10, 7.4)
    print(plot_list_no_means[[133 + i]])
    graphics.off()
}

graphics.off()

