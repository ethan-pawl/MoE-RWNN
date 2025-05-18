library(flowmix)
library(magrittr)
library(ggplot2)
library(parallel)
library(tibble)

sum_dir <- "cv_sums"
cv_sums <- list()

# Load CV Summaries from RDS files
k <- 1
for(i in c("half", "tenth", "linear")) { # First 5 are linear, last 5 are nonlinear
    for (j in 1:5) {
        cv_sums[[k]] <- readRDS(file.path(sum_dir, 
                                          paste0("oos_eval_v3_", i, 
                                          "_", j, ".RDS")))

        k <- k + 1
    }
}

# Load all data
datobj = readRDS(file = "~/01.Cyto/data/paper-data-v2/MGL1704-hourly-paper.RDS")
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()

# Get indices of out-of-sample data
load("/home/ethan/01.Cyto/scripts/2024_02_28_cv/out_sample_folds.RData")

#################

# ELM Transformation Functions

make_W_unif <- function(p, n.h, a) {
    return(matrix(runif((p + 1) * n.h, - a, a), nrow = p + 1))
}

activate <- function(V) {
    return(1 / (1 + exp(-V)))
}

X_hidden <- function(X, p, n.h, a) {
  W <- make_W_unif(p, n.h, a)
  return(activate(cbind(1, X) %*% W))
}

###################

# Re-create the ELM-transformed covariates
n.h <- 70
p <- ncol(X)

set.seed(0)
X_half <- X_hidden(X, p, n.h, 0.5)
set.seed(0)
X_tenth <- X_hidden(X, p, n.h, 0.1)
##############

# Check the hidden layer weights aren't too squished to 0 or 1

hist(as.vector(X_half), main = "Hidden Layer Weights: Unif(-0.5, 0.5)", xlab = "Weight")
hist(as.vector(X_tenth), main = "Hidden Layer Weights: Unif(-0.1, 0.1)", 
    xlab = "Weight", xlim = c(0, 1))
hist(as.vector(X), main = "Original (Centered and Scaled) Covariates", xlab = "Weight")

######################

# TODO:
# - for each model, get the model which minimized oos loglik and refit on all data. 
# Then,
#   - get in-sample loglikelihood
#   - plot clustering for each model
#   - plot partial responses

#######################

# Log-likelihoods

# subset out-of-sample data for each CV run
oos_y <- lapply(out_sample_folds, function(i) {
    ylist[i]
})

oos_counts <- lapply(out_sample_folds, function(i) {
    countslist[i]
})

# Linear X
oos_X <- lapply(out_sample_folds, function(i) {
    X[i,]
})

oos_X_half <- lapply(out_sample_folds, function(i) {
    X_half[i,]
})

oos_X_tenth <- lapply(out_sample_folds, function(i) {
    X_tenth[i,]
})

half_sums <- cv_sums[1:5]
tenth_sums <- cv_sums[6:10]
lin_sums <- cv_sums[11:15]

pred_lin <- list()
pred_half <- list()
pred_tenth <- list()
for(i in 1:5) {
    pred_lin[[i]] <- predict(lin_sums[[i]]$bestres, 
                             newx = oos_X[[i]])

    pred_half[[i]] <- predict(half_sums[[i]]$bestres, 
                            newx = oos_X_half[[i]])

    pred_tenth[[i]] <- predict(tenth_sums[[i]]$bestres, 
                            newx = oos_X_tenth[[i]])
}

oos_loglik_lin <- vector("numeric", 5)
oos_loglik_half <- vector("numeric", 5)
oos_loglik_tenth <- vector("numeric", 5)
for(i in 1:5) {
    oos_loglik_lin[i] <- objective_newdat(pred_lin[[i]], 
                                            oos_y[[i]],
                                            oos_counts[[i]])

    oos_loglik_half[i] <- objective_newdat(pred_half[[i]], 
                                            oos_y[[i]],
                                            oos_counts[[i]])

    oos_loglik_tenth[i] <- objective_newdat(pred_tenth[[i]], 
                                            oos_y[[i]],
                                            oos_counts[[i]])
}

in_loglik_lin <- vector("numeric", 5)
in_loglik_half <- vector("numeric", 5)
in_loglik_tenth <- vector("numeric", 5)
for(i in 1:5) {
    in_loglik_lin[i] <- objective_newdat(lin_sums[[i]]$bestres, 
                                            ylist,
                                            countslist)

    in_loglik_half[i] <- objective_newdat(half_sums[[i]]$bestres, 
                                            ylist,
                                            countslist)

    in_loglik_tenth[i] <- objective_newdat(tenth_sums[[i]]$bestres, 
                                            ylist,
                                            countslist)
}

in_loglik <- data.frame(linear = in_loglik_lin, 
                        half = in_loglik_half, 
                        tenth = in_loglik_tenth)

oos_loglik <- data.frame(linear = oos_loglik_lin, 
                        half = oos_loglik_half, 
                        tenth = oos_loglik_tenth)

apply(in_loglik, 2, mean)
apply(oos_loglik, 2, mean) # linear minimizes oos negative log likelihood
                                        # but the margin is very small

##############

# Checking chosen lambdas

lambda_heatmap <- function(res, title) {
    cvscore.mat <- res$cvscore.mat

    lam_list <- list()
    lam_list$prob_lambdas <- rownames(cvscore.mat)
    lam_list$mean_lambdas <- colnames(cvscore.mat)
    lam_grid <- expand.grid(lam_list)

    cvscore_df <- cbind(lam_grid, as.vector(cvscore.mat))
    colnames(cvscore_df) <- c("prob_lambda", "mean_lambda", "cvscore")

    ggplot(cvscore_df, aes(x = mean_lambda, y = prob_lambda, fill = cvscore)) + 
        geom_tile() + 
        scale_fill_viridis_c(option = "cividis") + 
        ggtitle(title)

}

dname <- "heatmaps"
if(!dir.exists(dname)) {
    dir.create(dname)

    for(i in c("linear", "half", "tenth")) { 
        for(j in 1:5) {
            prefix <- switch(i, 
                linear = "lin",
                half = "half", 
                tenth = "tenth")

            fname <- paste0("lam_tiles_", prefix, "_", j, ".png")
            fname <- file.path(dname, fname)

            res <- eval(str2lang(paste0(prefix, "_sums[[", j, "]]")))
            png(fname, width = 960, height = 960)
            print(lambda_heatmap(res, title = paste0(i, ", outer fold ", j)))
            graphics.off()
        }
    }
}

                                        # heatmaps pretty much look fine

out_sample_folds

flowtrend::plot_3d(ylist, lin_sums[[1]]$bestres, 21, countslist)
flowtrend::plot_3d(ylist, half_sums[[1]]$bestres, 21, countslist)
flowtrend::plot_3d(ylist, tenth_sums[[1]]$bestres, 21, countslist)
