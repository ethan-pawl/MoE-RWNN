library(flowmix)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(parallel)
library(tibble)
library(reshape2)
library(gridExtra)
library(gtools)
library(grid)
library(dplyr)
library(tidyr)

sum_dir <- "cv_sums"
cv_sums <- list()

models <- c("half", "linear")
for(i in models) {
  cv_sums[[i]] <- readRDS(file.path(sum_dir, 
                                    paste0("pcr_summary_", 
                                           i, 
                                           ".RDS")))
}

linear_cv <- cv_sums$linear 
half_cv <- cv_sums$half 

min(linear_cv$cvscore.mat) # 3.328
min(half_cv$cvscore.mat) # 3.322

rm(models)

linear_best <- linear_cv$bestres
half_best <- half_cv$bestres

# Load all data
datobj = readRDS(file = "~/00_Cyto/data/paper-data-v2/MGL1704-hourly-paper.RDS")
datobj %>% list2env(envir = .GlobalEnv) %>% invisible()
# rm(datobj)


###################################

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

# Re-create the Rotated and ELM-transformed covariates
n.h <- 70
p <- ncol(X)

#  Project X into the principal components space
X_no_cp <- X[,3:ncol(X)]

# X has already been centered and scaled
X_pca <- prcomp(X_no_cp, center = FALSE, scale = FALSE)
pca_var <- X_pca$sdev**2
cpve <- cumsum(pca_var / sum(pca_var)) 

d <- which(cpve > 0.95)[1]
# The first d principal components explain most of the variation

X_pc <- X_no_cp %*% X_pca$rotation[,1:d]

set.seed(0)
X_nl <- X_hidden(X_pc, 9, 70, 0.5)

# set.seed(0)
# predictions <- predict(half_best, newx = X_hidden(X_pc, 9, n.h, 0.5))
# str(predictions, max.level = 1)
# 
# plot(1:296, predictions$mn[,1,8])
# # TODO: hold sunlight constant and view predictions of diameter 
# # over time
# 
# lin_predictions <- predict(linear_best, newx = X_pc)
# plot(1:296, lin_predictions$mn[,1,10])
# 
# plot(lat, lin_predictions$mn[,1,10], type = "l")

##########################

# Partial Dependence Grid

my_kde2d <- function (x, y, h, gx, gy, lims = c(range(x), range(y))) {
    nx <- length(x)
    if (length(y) != nx) 
        stop("data vectors must be the same length")
    if (any(!is.finite(x)) || any(!is.finite(y))) 
        stop("missing or infinite values in the data are not allowed")
    if (any(!is.finite(lims))) 
        stop("only finite values are allowed in 'lims'")
    # n <- rep(n, length.out = 2L)
    # gx <- seq.int(lims[1L], lims[2L], length.out = n[1L])
    # gy <- seq.int(lims[3L], lims[4L], length.out = n[2L])
    h <- if (missing(h)) 
        c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))
    else rep(h, length.out = 2L)
    if (any(h <= 0)) 
        stop("bandwidths must be strictly positive")
    h <- h/4
    ax <- outer(gx, x, "-")/h[1L]
    ay <- outer(gy, y, "-")/h[2L]
    z <- tcrossprod(matrix(dnorm(ax), , nx), matrix(dnorm(ay), 
        , nx))/(nx * h[1L] * h[2L])
    list(x = gx, y = gy, z = z)
}

# Helper function to plot exactly one 1D partial dependence plot and ICE curves.
# Requires the ICE and PD curves to be previously calculated
plot_pdp_1d <- function(ice_pd, vary, resp, lim) {
  ggplot(ice_pd) + 
    geom_smooth(data = ice_pd[ice_pd$ice_num != "0",], 
      aes(.data[[vary]], .data[[resp]], group = ice_num), 
      color = "gray", method = "loess", formula = y ~ x, 
      se = FALSE, linewidth = 0.1) + 
    geom_smooth(data = ice_pd[ice_pd$ice_num == "0",], 
      aes(.data[[vary]], .data[[resp]]), 
      color = "black", method = "loess", formula = y ~ x, 
      se = FALSE, linewidth = 0.1) + 
    scale_y_continuous(limits = lim, n.breaks = 5) + 
    theme_bw() + 
    theme(legend.position = "none")
}

# Helper function to plot exactly one 2D partial dependence plot.
# Requires the PD contour map to be previously calculated.
plot_pdp_2d <- function(pd_df, vary1, vary2, resp, lims, alpha_vary) {
  if(alpha_vary) {
    ggplot(pd_df, aes(.data[[vary1]], .data[[vary2]])) + 
      geom_raster(aes(fill = .data[[resp]], alpha = kde), interpolate = TRUE) + 
      # geom_contour(aes(z = .data[[resp]]), breaks = seq(lims[1], lims[2], by = breakBy), 
      #   color = "red", linetype = "dotted", alpha = 0.85) + 
      scale_fill_viridis_c(option = "turbo", limits = lims) +
      theme(legend.position = "none")
  } else {
    ggplot(pd_df, aes(.data[[vary1]], .data[[vary2]])) + 
      geom_raster(aes(fill = .data[[resp]]), interpolate = TRUE) + 
      # geom_contour(aes(z = .data[[resp]]), breaks = seq(lims[1], lims[2], by = breakBy), 
      #   color = "red", linetype = "dotted", alpha = 0.85) + 
      scale_fill_viridis_c(option = "turbo", limits = lims) + 
      theme(legend.position = "none")
  }
  
}

pdp_make_pred <- function(X, fit, gridsize = 30, nIce = 30, n.h = NULL, logits = TRUE) {
  uni <- colnames(X)
  n <- nrow(X)
  p <- ncol(X)

  ice_rownums <- sample.int(n, size = nIce)
  ice_X <- X[ice_rownums,]

  ice_grid <- diag(nIce) %x% rep(1, gridsize) %*% as.matrix(ice_X)

  # Univariate
  uni_res <- lapply(uni, function(vary) {
    vary_min <- min(X[,vary])
    vary_max <- max(X[,vary])

    ice_grid_vary <- seq(vary_min, vary_max, length.out = gridsize)
    ice_grid[,vary] <- rep(ice_grid_vary, times = nIce)

    if(!is.null(n.h)) {
      set.seed(0)
      ice_grid_half <- X_hidden(ice_grid, p, n.h, 0.5)
      ice_pred <- predict(fit, logits = logits, newx = ice_grid_half)
    } else {
      ice_pred <- predict(fit, logits = logits, newx = ice_grid)
    }

    return(list(ice_grid = ice_grid[,vary,drop = FALSE], ice_pred = ice_pred))

  })

  # Bivariate
  ice_grid_2d <- diag(nIce) %x% rep(1, gridsize**2) %*% as.matrix(ice_X)
  bi_res <- lapply(combn(uni, 2, simplify = FALSE), function(cur_pair) {
    vary1 <- cur_pair[1]
    vary2 <- cur_pair[2]

    vary1_min <- min(X[,vary1])
    vary1_max <- max(X[,vary1])
    vary2_min <- min(X[,vary2])
    vary2_max <- max(X[,vary2])

    # TODO: pre-specify a grid for each covariate first to speed things up if 
    # necessary
    ice_grid_vary1 <- seq(vary1_min, vary1_max, length.out = gridsize)
    ice_grid_vary2 <- seq(vary2_min, vary2_max, length.out = gridsize)
    ice_grid_vary <- expand.grid(ice_grid_vary1, ice_grid_vary2)
    ice_grid_2d[,c(vary1, vary2)] <- rep(1, nIce) %x% as.matrix(ice_grid_vary)

    if(!is.null(n.h)) {
      set.seed(0)
      ice_grid_half <- X_hidden(ice_grid_2d, p, n.h, 0.5)
      ice_pred_2d <- predict(fit, logits = logits, newx = ice_grid_half)
    } else {
      ice_pred_2d <- predict(fit, logits = logits, newx = ice_grid_2d)
    }

    kde <- my_kde2d(X[,vary1], X[,vary2], gx = ice_grid_vary1, gy = ice_grid_vary2)
    kde <- as.vector(kde$z)

    return(list(ice_grid = ice_grid_2d[,c(vary1, vary2)], 
                ice_pred = ice_pred_2d, 
                kde = kde,
                cur_pair = cur_pair))
  })

  # Predictions for training data
  if(!is.null(n.h)) {
      set.seed(0)
      X_train <- X_hidden(X, p, n.h, 0.5)
      all_pred <- predict(fit, logits = logits, newx = X_train)
  } else {
      all_pred <- predict(fit, logits = logits, newx = X)
  }
  
  return(list(uni_res = uni_res, bi_res = bi_res, ice_rownums = ice_rownums, 
              all_pred = all_pred, orig_X = X))
}

plot_pdps <- function(grid_pred, clust, resp, lim, gridsize, alpha_vary) {
  uni_plot_list <- lapply(grid_pred$uni_res, function(res) {
    resp_data <- switch(resp,
      logitProb = res$ice_pred$prob[,clust],
      diam = res$ice_pred$mn[,1,clust],
      chl = res$ice_pred$mn[,2,clust],
      pe = res$ice_pred$mn[,3,clust]
    )
    ice_df <- cbind(res$ice_grid, resp_data)
    colnames(ice_df)[2] <- resp
    ice_df <- as.data.frame(ice_df)
    ice_df$ice_num <- as.character(rep(grid_pred$ice_rownums, each = gridsize))

    pd_df <- ice_df[,2] %>% 
      matrix(ncol = length(grid_pred$ice_rownums)) %>%
      rowMeans() %>% 
      cbind(ice_df[seq_len(gridsize),1], .) %>% # this is ice_grid_vary (CHECK)
      as.data.frame()

    pd_df$ice_num <- "0"
    colnames(pd_df) <- colnames(ice_df) # TODO: remove this if unnecessary
    ice_pd <- rbind(ice_df, pd_df)

    plot_pdp_1d(ice_pd, colnames(pd_df)[1], resp, lim)
  })

  bi_plot_list <- lapply(grid_pred$bi_res, function(res) {
    resp_data <- switch(resp,
      logitProb = res$ice_pred$prob[,clust],
      diam = res$ice_pred$mn[,1,clust],
      chl = res$ice_pred$mn[,2,clust],
      pe = res$ice_pred$mn[,3,clust])

    ice_df <- cbind(res$ice_grid, resp_data)
    colnames(ice_df)[3] <- resp
    ice_df <- as.data.frame(ice_df)

    pd_df <- ice_df[,3] %>% 
    matrix(ncol = length(grid_pred$ice_rownums)) %>%
    rowMeans() %>% 
    cbind(ice_df[seq_len(gridsize**2),c(1, 2)], .) %>% # this is ice_grid_vary (CHECK)
    as.data.frame()

    colnames(pd_df) <- colnames(ice_df)

    pd_df$kde <- res$kde
    plot_pdp_2d(pd_df, colnames(pd_df)[1], colnames(pd_df)[2], resp, lim, alpha_vary)
  })

  scatterplot_resp <- switch(resp,
    logitProb = grid_pred$all_pred$prob[,clust],
    diam = grid_pred$all_pred$mn[,1,clust],
    chl = grid_pred$all_pred$mn[,2,clust],
    pe = grid_pred$all_pred$mn[,3,clust]
  ) %>% as.data.frame()

  colnames(scatterplot_resp)[1] <- resp

  scatterplot_data <- cbind(grid_pred$orig_X, scatterplot_resp)

  scatterplot_list <- list()
  scatterplot_list[[1]] <- ggplot(scatterplot_data, 
        aes(.data[[colnames(grid_pred$bi_res[[1]]$ice_grid)[1]]], 
            .data[[colnames(grid_pred$bi_res[[1]]$ice_grid)[2]]], 
            colour = .data[[resp]])) + 
    geom_point(shape = ".") + 
    scale_color_viridis_c(option = "turbo", lim = lim, n.breaks = 5)
  
  leg <- ggpubr::get_legend(scatterplot_list[[1]])
  leg <- ggpubr::as_ggplot(leg)

  scatterplot_list[[1]] <- scatterplot_list[[1]] + theme(legend.position = "none")

  scatterplot_list <- lapply(grid_pred$bi_res, function(res) {
    ggplot(scatterplot_data, 
          aes(.data[[res$cur_pair[1]]], 
              .data[[res$cur_pair[2]]], 
              colour = .data[[resp]])) + 
      geom_point(shape = ".") + 
      scale_color_viridis_c(option = "turbo", lim = lim) + 
      theme(legend.position = "none")
  })

  return(list(uni_plot_list = uni_plot_list, bi_plot_list = bi_plot_list, 
              scatterplot_list = scatterplot_list, leg = list(leg)))
}

my_pdpPairs <- function(plot_list) {
  all_plots <- c(unlist(plot_list, recursive = FALSE))

  n_uni_plots <- length(plot_list$uni_plot_list)
  n_bi_plots <- length(plot_list$bi_plot_list)

  layout_mat <- diag(seq_len(n_uni_plots))
  layout_mat[lower.tri(layout_mat)] <- seq_len(n_bi_plots) + n_uni_plots
  layout_mat <- t(layout_mat)
  layout_mat[lower.tri(layout_mat)] <- seq_along(plot_list$scatterplot_list) + n_uni_plots + n_bi_plots

  layout_mat <- cbind(layout_mat, rep(length(all_plots), nrow(layout_mat)))

  gridExtra::grid.arrange(grobs = all_plots, layout_matrix = layout_mat)
}

linear_clust <- c(10, 1, 9, 5, 7)
half_clust <- c(8, 3, 10, 7, 5)
clust_names <- c("pro2", "syn", "pro1", "pico1", "pico2")
resps <- c("logitProb", "diam", "chl", "pe")

set.seed(0)
lin_pred <- pdp_make_pred(X_pc, linear_best)
set.seed(0)
half_pred <- pdp_make_pred(X_pc, half_best, n.h = 70)

dname <- "PC_partial_responses"
if(!dir.exists(dname)) {
  dir.create(dname)
}

list_pop_limits <- function(chl, diam, logitProb, pe) {
  res <- list(chl = chl, diam = diam, logitProb = logitProb, 
            pe = pe)

  return(res) 
}

pop_limits <- list()
pop_limits$pro2 <- list_pop_limits(c(0.75, 1.75),
                                    c(0.5, 2), 
                                    c(-6, 6), 
                                    c(0.25, 0.5))
pop_limits$syn <- list_pop_limits(c(2, 3), 
                                    c(2.5, 3.5), 
                                    c(-3, 3), 
                                    c(2.5, 4))
pop_limits$pro1 <- list_pop_limits(c(0.5, 0.75),
                                    c(1.5, 3), 
                                    c(0, 3), 
                                    c(0.25, 0.5))
pop_limits$pico1 <- list_pop_limits(c(5.25, 6.25),
                                    c(4.5, 6.5), 
                                    c(-4, 4), 
                                    c(0.5, 1))
pop_limits$pico2 <- list_pop_limits(c(6, 8),
                                    c(6.25, 7.25), 
                                    c(0, 1.5), 
                                    c(1, 2))

k <- 1
for(i in seq_along(resps)) {
  for(j in seq_along(linear_clust)) {
    fname <- paste0(dname, "/", clust_names[j], "_", resps[i], 
                  "_linear.svg")
    CairoSVG(fname, width = 11, height = 9)
    plots <- plot_pdps(lin_pred, linear_clust[j], resps[i], 
                      pop_limits[[j]][[resps[i]]], 30, alpha_vary = TRUE)
    my_pdpPairs(plots)
    graphics.off()
    print(paste0(round(k / 40, 2), "; ", Sys.time()))
    k <- k + 1
  }

  for(j in seq_along(half_clust)) {
    fname <- paste0(dname, "/", clust_names[j], "_", resps[i], 
                  "_half.svg")
    CairoSVG(fname, width = 11, height = 9)
    plots <- plot_pdps(half_pred, half_clust[j], resps[i], 
                      pop_limits[[j]][[resps[i]]], 30, alpha_vary = TRUE)
    my_pdpPairs(plots)
    graphics.off()
    print(paste0(round(k / 40, 2), "; ", Sys.time()))
    k <- k + 1
  }
}

str(plots, max.level = 1)
print(plots$bi_plot_list[[1]])
plots
plots$uni_plot_list[[3]]

## W/o transparency

linear_clust <- c(10, 1, 5, 7)
half_clust <- c(8, 3, 7, 5)
clust_names <- c("pro", "syn", "pico1", "pico2")
resps <- c("logitProb", "diam", "chl", "pe")

set.seed(0)
lin_pred <- pdp_make_pred(X_pc, linear_best)
set.seed(0)
half_pred <- pdp_make_pred(X_pc, half_best, n.h = 70)

pop_limits <- list()
pop_limits$pro <- list_pop_limits(c(0.75, 1.75),
                                    c(0.5, 2), 
                                    c(-6, 6), 
                                    c(0.25, 0.5))
pop_limits$syn <- list_pop_limits(c(2, 3), 
                                    c(2.5, 3.5), 
                                    c(-3, 3), 
                                    c(2.5, 4))
pop_limits$pico1 <- list_pop_limits(c(5.25, 6.25),
                                    c(4.5, 6.5), 
                                    c(-4, 4), 
                                    c(0.5, 1))
pop_limits$pico2 <- list_pop_limits(c(6, 8),
                                    c(6.25, 7.25), 
                                    c(0, 1.5), 
                                    c(1, 2))

dname <- "PC_PDP_no_alpha"
k <- 1
for(i in seq_along(resps)) {
  for(j in seq_along(linear_clust)) {
    fname <- paste0(dname, "/", clust_names[j], "_", resps[i], 
                  "_linear.svg")
    CairoSVG(fname, width = 11, height = 9)
    plots <- plot_pdps(lin_pred, linear_clust[j], resps[i], 
                      pop_limits[[j]][[resps[i]]], 30, alpha_vary = FALSE)
    my_pdpPairs(plots)
    graphics.off()
    print(paste0(round(k / 40, 2), "; ", Sys.time()))
    k <- k + 1
  }

  for(j in seq_along(half_clust)) {
    fname <- paste0(dname, "/", clust_names[j], "_", resps[i], 
                  "_half.svg")
    CairoSVG(fname, width = 11, height = 9)
    plots <- plot_pdps(half_pred, half_clust[j], resps[i], 
                      pop_limits[[j]][[resps[i]]], 30, alpha_vary = FALSE)
    my_pdpPairs(plots)
    graphics.off()
    print(paste0(round(k / 40, 2), "; ", Sys.time()))
    k <- k + 1
  }
}

#########################

X_pca$rotation[,1:3] %>% signif(3)

plot_ly(x = X_pc[,1], y = X_pc[,2], z = X_pc[,3], type = "scatter3d", mode = "markers", linetype = "solid")

pc_mat <- X_pca$rotation[,1:d]

# Interpreting the principal components
pc_list <- asplit(pc_mat, MARGIN = 2)
pc_each_sorted <- lapply(pc_list, function(pc) {
  sort(pc, decreasing = TRUE)  %>% 
    as.matrix() %>% 
    signif(2)
})

pc_str <- lapply(pc_each_sorted, function(pc) {
  plus_pos <- pc >= 0.1
  minus_pos <- pc <= -0.1

  plus_str <- paste0(pc[plus_pos], "(", rownames(pc)[plus_pos], ")", collapse = " + ")
  minus_str <- paste0(pc[minus_pos], "(", rownames(pc)[minus_pos], ")", collapse = " ")
  str2 <- paste0(plus_str, " + ... ", minus_str)

  str2
})

sink("pc_str.txt")
pc_str
sink()

dim(X)
length(lat)

cors_w_lat <- apply(X[,3:39], MARGIN = 2, function(x) {
  cor(x, lat)
}) %>% sort(decreasing = TRUE) %>% as.matrix()

sink("cors_w_lat.txt")
cors_w_lat
sink()

set.seed(0)
X_half <- X_hidden(X_pc, d, n.h, 0.5)

##############

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

    for(i in c("linear", "half")) { 
        fname <- paste0("lam_tiles_", i, ".png")
        fname <- file.path(dname, fname)

        res <- eval(str2lang(paste0(i, "_cv")))

        png(fname, width = 960, height = 960)
        print(lambda_heatmap(res, title = i))
        graphics.off()
    }
}

########################

# PAPER PLOT

# Plot clustering
plot_3d_model_comp <- function(tt) {
  plist <- flowtrend::plot_3d(ylist, linear_best, tt, countslist, 
                              return_list_of_plots = TRUE)
  plist_2 <- flowtrend::plot_3d(ylist, half_best, tt, countslist, 
                                return_list_of_plots = TRUE)

  xy_breaks <- seq(0, 8, 2)
  grid.arrange(grobs = list(plist[[1]] + 
                              scale_x_continuous(breaks = xy_breaks) + 
                              scale_y_continuous(breaks = xy_breaks) + 
                              labs(title = "Linear"), 
                            plist[[2]] + 
                              scale_x_continuous(breaks = xy_breaks) + 
                              scale_y_continuous(breaks = xy_breaks),
                            plist[[3]] + 
                              scale_x_continuous(breaks = xy_breaks) + 
                              scale_y_continuous(breaks = xy_breaks), 
                            plist_2[[1]] + 
                              scale_x_continuous(breaks = xy_breaks) + 
                              scale_y_continuous(breaks = xy_breaks) + 
                              labs(title = "Nonlinear"), 
                            plist_2[[2]] + 
                              scale_x_continuous(breaks = xy_breaks) + 
                              scale_y_continuous(breaks = xy_breaks),
                            plist_2[[3]] + 
                              scale_x_continuous(breaks = xy_breaks) + 
                              scale_y_continuous(breaks = xy_breaks), 
                            ggpubr::text_grob(paste0("Time = ", tt), size = 14)), 
                            layout_matrix = matrix(c(7, 7, 7, 1, 2, 3, 4, 5, 6), 
                                                   byrow = TRUE, nrow = 3), 
                            heights = c(0.1, 1, 1))
}

for(tt in c(6, 33, 92, 111, 255)) {
  pdf(paste0("paper_plots/clusters_", tt, ".pdf"), 12.75, 9.36)
  print(plot_3d_model_comp(tt))
  graphics.off()
}

# PAPER PLOT
plot_list <- flowtrend::plot_3d(ylist, linear_best, 33, countslist, 
                                return_list_of_plots = TRUE, labels = c("Syn", 
                                                                        "Other1", 
                                                                        "Bead", 
                                                                        "Other2", 
                                                                        "Pico1", 
                                                                        "Other3", 
                                                                        "Pico2", 
                                                                        "Other4", 
                                                                        "Other5", 
                                                                        "Pro"))
p1 <- plot_list[[1]]

plot_list_nl <- flowtrend::plot_3d(ylist, half_best, 33, countslist, 
                                return_list_of_plots = TRUE, labels = c("Other6", 
                                                                        "Other3", 
                                                                        "Syn", 
                                                                        "Other4", 
                                                                        "Pico2", 
                                                                        "Other1", 
                                                                        "Pico1", 
                                                                        "Pro", 
                                                                        "Other2", 
                                                                        "Other5"))
p2 <- plot_list_nl[[1]]

p1 <- p1 + 
  labs(title = "Linear Model Fit", y = "Chlorophyll", x = "") + 
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) + 
  scale_x_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) + 
  scale_y_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) 

p2 <- p2 + 
  labs(title = "Nonlinear Model Fit", y = "", x = "") + 
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) + 
  scale_x_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) + 
  scale_y_continuous(limits = c(-0.75, 9.25), breaks = seq(0, 8, 2)) 

x_lab <- text_grob("Log Diameter", size = 12, x = 0.515, y = 1)
overall_title <- text_grob("Model Fit on July 2nd, 2017, 8 AM", size = 16)

pdf("paper_plots/clusters_33_onerow.pdf", 10.5, 6)
grid.arrange(p1, p2, x_lab, overall_title, layout_matrix = matrix(c(4, 4, 
                                                                    1, 2, 
                                                                    3, 3), 
                                                                    byrow = TRUE, 
                                                                    nrow = 3), 
             heights = c(0.1, 1, 0.05), widths = c(1, 1, 0.05))
graphics.off()

trace(flowtrend::plot_3d, edit = TRUE)

# PAPER PLOT
lin_probs <- linear_best$prob[,c(10, 1, 5, 7)]
half_probs <- half_best$prob[,c(8, 3, 7, 5)]

time_dt <- as.POSIXct(time, format = "%Y-%m-%d-")
dt <- as.POSIXct(time, format = "%Y-%m-%dT%T")

probs <- data.frame(Relative_Abundance = c(lin_probs, half_probs), 
                    Cluster = rep(c("Pro", "Syn", "Pico1", "Pico2", 
                                  "Pro", "Syn", "Pico1", "Pico2"), each = 296), 
                    Model = rep(c("Linear", "Nonlinear"), each = 296 * 4), 
                    Time = rep(dt, 8))
probs$Cluster <- factor(probs$Cluster, levels = c("Pro", "Syn", "Pico1", "Pico2"))



pdf("paper_plots/prob_time.pdf", 12.8, 5.3)
ggplot(probs) + 
  geom_line(aes(Time, Relative_Abundance, color = Model, linetype = Model), linewidth = 0.75) + 
  facet_wrap(~ Cluster, scales = "free") + 
  scale_linetype_manual(values = c("Linear" = "longdash", "Nonlinear" = "solid")) + 
  # scale_y_continuous(limits = c(0, 1)) + 
  scale_x_datetime(date_breaks = "2 days", date_labels = "%b %d") + 
  labs(y = "Relative Abundance") + 
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 30, hjust = 1))
graphics.off()

pdf("paper_plots/prob_time_fixed_scales.pdf", 10.7, 8.8)
ggplot(probs) + 
  geom_line(aes(Time, Relative_Abundance, color = Model)) + 
  facet_wrap(~ Cluster, scales = "fixed") + 
  labs(y = "Relative Abundance") + 
  theme(text = element_text(size = 14))
graphics.off()

#########

dname <- "clustering"
if(!dir.exists(dname)) {
  dir.create(dname)

  sapply(seq_len(296), function(tt) {
    file_lin <- paste0("linear_", tt, ".png")
    path_lin <- file.path(dname, file_lin)

    file_half <- paste0("half_", tt, ".png")
    path_half <- file.path(dname, file_half)

    png(path_lin, width = 1670, height = 556)
    print(flowtrend::plot_3d(ylist, linear_best, tt, countslist))
    graphics.off()

    png(path_half, width = 1670, height = 556)
    print(flowtrend::plot_3d(ylist, half_best, tt, countslist))
    graphics.off()
  })
}

my_plot_3d <- function (ylist, obj = NULL, tt, countslist = NULL, labels = NULL, 
    bin = TRUE, plot_title = NULL, return_list_of_plots = FALSE, cex = 3) 
{
    stopifnot(ncol(ylist[[1]]) == 3)
    if (!is.null(labels)) 
        assertthat::assert_that(length(labels) == obj$numclust)
    labs = colnames(ylist[[1]])
    if (is.null(labs)) 
        labs = paste0("dim", 1:3)
    if (is.null(countslist)) {
        counts = rep(1, nrow(ylist[[1]]))
    }
    if (!is.null(countslist)) {
        counts = countslist[[tt]]
    }
    y2d_list = list()
    counts2d_list = list()
    for (ii in 1:3) {
        dims = list(c(1:2), c(2:3), c(3, 1))[[ii]]
        if (bin) {
            yy = flowmix::collapse_3d_to_2d(ylist[[tt]], counts, 
                dims)
            y2d = yy[, 1:2]
            colnames(y2d) = labs[dims]
            one_counts = yy[, 3]
        }
        else {
            y2d = ylist[[tt]][, dims]
            colnames(y2d) = labs[dims]
            one_counts = counts
        }
        y2d_list[[ii]] = y2d
        counts2d_list[[ii]] = one_counts
    }
    total_range = sapply(counts2d_list, range) %>% range()
    plotlist = list()
    for (ii in 1:3) {
        dims = list(c(1:2), c(2:3), c(3, 1))[[ii]]
        p = flowtrend::plot_2d(list(y2d_list[[ii]]), list(counts2d_list[[ii]]), 
            obj = NULL, tt = 1, bin = bin)
        if (bin) {
            p$scales$scales <- list()
            p = p + scale_fill_gradientn(colors = c("white", 
                "blue", "yellow"), limits = total_range, guide = "none")
        }
        if (!is.null(obj)) {
            mnlist = lapply(1:obj$numclust, function(iclust) {
                one_mnmat = obj$mn[tt, dims, iclust] %>% t()
                colnames(one_mnmat) = paste0("dim", 1:2)
                one_mnmat %>% as_tibble() %>% add_column(cluster = iclust)
            })
            mnmat = do.call(rbind, mnlist)
            mn_colours = rep("red", obj$numclust)
            for (iclust in 1:obj$numclust) {
                el = ellipse::ellipse(x = obj$sigma[iclust, dims, 
                  dims], centre = obj$mn[tt, dims, iclust]) %>% 
                  as_tibble()
                p = p + geom_path(aes(x = x, y = y), data = el, 
                  colour = mn_colours[iclust], lty = 2, linewidth = pmin(obj$prob[tt, 
                    iclust] * 8, 0.8))
                p = p + geom_point(aes(x = dim1, y = dim2), data = mnmat %>% 
                  subset(cluster == iclust), colour = mn_colours[iclust], 
                  size = obj$prob[tt, iclust] * 10)
                cex = rel(cex)
                fac = 10
                if (is.null(labels)) {
                  labels = 1:(obj$numclust)
                }
                dt = data.frame(dim1 = mnmat[, 1], dim2 = mnmat[, 
                  2], prob = obj$prob[tt, ] * fac)
                p = p + ggrepel::geom_text_repel(aes(x = dim1, 
                  y = dim2, label = labels, point.size = sqrt(prob)), 
                  col = mn_colours, cex = cex, bg.color = "white", 
                  bg.r = 0.1, fontface = "bold", force_pull = 5, 
                  data = dt, seed = 1)
            }
        }
        if (ii == 1) {
            if (is.null(plot_title)) 
                plot_title = paste0("Time=", tt)
            p = p + ggtitle(plot_title)
        }
        else {
            p = p + ggtitle("")
        }
        p = p + theme(legend.position = "none")
        plotlist[[ii]] = p
    }
    p_combined = flowtrend::my_mfrow(plotlist)
    if (return_list_of_plots) 
        return(plotlist)
    return(p_combined)
}

my_plot_3d(ylist, half_best, 27, countslist, return_list_of_plots = TRUE, 
  plot_title = "", labels = c("", "", "syn", "", "pico2", "", "pico1", "pro", "", ""), 
  cex = 8)[[1]]

########################

pop_mat <- matrix(c(10, 8, 
                    1, 3, 
                    9, 10,
                    5, 7,
                    7, 5), ncol = 2, byrow = TRUE)
# 1d response plots
models <- c("linear", "half")
pops <- c("pro2", "syn", "pro1", "pico1",
          "pico2")

colnames(pop_mat) <- models
rownames(pop_mat) <- pops

# 1d response plots

# response vs time
plot_1d_resp <- function(pop_mat, pop, resp) {
  resp_dat <- switch(resp, 
                      prob = "$prob[,", 
                      diam = "$mn[,1,", 
                      chl = "$mn[,2,", 
                      pe = "$mn[,3,")
  
  clusts <- pop_mat[pop,]

  lin_resp <- linear_best %>% 
                substitute() %>% 
                deparse() %>%
                paste0(resp_dat, clusts["linear"], "]") %>% 
                str2lang() %>%
                eval()
  half_resp <- half_best %>% 
                substitute() %>% 
                deparse() %>%
                paste0(resp_dat, clusts["half"], "]") %>% 
                str2lang() %>%
                eval()
  
  lin_df <- data.frame(time = time, lin_resp, model = "linear")
  half_df <- data.frame(time = time, half_resp, model = "half")
  colnames(lin_df)[2] <- resp 
  colnames(half_df)[2] <- resp 

  resp_df <- rbind(lin_df, half_df)

  plot_title <- paste0(pop, ": Linear-", clusts["linear"], 
                      ", Half-", clusts["half"], ": ", resp)
  ggplot(resp_df, aes(x = time, y = .data[[resp]], group = model, color = model)) + 
    geom_line() + 
    ggtitle(plot_title) + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title = element_text(size = 16), 
          legend.key.height = unit(1, 'cm'), 
          legend.key.width = unit(1.5, "cm"), 
          legend.text = element_text(size=16), 
          legend.title = element_text(size=16), 
          title = element_text(size = 16))
}

dname <- "response_time"
if(!dir.exists(dname)) {
  dir.create(dname)
  width <- round(480 * 16/9)
  height <- 480
  for(pop in pops) {
    fname <- paste0(pop, ".png")
    fpath <- file.path(dname, fname)
    png(fpath, width * 2, height * 2)
    p1 <- plot_1d_resp(pop_mat, pop, "prob")
    p2 <- plot_1d_resp(pop_mat, pop, "diam")
    p3 <- plot_1d_resp(pop_mat, pop, "chl")
    p4 <- plot_1d_resp(pop_mat, pop, "pe")

    print(grid.arrange(p1, p2, p3, p4, nrow = 2))
    graphics.off()
  }
}

#################################################

# TODO: continue here after finishing the pdpPairs function

# 2d partial response

# kv is a vector of old values with names(kv) equal to the new values
kv_rep <- function(dat, kv) {
  replaced <- names(kv)[match(dat, kv)]
  na_inds <- is.na(replaced)
  replaced[na_inds] <- dat[na_inds]
  return(replaced)
}

get_prfs <- function(cov1, cov2, gridsize, defaults) {
  grid_1 <- seq.int(min(X[,cov1]), max(X[,cov1]), length.out = gridsize)
  grid_2 <- seq.int(min(X[,cov2]), max(X[,cov2]), length.out = gridsize)

  prf_grid <- expand.grid(grid_1, grid_2)
  
  prob_responses <- list()
  prob_responses$linear <- matrix(nrow = gridsize ** 2, ncol = linear_best$numclust)
  prob_responses$half <- matrix(nrow = gridsize ** 2, ncol = half_best$numclust)

  mean_responses <- list()
  mean_responses$linear <- array(dim = c(3, gridsize ** 2, linear_best$numclust))
  mean_responses$half <- array(dim = c(3, gridsize ** 2, half_best$numclust))

  X_new <- tcrossprod(rep(1, gridsize ** 2), defaults)
  colnames(X_new) <- colnames(X_no_cp)
  X_new[,cov1] <- prf_grid$Var1
  X_new[,cov2] <- prf_grid$Var2

  X_lin <- cbind(1, X_new %*% X_pca$rotation[,1:d])
  set.seed(0)
  X_half <- X_hidden(X_new %*% X_pca$rotation[,1:d], d, n.h, 0.5) %>% cbind(1, .)

  for(i in 1:10) {
    prob_responses$linear[,i] <- X_lin %*% linear_best$alpha[i,] 
    prob_responses$half[,i] <- X_half %*% half_best$alpha[i,] 

    mean_responses$linear[,,i] <- t(X_lin %*% linear_best$beta[[i]])
    mean_responses$half[,,i] <- t(X_half %*% half_best$beta[[i]])
  }

  numclust <- linear_best$numclust
  lin_res <- data.frame(prf_grid$Var1, prf_grid$Var2, cluster = rep(1:numclust, each = gridsize ** 2), 
                        model = "linear", 
                        prob = as.vector(prob_responses$linear), 
                        diam = as.vector(mean_responses$linear[1,,]), 
                        chl = as.vector(mean_responses$linear[2,,]), 
                        pe = as.vector(mean_responses$linear[3,,]))

  lin_res$cluster <- kv_rep(lin_res$cluster, pop_mat[,"linear"])

  half_res <- data.frame(prf_grid$Var1, prf_grid$Var2, cluster = rep(1:numclust, each = gridsize ** 2), 
                        model = "half", 
                        prob = as.vector(prob_responses$half), 
                        diam = as.vector(mean_responses$half[1,,]), 
                        chl = as.vector(mean_responses$half[2,,]), 
                        pe = as.vector(mean_responses$half[3,,]))

  half_res$cluster <- kv_rep(half_res$cluster, pop_mat[,"half"])

  res <- rbind(lin_res, half_res)
  colnames(res)[1:2] <- c(cov1, cov2)

  return(res)
} 

get_prfs_2 <- function(cov1, cov2, gridsize, defaults) {
  grid_1 <- seq.int(min(X[,cov1]), max(X[,cov1]), length.out = gridsize)
  grid_2 <- seq.int(min(X[,cov2]), max(X[,cov2]), length.out = gridsize)

  prf_grid <- expand.grid(grid_1, grid_2)
  TT <- nrow(X)
  
  prob_responses <- list()
  prob_responses$linear <- array(dim = c(TT, gridsize ** 2, linear_best$numclust))
  prob_responses$half <- array(dim = c(TT, gridsize ** 2, half_best$numclust))

  mean_responses <- list()
  mean_responses$linear <- array(dim = c(3, TT, gridsize ** 2, linear_best$numclust))
  mean_responses$half <- array(dim = c(3, TT, gridsize ** 2, half_best$numclust))

  for(i in 1:10) {
    for(tt in 1:TT) {
      X_new <- tcrossprod(rep(1, gridsize**2), X_no_cp[tt,])
      colnames(X_new) <- colnames(X_no_cp)
      X_new[,cov1] <- prf_grid$Var1
      X_new[,cov2] <- prf_grid$Var2

      X_lin <- cbind(1, X_new %*% X_pca$rotation[,1:d])
      set.seed(0)
      X_half <- X_hidden(X_new %*% X_pca$rotation[,1:d], d, n.h, 0.5) %>% cbind(1, .)

      prob_responses$linear[tt,,i] <- X_lin %*% linear_best$alpha[i,] 
      prob_responses$half[tt,,i] <- X_half %*% half_best$alpha[i,] 

      mean_responses$linear[,tt,,i] <- t(X_lin %*% linear_best$beta[[i]])
      mean_responses$half[,tt,,i] <- t(X_half %*% half_best$beta[[i]])
    }
  }

  prob_responses <- lapply(prob_responses, function(resp) {
    apply(resp, MARGIN = c(2, 3), mean)
  })

  mean_responses <- lapply(mean_responses, function(resp) {
    apply(resp, MARGIN = c(1, 3, 4), mean)
  })

  numclust <- linear_best$numclust
  lin_res <- data.frame(prf_grid$Var1, prf_grid$Var2, cluster = rep(1:numclust, each = gridsize ** 2), 
                        model = "linear", 
                        prob = as.vector(prob_responses$linear), 
                        diam = as.vector(mean_responses$linear[1,,]), 
                        chl = as.vector(mean_responses$linear[2,,]), 
                        pe = as.vector(mean_responses$linear[3,,]))

  lin_res$cluster <- kv_rep(lin_res$cluster, pop_mat[,"linear"])

  half_res <- data.frame(prf_grid$Var1, prf_grid$Var2, cluster = rep(1:numclust, each = gridsize ** 2), 
                        model = "half", 
                        prob = as.vector(prob_responses$half), 
                        diam = as.vector(mean_responses$half[1,,]), 
                        chl = as.vector(mean_responses$half[2,,]), 
                        pe = as.vector(mean_responses$half[3,,]))

  half_res$cluster <- kv_rep(half_res$cluster, pop_mat[,"half"])

  res <- rbind(lin_res, half_res)
  colnames(res)[1:2] <- c(cov1, cov2)

  return(res)
} 

defaults <- colMeans(X[,setdiff(colnames(X), c("b1", "b2"))])

prp_covs <- c("sst", "sss", "O2", "NO3", "PO4")

dname <- "2d_prps_new"
if(!dir.exists(dname)) {
  dir.create(dname)

  # This creates an expand.grid with no "flipped rows" and no rows where 
  # the values of col1 and cov2 are the same
  cov_grid <- t(combn(prp_covs, 2))
  # layout_mat <- matrix(c(1, 2, 2, 2), nrow = 1)

  apply(cov_grid, 1, function(cur_row) {
    cur_res <- get_prfs(cur_row[1], cur_row[2], 30, defaults)
    cur_res_pop <- cur_res[cur_res$cluster %in% rownames(pop_mat),]

    # TODO: remove the melt call and use the wide format directly.
    # cur_res_pop_long <- melt(cur_res_pop, 
    #                          measure.vars = c("prob", "diam", "chl", "pe"), 
    #                          value.name = "response")
    # probs <- cur_res_pop_long[cur_res_pop_long$variable == "prob",]
    # means <- cur_res_pop_long[cur_res_pop_long$variable != "prob",]
    for(pop in rownames(pop_mat)) {
      fname <- paste0(pop, "_", cur_row[1], "_", cur_row[2], ".png")
      fpath <- file.path(dname, fname)

      pop_clusts <- pop_mat[pop,]
      # TODO: give separate ggtitles for each plot
      title_str <- paste0(pop, ": Linear-", pop_clusts["linear"], 
                          " Half-", pop_clusts["half"], 
                          " Tenth-", pop_clusts["tenth"])
      title_grob <- textGrob(title_str, gp = gpar(fontsize = 18))

      pop_prps <- cur_res_pop[cur_res_pop$cluster == pop,]

      p1 <- ggplot(pop_prps,
                   aes(x = .data[[cur_row[1]]], y = .data[[cur_row[2]]])) +
          geom_tile(aes(fill = prob)) +
          geom_contour(aes(z = prob), bins = 20, color = "red",
                       linetype = "dotted", alpha = 0.85) +
          facet_wrap(~ model, nrow = 3) +
          scale_fill_viridis_c(option = "viridis") +
          theme(legend.position = "right")

      p2 <- ggplot(pop_prps,
                   aes(x = .data[[cur_row[1]]], y = .data[[cur_row[2]]])) +
          geom_tile(aes(fill = diam)) +
          geom_contour(aes(z = diam), bins = 20, color = "black",
                       linetype = "dotted", alpha = 0.85) +
          facet_wrap(~ model, nrow = 3) +
          scale_fill_viridis_c(option = "magma") +
          theme(legend.position = "right")

      p3 <- ggplot(pop_prps,
                   aes(x = .data[[cur_row[1]]], y = .data[[cur_row[2]]])) +
          geom_tile(aes(fill = chl)) +
          geom_contour(aes(z = chl), bins = 20, color = "red",
                       linetype = "dotted", alpha = 0.85) +
          facet_wrap(~ model, nrow = 3) +
          scale_fill_viridis_c(option = "cividis") +
          theme(legend.position = "right")

      p4 <- ggplot(pop_prps,
                   aes(x = .data[[cur_row[1]]], y = .data[[cur_row[2]]])) +
          geom_tile(aes(fill = pe)) +
          geom_contour(aes(z = pe), bins = 20, color = "red",
                       linetype = "dotted", alpha = 0.85) +
          facet_wrap(~ model, nrow = 3) +
          scale_fill_viridis_c(option = "mako") +
          theme(legend.position = "right")

      png(fpath, 1920, 1080)
      print(grid.arrange(p1, p2, p3, p4, top = title_grob, ncol = 4))
      graphics.off()
    }
  })
}

###########################

# Looking more closely at the partial responses, compared with the linear 
# model results

# First, make a new prp plotting function which uses a separate scale for 
# each model and response. Additionally, it uses a different color scheme
# for each response

plot_prps <- function(prps, pop, pop_mat, dname = NULL) {
  pop_clusts <- pop_mat[pop,]
  title_str <- paste0(pop, ": Linear-", pop_clusts["linear"], 
                      " Half-", pop_clusts["half"])
  title_grob <- textGrob(title_str, gp = gpar(fontsize = 18))

  pop_prps <- prps[prps$cluster == pop,]
  lin_prp <- pop_prps %>% 
    filter(model == "linear")
  half_prp <- pop_prps %>% 
    filter(model == "half")

  cov1 <- colnames(pop_prps)[1]
  cov2 <- colnames(pop_prps)[2]

  p1 <- ggplot(lin_prp,
              aes(x = .data[[cov1]], y = .data[[cov2]])) +
    geom_tile(aes(fill = prob)) +
    geom_contour(aes(z = prob), bins = 20, color = "red",
                 linetype = "dotted", alpha = 0.85) + 
    facet_wrap(~ model) +
    scale_fill_viridis_c(option = "viridis") +
    theme(legend.position = "right")

  p2 <- ggplot(half_prp,
             aes(x = .data[[cov1]], y = .data[[cov2]])) +
    geom_tile(aes(fill = prob)) +
    geom_contour(aes(z = prob), bins = 20, color = "red",
                 linetype = "dotted", alpha = 0.85) +
    facet_wrap(~ model) +
    scale_fill_viridis_c(option = "viridis") +
    theme(legend.position = "right")

  p3 <- ggplot(lin_prp,
               aes(x = .data[[cov1]], y = .data[[cov2]])) +
      geom_tile(aes(fill = diam)) +
      geom_contour(aes(z = diam), bins = 20, color = "black",
                   linetype = "dotted", alpha = 0.85) +
      facet_wrap(~ model) +
      scale_fill_viridis_c(option = "magma") +
      theme(legend.position = "right")

  p4 <- ggplot(half_prp,
             aes(x = .data[[cov1]], y = .data[[cov2]])) +
    geom_tile(aes(fill = diam)) +
    geom_contour(aes(z = diam), bins = 20, color = "black",
                 linetype = "dotted", alpha = 0.85) +
    facet_wrap(~ model) +
    scale_fill_viridis_c(option = "magma") +
    theme(legend.position = "right")

  p5 <- ggplot(lin_prp,
             aes(x = .data[[cov1]], y = .data[[cov2]])) +
    geom_tile(aes(fill = chl)) +
    geom_contour(aes(z = chl), bins = 20, color = "red",
                 linetype = "dotted", alpha = 0.85) +
    facet_wrap(~ model) +
    scale_fill_viridis_c(option = "cividis") +
    theme(legend.position = "right")

  p6 <- ggplot(half_prp,
             aes(x = .data[[cov1]], y = .data[[cov2]])) +
    geom_tile(aes(fill = chl)) +
    geom_contour(aes(z = chl), bins = 20, color = "red",
                 linetype = "dotted", alpha = 0.85) +
    facet_wrap(~ model) +
    scale_fill_viridis_c(option = "cividis") +
    theme(legend.position = "right")

  p7 <- ggplot(lin_prp,
             aes(x = .data[[cov1]], y = .data[[cov2]])) +
    geom_tile(aes(fill = pe)) +
    geom_contour(aes(z = pe), bins = 20, color = "red",
                 linetype = "dotted", alpha = 0.85) +
    facet_wrap(~ model) +
    scale_fill_viridis_c(option = "mako") +
    theme(legend.position = "right")

  p8 <- ggplot(half_prp,
             aes(x = .data[[cov1]], y = .data[[cov2]])) +
    geom_tile(aes(fill = pe)) +
    geom_contour(aes(z = pe), bins = 20, color = "red",
                 linetype = "dotted", alpha = 0.85) +
    facet_wrap(~ model) +
    scale_fill_viridis_c(option = "mako") +
    theme(legend.position = "right")

  if(!is.null(dname)) {
    fname <- paste0(pop, "_", cov1, "_", cov2, ".png")
    fname <- file.path(dname, fname)
    png(fname, 1920, 853)
  }
  print(grid.arrange(p1, p3, p5, p7, 
                     p2, p4, p6, p8, top = title_grob, ncol = 4))
  if(!is.null(fname)) graphics.off()
}

#       supplement   half
# pro1       9         10 
# pro2      10         5 
# syn        8         3 
# pico1      5         8/9
# pcio2      2         2

gridsize <- 30
sst_phos <- get_prfs("sst", "phosphate_WOA_clim", gridsize, defaults)
sst_nit <- get_prfs("sst", "nitrate_WOA_clim", gridsize, defaults)
sdes_wind <- get_prfs("surface_downward_eastward_stress", "wind_stress", gridsize, defaults)
east_speed <- get_prfs("eastward_wind", "wind_speed", gridsize, defaults)
o2sat_nit <- get_prfs("o2sat_WOA_clim", "nitrate_WOA_clim", gridsize, defaults)
si_chl <- get_prfs("Si", "CHL", gridsize, defaults)
wind_AOU <- get_prfs("AOU_WOA_clim", "wind_speed", gridsize, defaults)
dens_sla <- get_prfs("density_WOA_clim", "sla", gridsize, defaults)
ugos_disp <- get_prfs("ugos", "disp_bw_sla", gridsize, defaults)
o2sat_phos <- get_prfs("o2sat_WOA_clim", "phosphate_WOA_clim", gridsize, defaults)
pp_chl <- get_prfs("PP", "CHL", gridsize, defaults)
dens_pp <- get_prfs("density_WOA_clim", "PP", gridsize, defaults)

prp_list <- function(prp, pop) {
  list(prp = prp, pop = pop)
}

plots <- list()

plots$plot1 <- prp_list(sst_phos, "pro2")
plots$plot2 <- prp_list(sst_nit, "pro2")
plots$plot3 <- prp_list(sdes_wind, "pro1")
plots$plot4 <- prp_list(east_speed, "pro1")
plots$plot5 <- prp_list(o2sat_nit, "syn")
plots$plot6 <- prp_list(sdes_wind, "syn")
plots$plot7 <- prp_list(si_chl, "syn")
plots$plot8 <- prp_list(wind_AOU, "syn")
plots$plot9 <- prp_list(dens_sla, "pico1")
plots$plot10 <- prp_list(ugos_disp, "pico1")
plots$plot11 <- prp_list(o2sat_phos, "pico1")
plots$plot12 <- prp_list(pp_chl, "pico1")
plots$plot13 <- prp_list(o2sat_phos, "pico2")
plots$plot14 <- prp_list(dens_pp, "pico2")
plots$plot15 <- prp_list(si_chl, "pico2")

lapply(plots, function(cur_plot) {
  plot_prps(cur_plot$prp, cur_plot$pop, pop_mat, "prps_free_scales")
})

####### Testing Friedman's PDP

gridsize <- 30
sst_phos <- get_prfs_2("sst", "phosphate_WOA_clim", gridsize, defaults)
plot_prps(sst_phos, "pro2", pop_mat, "friedman_pdp")

sst_nit <- get_prfs_2("sst", "nitrate_WOA_clim", gridsize, defaults)
sdes_wind <- get_prfs_2("surface_downward_eastward_stress", "wind_stress", gridsize, defaults)
east_speed <- get_prfs_2("eastward_wind", "wind_speed", gridsize, defaults)
o2sat_nit <- get_prfs_2("o2sat_WOA_clim", "nitrate_WOA_clim", gridsize, defaults)
si_chl <- get_prfs_2("Si", "CHL", gridsize, defaults)
wind_AOU <- get_prfs_2("AOU_WOA_clim", "wind_speed", gridsize, defaults)
dens_sla <- get_prfs_2("density_WOA_clim", "sla", gridsize, defaults)
ugos_disp <- get_prfs_2("ugos", "disp_bw_sla", gridsize, defaults)
o2sat_phos <- get_prfs_2("o2sat_WOA_clim", "phosphate_WOA_clim", gridsize, defaults)
pp_chl <- get_prfs_2("PP", "CHL", gridsize, defaults)
dens_pp <- get_prfs_2("density_WOA_clim", "PP", gridsize, defaults)

prp_list <- function(prp, pop) {
  list(prp = prp, pop = pop)
}

plots <- list()

plots$plot1 <- prp_list(sst_phos, "pro2")
plots$plot2 <- prp_list(sst_nit, "pro2")
plots$plot3 <- prp_list(sdes_wind, "pro1")
plots$plot4 <- prp_list(east_speed, "pro1")
plots$plot5 <- prp_list(o2sat_nit, "syn")
plots$plot6 <- prp_list(sdes_wind, "syn")
plots$plot7 <- prp_list(si_chl, "syn")
plots$plot8 <- prp_list(wind_AOU, "syn")
# plots$plot9 <- prp_list(dens_sla, "pico1")
plots$plot10 <- prp_list(ugos_disp, "pico1")
plots$plot11 <- prp_list(o2sat_phos, "pico1")
plots$plot12 <- prp_list(pp_chl, "pico1")
plots$plot13 <- prp_list(o2sat_phos, "pico2")
plots$plot14 <- prp_list(dens_pp, "pico2")
plots$plot15 <- prp_list(si_chl, "pico2")

lapply(plots, function(cur_plot) {
  plot_prps(cur_plot$prp, cur_plot$pop, pop_mat, "friedman_pdp")
})

########################

# Check linear model coefficients

get_alphas <- function(pop, res = lin_best, mat = pop_mat, covs = prp_covs, model = "linear") {
    result <- res$alpha[pop_mat[pop, model], covs]
    t(t(result))
}

get_betas <- function(pop, res = lin_best, mat = pop_mat, covs = prp_covs, model = "linear") {
    result <- res$beta[[pop_mat[pop, model]]][covs,]
    colnames(result) <- c("diam", "chl", "pe")
    result
}

# This is probably terrible coding practice but I won't fix it unless or until I need to share
# this code with someone
get_ab <- function(...) {
    alphas <- get_alphas(...)
    betas <- get_betas(...)

    colnames(alphas) <- "prob"
    print(pop)
    
    as(cbind(alphas, betas), "sparseMatrix")
}

for(pop in pops) {
    print(get_ab(pop))
}

#############################

# Check the behavior of the "dip" in pro2 for all the refits (all entries in bestreslist)

l2_dist <- function(x, y) {
    crossprod(x-y)[1]
}

ref <- cv_sums[["linear"]]$bestreslist[[1]]$mn[141,,8]

# For each type of model and each model in bestreslist, get the cluster with the center
# closest in L2 distance to ref

clusts <- matrix(nrow = 100, ncol = 3)
colnames(clusts) <- models
for(model in models) {
    bestreslist <- cv_sums[[model]]$bestreslist
    for(i in seq_along(bestreslist)) {
        mns <- bestreslist[[i]]$mn[141,,]
        
        clusts[i, model] <- apply(mns, 2, function(mn) {
            l2_dist(mn, ref)
        }) %>% which.min()
    }
}

sink("pro2_bestreslist_clusts.txt")
clusts
sink()

#########################################

# Use the cluster numbers to make the plots 

dname <- "bestreslist_response_time"
if(!dir.exists(dname)) {
  dir.create(dname)
  width <- round(480 * 16/9)
  height <- 480
  for(pop in 1:100) {
    fname <- paste0(pop, ".png")
    fpath <- file.path(dname, fname)
    Cairo(width * 2, height * 2, fpath)
    p1 <- plot_1d_resp(clusts, pop, "prob")
    p2 <- plot_1d_resp(clusts, pop, "diam")
    p3 <- plot_1d_resp(clusts, pop, "chl")
    p4 <- plot_1d_resp(clusts, pop, "pe")

    print(grid.arrange(p1, p2, p3, p4, nrow = 2))
    graphics.off()
  }
}

graphics.off()

plot_1d_resp(clusts, 1, "prob")

# These results are garbage - do something else

# Try pre-specifying the lambdas and fitting each model with different randomly generated
# hidden layer weights and maybe different numbers of hidden nodes

# Also need to go back and redo the three models cross-validation with minMax normalization

####### For poster: plot two of the dimensions

new_ylist <- lapply(ylist, function(x) {
  colnames(x)[1:2] <- c("Log Diameter", "Log Chlorophyll Fluorescence")
  x
})

plots_6 <- flowtrend::plot_3d(new_ylist, half_best, 6, countslist, 
  labels = c("", "", "syn", "", "pico2", "", "pico1", "pro", "", ""), return_list_of_plots = TRUE)

plots_27 <- flowtrend::plot_3d(new_ylist, half_best, 27, countslist, 
  labels = c("", "", "syn", "", "pico2", "", "pico1", "pro", "", ""), return_list_of_plots = TRUE)

plots_6[[1]]
plots_27[[1]]

############# For poster; area plot of probability over time


# linear_clust <- c(10, 1, 9, 5, 7)
# half_clust <- c(8, 3, 10, 7, 5)
# clust_names <- c("pro2", "syn", "pro1", "pico1", "pico2")
# resps <- c("logitProb", "diam", "chl", "pe")

# lin_pred <- pdp_make_pred(X_pc, linear_best)
# half_pred <- pdp_make_pred(X_pc, half_best, n.h = 70)

# for(i in seq_along(resps)) {
  # for(j in seq_along(linear_clust))
# plots <- plot_pdps(lin_pred, linear_clust[j], resps[i], 
                      # pop_limits[[j]][[resps[i]]], 30)
# 
# my_pdpPairs(plots)

## NEW PLOT


# on X-axis should be *_pred$uni_res[[1]]$ice_grid[1:nIce] (nIce = 30 in this case)
# on y-axis should be *_pred$uni_res[[1]]$ice_pred$prob[1:nIce,]

# TODO: currently I'm plotting one ICE. Either plot the ICE at a reasonable point 
# (transition line/median) or compute the PD (preferred) and then plot.

pc1_grid <- half_pred$uni_res[[1]]$ice_grid[1:30]
pc1_pred_logits <- half_pred$uni_res[[1]]$ice_pred$prob[1:30,]

# need to transform linear predictors to abundances

pc1_pred_prob <- apply(pc1_pred_logits, 1, function(cur_row) {
  exp_cur_row <- exp(cur_row)
  exp_cur_row / sum(exp_cur_row)
}) %>% t()

pc1_pred_prob <- pc1_pred_prob[,c(8, 3, 7, 5)]

pc1_pred_pico_prob_half <- pc1_pred_prob[,c(3, 4)]

# pc1_pred_prob[,3] <- pc1_pred_prob[,3] + pc1_pred_prob[,4]
# pc1_pred_prob <- pc1_pred_prob[,1:3]
colnames(pc1_pred_prob) <- c("pro", "syn", "pico1", "pico2")
pc1_pred_prob <- apply(pc1_pred_prob, MARGIN = 1, function(cur_row) cur_row / sum(cur_row)) %>% t()

pc1_pred_prob <- as.data.frame(pc1_pred_prob)
pc1_pred_prob$PC1 <- pc1_grid

pc1_pred_prob_melt <- melt(pc1_pred_prob, id.vars = "PC1", 
                                            variable.name = "Population", 
                                            value.name = "Relative Abundance")
pc1_pred_prob_melt_nl <- pc1_pred_prob_melt
pc1_pred_prob_melt_nl$model <- "nonlinear"

p1 <- ggplot(pc1_pred_prob_melt) + 
  geom_area(aes(x = -PC1, y = .data[["Relative Abundance"]], fill = Population)) + 
  ggtitle("Nonlinear")

p1 <- ggplot(pc1_pred_prob_melt) + 
  geom_line(aes(x = -PC1, y = .data[["Relative Abundance"]], col = Population, group = Population)) + 
  ggtitle("Nonlinear")

leg1 <- ggpubr::get_legend(p1)
leg1 <- ggpubr::as_ggplot(leg1)
p1 <- p1 + theme(legend.position = "none")

##### LINEAR

pc1_grid <- lin_pred$uni_res[[1]]$ice_grid[1:30]
pc1_pred_logits <- lin_pred$uni_res[[1]]$ice_pred$prob[1:30,]

# need to transform linear predictors to abundances

pc1_pred_prob <- apply(pc1_pred_logits, 1, function(cur_row) {
  exp_cur_row <- exp(cur_row)
  exp_cur_row / sum(exp_cur_row)
}) %>% t()

pc1_pred_prob <- pc1_pred_prob[,c(10, 1, 5, 7)]

pc1_pred_pico_prob_lin <- pc1_pred_prob[,c(3, 4)]

# pc1_pred_prob[,3] <- pc1_pred_prob[,3] + pc1_pred_prob[,4]
# pc1_pred_prob <- pc1_pred_prob[,1:3]
colnames(pc1_pred_prob) <- c("pro", "syn", "pico1", "pico2")
pc1_pred_prob <- apply(pc1_pred_prob, MARGIN = 1, function(cur_row) cur_row / sum(cur_row)) %>% t()

pc1_pred_prob <- as.data.frame(pc1_pred_prob)
pc1_pred_prob$PC1 <- pc1_grid

pc1_pred_prob_melt <- melt(pc1_pred_prob, id.vars = "PC1", 
                                            variable.name = "Population", 
                                            value.name = "Relative Abundance")
pc1_pred_prob_melt_lin <- pc1_pred_prob_melt
pc1_pred_prob_melt_lin$model <- "linear"

pc1_pred_prob_melt <- rbind(pc1_pred_prob_melt_lin, pc1_pred_prob_melt_nl)

ggplot(pc1_pred_prob_melt) + 
  geom_line(aes(x = -PC1, y = .data[["Relative Abundance"]], 
            col = Population, linetype = model))

ggplot(pc1_pred_prob_melt) + 
  geom_line(aes(x = -PC1, y = .data[["Relative Abundance"]], linetype = model)) + 
  facet_wrap(~ Population, scales = "free", nrow = 2)


p2 <- ggplot(pc1_pred_prob_melt) + 
  geom_area(aes(x = -PC1, y = .data[["Relative Abundance"]], fill = Population)) + 
  ggtitle("Linear") + 
  theme(legend.position = "none")

p2 <- ggplot(pc1_pred_prob_melt) + 
  geom_line(aes(x = -PC1, y = .data[["Relative Abundance"]], col = Population, group = Population)) + 
  ggtitle("Linear")

grid.arrange(p2, p1, leg1, nrow = 1, widths = c(3, 3, 1), 
            top = textGrob("Predicted Relative Abundances of Populations Across \"Latitude\"", 
                  gp = gpar(fontsize = 18)))

# remaking these graphs with ICEs with the other PCs fixed at "middle" values (mean)

ice_X <- tcrossprod(rep(1, 30), colMeans(X_pc))
ice_X[,1] <- PC1_seq <- seq(min(X_pc[,1]), max(X_pc[,1]), length.out = 30)

set.seed(0)
ice_X_nl <- X_hidden(ice_X, 9, n.h, 0.5)

lin_ice <- predict(linear_best, newx = ice_X, logits = FALSE)
nl_ice <- predict(half_best, newx = ice_X_nl, logits = FALSE)

lin_probs <- lin_ice$prob[,c(10, 1, 5, 7)]
nl_probs <- nl_ice$prob[,c(8, 3, 7, 5)]
colnames(lin_probs) <- c("pro", "syn", "pico1", "pico2")
colnames(nl_probs) <- c("pro", "syn", "pico1", "pico2")
lin_probs <- as.data.frame(lin_probs)
nl_probs <- as.data.frame(nl_probs)
lin_probs$PC1 <- nl_probs$PC1 <- PC1_seq

lin_probs$Response <- nl_probs$Response <- "Relative_Abundance"

lin_probs_melt <- melt(lin_probs, id.vars = c("PC1", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")
nl_probs_melt <- melt(nl_probs, id.vars = c("PC1", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

lin_probs_melt$model <- "linear"
nl_probs_melt$model <- "nonlinear"

lin_diam <- lin_ice$mn[,1,c(10, 1, 5, 7)]
nl_diam <- nl_ice$mn[,1,c(8, 3, 7, 5)]
colnames(lin_diam) <- c("pro", "syn", "pico1", "pico2")
colnames(nl_diam) <- c("pro", "syn", "pico1", "pico2")
lin_diam <- as.data.frame(lin_diam)
nl_diam <- as.data.frame(nl_diam)
lin_diam$PC1 <- nl_diam$PC1 <- PC1_seq

lin_diam$Response <- nl_diam$Response <- "Log_Diam"

lin_diam_melt <- melt(lin_diam, id.vars = c("PC1", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")
nl_diam_melt <- melt(nl_diam, id.vars = c("PC1", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

lin_diam_melt$model <- "linear"
nl_diam_melt$model <- "nonlinear"

lin_chl <- lin_ice$mn[,2,c(10, 1, 5, 7)]
nl_chl <- nl_ice$mn[,2,c(8, 3, 7, 5)]
colnames(lin_chl) <- c("pro", "syn", "pico1", "pico2")
colnames(nl_chl) <- c("pro", "syn", "pico1", "pico2")
lin_chl <- as.data.frame(lin_chl)
nl_chl <- as.data.frame(nl_chl)
lin_chl$PC1 <- nl_chl$PC1 <- PC1_seq

lin_chl$Response <- nl_chl$Response <- "Log_chl"

lin_chl_melt <- melt(lin_chl, id.vars = c("PC1", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")
nl_chl_melt <- melt(nl_chl, id.vars = c("PC1", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

lin_chl_melt$model <- "linear"
nl_chl_melt$model <- "nonlinear"

lin_pe <- lin_ice$mn[,3,c(10, 1, 5, 7)]
nl_pe <- nl_ice$mn[,3,c(8, 3, 7, 5)]
colnames(lin_pe) <- c("pro", "syn", "pico1", "pico2")
colnames(nl_pe) <- c("pro", "syn", "pico1", "pico2")
lin_pe <- as.data.frame(lin_pe)
nl_pe <- as.data.frame(nl_pe)
lin_pe$PC1 <- nl_pe$PC1 <- PC1_seq

lin_pe$Response <- nl_pe$Response <- "Log_pe"

lin_pe_melt <- melt(lin_pe, id.vars = c("PC1", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")
nl_pe_melt <- melt(nl_pe, id.vars = c("PC1", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

lin_pe_melt$model <- "linear"
nl_pe_melt$model <- "nonlinear"

pred_df <- rbind(lin_probs_melt, nl_probs_melt, lin_diam_melt, nl_diam_melt, 
                    lin_chl_melt, nl_chl_melt, lin_pe_melt, nl_pe_melt)

pred_df$Response <- factor(pred_df$Response, levels = c("Relative_Abundance", 
                                                        "Log_Diam", 
                                                        "Log_chl", 
                                                        "Log_pe"))
pred_df$Population <- factor(pred_df$Population, levels = c("pro", "syn", "pico1", 
                                                              "pico2"))

positions <- data.frame(Population = rep(c("pro", "syn", "pico1", "pico2"), each = 4), 
                        PC1 = rep(c(9, -0.5) - 1.75, times = 4), 
                        Relative_Abundance = c(rep(0.375, 2), rep(0.15, 2), rep(0.2, 2), 0.15, 0.125),
                        label = rep(c("31.9°N", "34°N"), times = 4))

positions_2 <- data.frame(Population = "pro", Response = "Relative_Abundance", 
                          PC1 = c(2.5, -2), Prediction = c(0.75, 0.5), 
                          label = c("31.9°N", "34°N"))

positions_2$Response <- factor(positions_2$Response, levels = c("Relative_Abundance", 
                                                        "Log_Diam", 
                                                        "Log_chl", 
                                                        "Log_pe"))
positions_2$Population <- factor(positions_2$Population, levels = c("pro", "syn", "pico1", 
                                                              "pico2"))

ggplot(pred_df) + 
  geom_line(aes(x = -PC1, y = Prediction, linetype = model)) + 
  facet_wrap(~ Population + Response, scales = "free") + 
  scale_linetype_manual(values = c("linear" = "dashed", "nonlinear" = "solid")) + 
  geom_vline(xintercept = -X_pc[203,1], linetype = "dotted") + 
  geom_vline(xintercept = -X_pc[221,1], linetype = "dotted") + 
  geom_text(data = positions_2, aes(-PC1, Prediction, label = label))

ggplot(pc1_pred_prob_melt) + 
  geom_line(aes(x = -PC1, y = Relative_Abundance, 
            col = Population, linetype = model))

ggplot(pc1_pred_prob_melt) + 
  geom_line(aes(x = -PC1, y = .data[["Relative Abundance"]], linetype = model)) + 
  facet_wrap(~ Population, scales = "free", nrow = 2)
###############

# only pro and pico1 (POSTER)

old_probs_melt <- probs_melt 

probs_melt <- old_probs_melt %>% 
  rename(Relative_Abundance = Prediction) %>% 
  select(-Response)

ggplot(probs_melt[probs_melt$Population %in% c("pro", "pico1"),]) + 
  geom_line(aes(x = -PC1, y = Relative_Abundance, linetype = model, color = Population), linewidth = 1) + 
  scale_color_discrete(breaks = NULL) + 
  scale_linetype_manual(values = c("dotdash", "solid")) + 
  facet_wrap(~ Population, scales = "free", nrow = 1, as.table = FALSE) + 
  geom_vline(xintercept = -X_pc[203,1], linetype = "dotted", linewidth = 1) + 
  geom_vline(xintercept = -X_pc[221,1], linetype = "dotted", linewidth = 1) + 
  geom_text(data = positions[positions$Population %in% c("pro", "pico1"),], aes(-PC1, Relative_Abundance, label = label), 
            size = 8) + 
  theme(text = element_text(size = 22), legend.key.size = unit(1, "cm"))

# POSTER/JSM: JUST PRO

positions <- data.frame(Population = "pro", 
                        PC1 = c(8, 0.25) - 1.75,
                        Relative_Abundance = c(0.55, 0.12),
                        label = c("31.9°N", "34°N"))

ggplot(probs_melt[probs_melt$Population == "pro",]) + 
  geom_line(aes(x = -PC1, y = Relative_Abundance, linetype = model, color = Population), linewidth = 1) + 
  scale_color_discrete(breaks = NULL) + 
  scale_linetype_manual(values = c("linear" = "dotdash", "nonlinear" = "solid")) +  
  geom_vline(xintercept = -X_pc[203,1], linetype = "dotted", linewidth = 1) + 
  geom_vline(xintercept = -X_pc[221,1], linetype = "dotted", linewidth = 1) + 
  geom_text(data = positions[positions$Population %in% c("pro", "pico1"),], aes(-PC1, Relative_Abundance, label = label), 
            size = 8) + 
  theme(text = element_text(size = 22), legend.key.size = unit(1, "cm"))

# post-poster (4x4 plot)

ggplot(probs_melt) + 
  geom_line(aes(x = -PC1, y = Relative_Abundance, linetype = model, color = Population), linewidth = 1) + 
  scale_color_discrete(breaks = NULL) + 
  scale_linetype_manual(values = c("dotdash", "solid")) + 
  facet_wrap(~ Population, scales = "free", nrow = 1, as.table = FALSE) + 
  geom_vline(xintercept = -X_pc[203,1], linetype = "dotted", linewidth = 1) + 
  geom_vline(xintercept = -X_pc[221,1], linetype = "dotted", linewidth = 1)

####################

# 

#### Diameters

pc1_grid <- half_pred$uni_res[[1]]$ice_grid[1:30]
pc1_pred_diam <- half_pred$uni_res[[1]]$ice_pred$mn[1:30,1,]

pc1_pred_diam <- pc1_pred_diam[,c(8, 3, 7, 5)]


colnames(pc1_pred_diam) <- c("pro", "syn", "pico1", "pico2")

pc1_pred_diam <- as.data.frame(pc1_pred_diam)
pc1_pred_diam$PC1 <- pc1_grid

pc1_pred_diam_melt <- melt(pc1_pred_diam, id.vars = "PC1", 
                                            variable.name = "Population", 
                                            value.name = "Log Diameter")
pc1_diam_half <- pc1_pred_diam_melt
pc1_diam_half$model <- "nonlinear"

# TODO: make this look nicer

p1 <- ggplot(pc1_pred_diam_melt[pc1_pred_diam_melt$Population == "pro",]) + 
  geom_smooth(aes(x = -PC1, y = .data[["Log Diameter"]]), se = FALSE) +
  ylim(c(1, 2)) + 
  ggtitle("pro")

# ggplot(pc1_pred_diam_melt[pc1_pred_diam_melt$Population %in% c("pico1", "pico2", "aggregatePico"),]) + 
#   geom_smooth(aes(x = -PC1, y = .data[["Log Diameter"]], group = Population, color = Population))

p2 <- ggplot(pc1_pred_diam_melt[pc1_pred_diam_melt$Population == "syn",]) + 
  geom_smooth(aes(x = -PC1, y = .data[["Log Diameter"]]), se = FALSE) + 
  ylab("") +
  ylim(c(2.4, 3.4)) + 
  ggtitle("syn")

p3 <- ggplot(pc1_pred_diam_melt[pc1_pred_diam_melt$Population == "pico1",]) + 
  geom_smooth(aes(x = -PC1, y = .data[["Log Diameter"]]), se = FALSE) + 
  ylab("") +
  ggtitle("pico1")

p4 <- ggplot(pc1_pred_diam_melt[pc1_pred_diam_melt$Population == "pico2",]) + 
  geom_smooth(aes(x = -PC1, y = .data[["Log Diameter"]]), se = FALSE) + 
  ylab("") +
  ggtitle("pico2")

grid.arrange(p1, p2, p3, p4, ncol = 2, top = textGrob("Predicted Log Diameter Across Latitude: Nonlinear Flowmix", 
                                                  gp = gpar(fontsize = 16)))

##### LINEAR

pc1_grid <- lin_pred$uni_res[[1]]$ice_grid[1:30]
pc1_pred_diam <- lin_pred$uni_res[[1]]$ice_pred$mn[1:30,1,]

pc1_pred_diam <- pc1_pred_diam[,c(10, 1, 5, 7)]
# pc1_pred_diam[,3] <- pc1_pred_diam[,3] + pc1_pred_diam[,4]
# pc1_pred_diam <- pc1_pred_diam[,1:3]
colnames(pc1_pred_diam) <- c("pro", "syn", "pico1", "pico2")

pc1_pred_diam <- as.data.frame(pc1_pred_diam)
pc1_pred_diam$PC1 <- pc1_grid

pc1_pred_diam_melt <- melt(pc1_pred_diam, id.vars = "PC1", 
                                            variable.name = "Population", 
                                            value.name = "Log Diameter")
pc1_diam_lin <- pc1_pred_diam_melt
pc1_diam_lin$model <- "linear"


pc1_diam <- rbind(pc1_diam_lin, pc1_diam_half)

# TODO: figure out why they're on such different scales

p1 <- ggplot(pc1_diam[pc1_diam$Population == "pro",]) + 
  geom_smooth(aes(x = -PC1, y = .data[["Log Diameter"]], group = model, color = model), se = FALSE) +
  ggtitle("pro")

leg <- ggpubr::get_legend(p1)
leg <- ggpubr::as_ggplot(leg)

p1 <- p1 + theme(legend.position = "none")

p2 <- ggplot(pc1_diam[pc1_diam$Population == "syn",]) + 
  geom_smooth(aes(x = -PC1, y = .data[["Log Diameter"]], group = model, color = model), se = FALSE) +
  ggtitle("syn") + 
  theme(legend.position = "none")

p3 <- ggplot(pc1_diam[pc1_diam$Population == "pico1",]) + 
  geom_smooth(aes(x = -PC1, y = .data[["Log Diameter"]], group = model, color = model), se = FALSE) +
  ggtitle("pico1") + 
  theme(legend.position = "none")

p4 <- ggplot(pc1_diam[pc1_diam$Population == "pico2",]) + 
  geom_smooth(aes(x = -PC1, y = .data[["Log Diameter"]], group = model, color = model), se = FALSE) +
  ggtitle("pico2") + 
  theme(legend.position = "none")

grid.arrange(p1, p2, p3, p4, leg, layout_matrix = matrix(c(1, 2, 5, 3, 4, 5), nrow = 2, byrow = TRUE), 
              widths = c(3, 3, 1), top = textGrob("Predicted Log Diameter Across Latitude", 
                            gp = gpar(fontsize = 16)))

###### Line Plot of Covariates -> PCs

# FMI 

plot(seq_along(lat), lat)
which.max(lat) # 55
max(lat) # 42.39

which.min(lat) # 296
min(lat) # 22.15

lat[1] # 33.25


# t = 55

# syn
half_best$mn[55,1,3] # 2.85 (slightly better)
linear_best$mn[55,1,1] # 2.96

# syn
half_best$mn[55,1,3] # 2.85 (slightly better)
linear_best$mn[55,1,1] # 2.96

# pico1 (Good fit!)
half_best$mn[55,1,7] # 5.44
linear_best$mn[55,1,5] # 5.45

# pico2
half_best$mn[55,1,5] # 6.58
linear_best$mn[55,1,7] # 6.61

###### New thing: instead of plotting predictions, I;ll plot actual in-sample estimates 

lin_probs <- linear_best$prob[,c(10, 1, 5, 7)]
colnames(lin_probs) <- c("pro", "syn", "pico1", "pico2")
lin_probs <- apply(lin_probs, 1, function(cur_row) cur_row / sum(cur_row)) %>% t()
lin_probs <- as.data.frame(lin_probs)
lin_probs$PC1 <- X_pc[,1]
lin_probs <- melt(lin_probs, id.vars = "PC1", variable.name = "Population", 
                  value.name = "Relative_Abundance")

lin_probs$model <- "linear"

p1 <- ggplot(lin_probs) + 
  geom_area(aes(-PC1, Relative_Abundance, fill = Population)) + 
  ggtitle("Linear")

# and then overlay the PC curves

half_probs <- half_best$prob[,c(8, 3, 7, 5)]
colnames(half_probs) <- c("pro", "syn", "pico1", "pico2")
half_probs <- apply(half_probs, 1, function(cur_row) cur_row / sum(cur_row)) %>% t()
half_probs <- as.data.frame(half_probs)
half_probs$PC1 <- X_pc[,1]
half_probs <- melt(half_probs, id.vars = "PC1", variable.name = "Population", 
                  value.name = "Relative_Abundance")

half_probs$model <- "nonlinear"

p2 <- ggplot(half_probs) + 
  geom_area(aes(-PC1, Relative_Abundance, fill = Population)) + 
  ggtitle("Nonlinear")

grid.arrange(p1, p2, nrow = 1)

# TODO: plot side-by-side with PC curves

set.seed(0)
lin_pred_nlog <- pdp_make_pred(X_pc, linear_best, logits = FALSE)
set.seed(0)
half_pred_nlog <- pdp_make_pred(X_pc, half_best, n.h = 70, logits = FALSE)

lin_ice <- lin_pred_nlog$uni_res[[1]]$ice_pred$prob[,c(10, 1, 5, 7)]
colnames(lin_ice) <- c("pro", "syn", "pico1", "pico2")

lin_pd <- lin_ice %>% 
      array(dim = c(30, 30, 4)) %>% # nIce x gridsize x numclust
      apply(c(1, 3), mean)

colnames(lin_pd) <- c("pro", "syn", "pico1", "pico2")
lin_pd <- apply(lin_pd, 1, function(cur_row) cur_row / sum(cur_row)) %>% t()
lin_pd <- as.data.frame(lin_pd)
lin_pd$PC1 <- lin_pred_nlog$uni_res[[1]]$ice_grid[1:30,]

lin_pd <- melt(lin_pd, id.vars = "PC1", variable.name = "Population", 
                  value.name = "Relative_Abundance")

lin_pd$model <- "linear"

p1 <- ggplot(lin_pd) + 
  geom_area(aes(-PC1, Relative_Abundance, fill = Population)) + 
  ggtitle("Linear")

half_ice <- half_pred_nlog$uni_res[[1]]$ice_pred$prob[,c(8, 3, 7, 5)]
colnames(half_ice) <- c("pro", "syn", "pico1", "pico2")

half_pd <- half_ice %>% 
      array(dim = c(30, 30, 4)) %>% # nIce x gridsize x numclust
      apply(c(1, 3), mean)

colnames(half_pd) <- c("pro", "syn", "pico1", "pico2")
half_pd <- apply(half_pd, 1, function(cur_row) cur_row / sum(cur_row)) %>% t()
half_pd <- as.data.frame(half_pd)
half_pd$PC1 <- half_pred_nlog$uni_res[[1]]$ice_grid[1:30,]

half_pd <- melt(half_pd, id.vars = "PC1", variable.name = "Population", 
                  value.name = "Relative_Abundance")

half_pd$model <- "nonlinear"

p2 <- ggplot(half_pd) + 
  geom_area(aes(-PC1, Relative_Abundance, fill = Population)) + 
  ggtitle("Nonlinear")

grid.arrange(p1, p2, nrow = 1)

#######################

# Now trying ICE where PC1 = 0

# which(abs(X_pc[,1]) < 0.1 ) # 191
# 
# X_ice <- rep(1, 30) %*% X_pc[191,,drop = FALSE]
# X_ice[,1] <- seq(min(X_pc[,1]), max(X_pc[,1]), length.out = 30)
# 
# set.seed(0)
# X_ice_half <- X_hidden(X_ice, 9, n.h, 0.5)
# 
# ice_half <- predict(half_best, logits = FALSE, newx = X_ice_half)
# ice_lin <- predict(linear_best, logits = FALSE, newx = X_ice)
# 
# lin_pred <- ice_lin$prob[,c(10, 1, 5, 7)]
# colnames(lin_pred) <- c("pro", "syn", "pico1", "pico2")
# lin_pred <- apply(lin_pred, 1, function(cur_row) cur_row / sum(cur_row)) %>% t()
# lin_pred <- as.data.frame(lin_pred)
# lin_pred$PC1 <- X_ice[,1]
# 
# lin_pred <- melt(lin_pred, id.vars = "PC1", variable.name = "Population", 
#                   value.name = "Relative_Abundance")
# 
# lin_pred$model <- "Linear"
# 
# p1 <- ggplot(lin_pred) + 
#   geom_area(aes(-PC1, Relative_Abundance, fill = Population)) + 
#   ggtitle("Linear")
# 
# half_pred <- ice_half$prob[,c(8, 3, 7, 5)]
# colnames(half_pred) <- c("pro", "syn", "pico1", "pico2")
# half_pred <- apply(half_pred, 1, function(cur_row) cur_row / sum(cur_row)) %>% t()
# half_pred <- as.data.frame(half_pred)
# half_pred$PC1 <- X_ice[,1]
# 
# half_pred <- melt(half_pred, id.vars = "PC1", variable.name = "Population", 
#                   value.name = "Relative_Abundance")
# 
# half_pred$model <- "nonlinear"
# 
# p2 <- ggplot(half_pred) + 
#   geom_area(aes(-PC1, Relative_Abundance, fill = Population)) + 
#   ggtitle("Nonlinear")
# 
# grid.arrange(p1, p2, nrow = 1)

################# PDP is best

# now for diameter

lin_ice <- lin_pred_nlog$uni_res[[1]]$ice_pred$mn[,1,c(10, 1, 5, 7)]
colnames(lin_ice) <- c("pro", "syn", "pico1", "pico2")

lin_pd <- lin_ice %>% 
      array(dim = c(30, 30, 4)) %>% # nIce x gridsize x numclust
      apply(c(1, 3), mean)

colnames(lin_pd) <- c("pro", "syn", "pico1", "pico2")
lin_pd <- as.data.frame(lin_pd)
lin_pd$PC1 <- lin_pred_nlog$uni_res[[1]]$ice_grid[1:30,]

lin_pd <- melt(lin_pd, id.vars = "PC1", variable.name = "Population", 
                  value.name = "Log_Diameter")

lin_pd$model <- "linear"

# p1 <- ggplot(lin_pd) + 
#   geom_smooth(aes(-PC1, Log_Diameter)) +
#   facet_wrap(~ Population, scales = "free") 
#   ggtitle("Linear")

half_ice <- half_pred_nlog$uni_res[[1]]$ice_pred$mn[,1,c(8, 3, 7, 5)]
colnames(half_ice) <- c("pro", "syn", "pico1", "pico2")

half_pd <- half_ice %>% 
      array(dim = c(30, 30, 4)) %>% # nIce x gridsize x numclust
      apply(c(1, 3), mean)

colnames(half_pd) <- c("pro", "syn", "pico1", "pico2")
half_pd <- as.data.frame(half_pd)
half_pd$PC1 <- half_pred_nlog$uni_res[[1]]$ice_grid[1:30,]

half_pd <- melt(half_pd, id.vars = "PC1", variable.name = "Population", 
                  value.name = "Log_Diameter")

half_pd$model <- "nonlinear"

# p2 <- ggplot(half_pd) + 
#   geom_smooth(aes(-PC1, Log_Diameter)) +
#   facet_wrap(~ Population, scales = "free") 
#   ggtitle("Nonlinear")

diam_pd <- rbind(lin_pd, half_pd)

text_pos <- data.frame(Population = rep("pro", 2), 
                        PC1 = c(0, 4), 
                        Log_Diameter = rep(1.2, 2), 
                        label = c("34°N", "31.9°N"))
ggplot(diam_pd) + 
  geom_smooth(aes(-PC1, Log_Diameter, group = model, color = model)) +
  facet_wrap(~ Population, scales = "free") + 
  geom_vline(xintercept = X_pc[203,1], linetype = "dotdash", linewidth = 0.75) + 
  geom_vline(xintercept = X_pc[221,1], linetype = "dotdash", linewidth = 0.75)#  + 
  geom_text(data = text_pos, aes(-PC1, Log_Diameter, label = label))
  

grid.arrange(p1, p2, nrow = 1)

# Let's compare this with actual estimates: 

lin_diam <- linear_best$mn[,1,c(10, 1, 5, 7)]
colnames(lin_diam) <- c("pro", "syn", "pico1", "pico2")
lin_diam <- as.data.frame(lin_diam)
lin_diam$PC1 <- X_pc[,1]
lin_diam <- melt(lin_diam, id.vars = "PC1", variable.name = "Population", 
                  value.name = "Log_Diameter")

lin_diam$model <- "linear"

p1 <- ggplot(lin_diam) + 
  geom_area(aes(-PC1, Log_Diameter, fill = Population)) + 
  ggtitle("Linear")

# and then overlay the PC curves

half_diam <- half_best$mn[,1,c(8, 3, 7, 5)]
colnames(half_diam) <- c("pro", "syn", "pico1", "pico2")
half_diam <- as.data.frame(half_diam)
half_diam$PC1 <- X_pc[,1]
half_diam <- melt(half_diam, id.vars = "PC1", variable.name = "Population", 
                  value.name = "Log_Diameter")

half_diam$model <- "nonlinear"

diams <- rbind(lin_diam, half_diam)

ggplot(diams) + 
  geom_smooth(aes(-PC1, Log_Diameter, color = model, group = model), se = FALSE) +
  facet_wrap(~ Population, scales = "free") 

########################################

# Relative biomass plot

## Load the model
# datobj <- readRDS("~/repos/flowmix/paper-data/new/MGL1704-hourly-paper-new.RDS")

## Picking populations
pro <- linear_best$prob[,10]
syn <- linear_best$prob[,1]
pico1 <- linear_best$prob[,5]
pico2 <- linear_best$prob[,7]
pico_sum <- linear_best$prob[,c(5, 7)] %>% rowSums()
other1 <- linear_best$prob[,2]
other2 <- linear_best$prob[,3]
other3 <- linear_best$prob[,4]
other4 <- linear_best$prob[,6]
other5 <- linear_best$prob[,8]
other6 <- linear_best$prob[,9]

probmat = tibble(Pro = pro, Syn = syn, Pico1 = pico1, Pico2 = pico2, Total_Pico = pico_sum, 
                 Other1 = other1, Other2 = other2, Other3 = other3, Other4 = other4, 
                 Other5 = other5, Other6 = other6)

## Plot setup
all_cols = RColorBrewer::brewer.pal(8, "Set2")
myColors = c(all_cols[1], all_cols[2], rep(all_cols[3], 3),  rep("grey80", 6))
names(myColors) = c("Pro", "Syn", "Pico1", "Pico2", "Total_Pico", "Other1", "Other2", "Other3", 
                    "Other4", "Other5", "Other6")
color_catalog = tibble(name = names(myColors),
                       col = myColors)
linetypes = c("dashed", "solid")[c(1,1,2,2,1,rep(2,6))]
linetype_catalog = tibble(name = names(myColors),
                          linetype = factor(linetypes))

linewidths = c(1.3, 1)[c(1,1,2,2,1,rep(2,6))]
linewidth_catalog = tibble(name = names(myColors),
                          linewidth = linewidths)
TT = length(pro)
time_catalog = tibble(time=1:TT, date = rownames(X) %>% lubridate::as_datetime())
lat_catalog = tibble(time=1:TT, lat = datobj$lat) %>%
  mutate(transition = ifelse(lat > 30 & lat < 32, TRUE, FALSE)) %>% 
  mutate(lat = (lat-min(lat))) %>% mutate(lat = lat/max(lat))

long_probmat = probmat %>%
  ## add_column(time=datobj$time) %>%
  ## mutate(time = lubridate::as_datetime(time)) %>%
  mutate(time=row_number()) %>% 
  pivot_longer(-time)

## Create the plot
g1 = 
  long_probmat %>%
  left_join(time_catalog, by = "time") %>% 
  left_join(lat_catalog, by = "time") %>% 
  left_join(linetype_catalog, by = "name") %>% 
  left_join(linewidth_catalog, by = "name") %>% 
  ggplot() +
  theme_minimal() +
  geom_line(aes(x=date, y=value, group = name, col = name, linetype = linetype, linewidth = linewidth),
            lineend="round", alpha = .9)  +
  geom_line(aes(x=date, y=lat),
            data =.  %>% group_by(date) %>% summarize(date=unique(date),lat=unique(lat)), col = 'grey35', linetype='dashed') +
  annotate("text", x =  time_catalog %>% .[75,] %>% pull(date), y = 1.04, label = "Ship turns around\n at 42.4°N", lineheight = .8, col = 'grey35')+
  annotate("text", x =  time_catalog %>% .[24,] %>% pull(date), y = .82, label = "Latitude",
           angle = 44, col= 'grey35')+
  scale_color_manual(values=myColors) +
  scale_linewidth(range = c(.5,1)) +
  ylim(c(0,1.05)) +
  geom_vline(xintercept=167, linetype = "dashed", col = 'grey50') +
  ylab("Relative biomass") +
  xlab("") +
  scale_x_datetime(date_breaks = "days")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = time_catalog %>% .[167,] %>% pull(date), col = 'red', linetype=2) +
  ## geom_vline(xintercept = time_catalog %>% .[15,] %>% pull(date), col = 'red', linetype=2) +
  ## geom_vline(xintercept = lat_catalog %>% left_join(time_catalog) %>% subset(transition==TRUE) %>% pull(date), col = 'red', linetype=2) +
  ## geom_vline(xintercept = lat_catalog %>% left_join(time_catalog) %>% subset(transition==TRUE) %>% pull(date), col = 'red', linetype=2) +
  ggtitle("Automatic gating of Seaflow Flow Cytometry data (Gradients 2)")  +
  theme(legend.position="none") +
  annotate("text", x= time_catalog %>% .[167+23,] %>% pull(date), y=0.95, label = "North Pacific \nTransition Zone \n(NPTZ) starts", lineheight = .8)  +
  annotate("text", x= time_catalog %>% .[195,] %>% pull(date),
           y = as.numeric(probmat[195, "Syn"])+.1,
           label = "Synechococcus",
           col = color_catalog %>% subset(name=="Syn") %>% pull(col)) +
  annotate("text", x= time_catalog %>% .[250,] %>% pull(date),
           y = as.numeric(probmat[250, "Pro"])+.1, label = "Prochlorococcus",
           col = color_catalog %>% subset(name=="Pro") %>% pull(col)) +
  annotate("text", x= time_catalog %>% .[50,] %>% pull(date), y = as.numeric(probmat[50, "Total_Pico"])+.15,
           label = "PicoEukaryotes",
           col = color_catalog %>% subset(name=="Total_Pico") %>% pull(col))   +
  annotate("text", x= time_catalog %>% .[285,] %>% pull(date), 
           y = 0.8, label = "C", fontface = "bold", size = 10)
g1
ggsave(g1, file = file.path(plotdir, "rel-biomass.pdf"), width = 9, height = 4.5)

######

# Half

pro <- half_best$prob[,8]
syn <- half_best$prob[,3]
pico1 <- half_best$prob[,7]
pico2 <- half_best$prob[,5]
pico_sum <- half_best$prob[,c(7, 5)] %>% rowSums()
other1 <- half_best$prob[,1]
other2 <- half_best$prob[,2]
other3 <- half_best$prob[,4]
other4 <- half_best$prob[,6]
other5 <- half_best$prob[,9]
other6 <- half_best$prob[,10]

probmat = tibble(Pro = pro, Syn = syn, Pico1 = pico1, Pico2 = pico2, Total_Pico = pico_sum, 
                 Other1 = other1, Other2 = other2, Other3 = other3, Other4 = other4, 
                 Other5 = other5, Other6 = other6)

## Plot setup
all_cols = RColorBrewer::brewer.pal(8, "Set2")
myColors = c(all_cols[1], all_cols[2], rep(all_cols[3], 3),  rep("grey80", 6))
names(myColors) = c("Pro", "Syn", "Pico1", "Pico2", "Total_Pico", "Other1", "Other2", "Other3", 
                    "Other4", "Other5", "Other6")
color_catalog = tibble(name = names(myColors),
                       col = myColors)
linetypes = c("dashed", "solid")[c(1,1,2,2,1,rep(2,6))]
linetype_catalog = tibble(name = names(myColors),
                          linetype = factor(linetypes))

linewidths = c(1.3, 1)[c(1,1,2,2,1,rep(2,6))]
linewidth_catalog = tibble(name = names(myColors),
                          linewidth = linewidths)
TT = length(pro)
time_catalog = tibble(time=1:TT, date = rownames(X) %>% lubridate::as_datetime())
lat_catalog = tibble(time=1:TT, lat = datobj$lat) %>%
  mutate(transition = ifelse(lat > 30 & lat < 32, TRUE, FALSE)) %>% 
  mutate(lat = (lat-min(lat))) %>% mutate(lat = lat/max(lat))

long_probmat = probmat %>%
  ## add_column(time=datobj$time) %>%
  ## mutate(time = lubridate::as_datetime(time)) %>%
  mutate(time=row_number()) %>% 
  pivot_longer(-time)

## Create the plot
g1 = 
  long_probmat %>%
  left_join(time_catalog, by = "time") %>% 
  left_join(lat_catalog, by = "time") %>% 
  left_join(linetype_catalog, by = "name") %>% 
  left_join(linewidth_catalog, by = "name") %>% 
  ggplot() +
  theme_minimal() +
  geom_line(aes(x=date, y=value, group = name, col = name, linetype = linetype, linewidth = linewidth),
            lineend="round", alpha = .9)  +
  geom_line(aes(x=date, y=lat),
            data =.  %>% group_by(date) %>% summarize(date=unique(date),lat=unique(lat)), col = 'grey35', linetype='dashed') +
  annotate("text", x =  time_catalog %>% .[75,] %>% pull(date), y = 1.04, label = "Ship turns around\n at 42.4°N", lineheight = .8, col = 'grey35')+
  annotate("text", x =  time_catalog %>% .[24,] %>% pull(date), y = .82, label = "Latitude",
           angle = 44, col= 'grey35')+
  scale_color_manual(values=myColors) +
  scale_linewidth(range = c(.5,1)) +
  ylim(c(0,1.05)) +
  geom_vline(xintercept=167, linetype = "dashed", col = 'grey50') +
  ylab("Relative biomass") +
  xlab("") +
  scale_x_datetime(date_breaks = "days")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = time_catalog %>% .[167,] %>% pull(date), col = 'red', linetype=2) +
  ## geom_vline(xintercept = time_catalog %>% .[15,] %>% pull(date), col = 'red', linetype=2) +
  ## geom_vline(xintercept = lat_catalog %>% left_join(time_catalog) %>% subset(transition==TRUE) %>% pull(date), col = 'red', linetype=2) +
  ## geom_vline(xintercept = lat_catalog %>% left_join(time_catalog) %>% subset(transition==TRUE) %>% pull(date), col = 'red', linetype=2) +
  ggtitle("Automatic gating of Seaflow Flow Cytometry data (Gradients 2)")  +
  theme(legend.position="none") +
  annotate("text", x= time_catalog %>% .[167+23,] %>% pull(date), y=0.95, label = "North Pacific \nTransition Zone \n(NPTZ) starts", lineheight = .8)  +
  annotate("text", x= time_catalog %>% .[195,] %>% pull(date),
           y = as.numeric(probmat[195, "Syn"])+.1,
           label = "Synechococcus",
           col = color_catalog %>% subset(name=="Syn") %>% pull(col)) +
  annotate("text", x= time_catalog %>% .[250,] %>% pull(date),
           y = as.numeric(probmat[250, "Pro"])+.1, label = "Prochlorococcus",
           col = color_catalog %>% subset(name=="Pro") %>% pull(col)) +
  annotate("text", x= time_catalog %>% .[50,] %>% pull(date), y = as.numeric(probmat[50, "Total_Pico"])+.15,
           label = "PicoEukaryotes",
           col = color_catalog %>% subset(name=="Total_Pico") %>% pull(col))   +
  annotate("text", x= time_catalog %>% .[285,] %>% pull(date), 
           y = 0.8, label = "C", fontface = "bold", size = 10)
g1
ggsave(g1, file = file.path(plotdir, "rel-biomass.pdf"), width = 9, height = 4.5)

######

## Load flowmix object
load("~/01_Cyto/data/paper-data-v2/cvres-2-64-10.Rdata", verbose=TRUE) 
cvres_flowmix = cvres
res_flowmix = cvres_flowmix$bestres
res_flowmix$prob = res_flowmix$pie ## for back-compatilibility

## Picking populations
res_flowmix$prob[,c(10)] -> pro
res_flowmix$prob[,c(3)] -> syn1
res_flowmix$prob[,c(6)] -> syn2
res_flowmix$prob[,c(3,6)] %>% rowSums() -> syn_sum
res_flowmix$prob[,c(5)]->  picoeuk1
res_flowmix$prob[,c(7)]->  picoeuk2
res_flowmix$prob[,c(5,7)] %>% rowSums() ->  picoeuk_sum
res_flowmix$prob[,c(1)]->  bead
res_flowmix$prob[,c(2)]->  other1
res_flowmix$prob[,c(4)]->  other2
res_flowmix$prob[,c(9)]->  other3
res_flowmix$prob[,c(10)]->  other4
probmat = tibble(Pro = pro, Syn1 = syn1, Syn2 = syn2, Total_Syn = syn_sum,
                 Pico1 = picoeuk1, Pico2 = picoeuk2, Total_Pico = picoeuk_sum, Bead=bead,
                 Other1 = other1, Other2 = other2, Other3 = other3, Other4 = other4)

## Plot setup
all_cols = RColorBrewer::brewer.pal(8, "Set2")
myColors = c(all_cols[1], rep(all_cols[2],3), rep(all_cols[3], 3),  rep("grey80", 5))
names(myColors) = c("Pro", "Syn1", "Syn2", "Total_Syn", "Pico1", "Pico2", "Total_Pico", "Bead", "Other1", "Other2", "Other3", "Other4")
color_catalog = tibble(name =
                         c("Pro", "Syn1", "Syn2", "Total_Syn", "Pico1", "Pico2",
                           "Total_Pico", "Bead", "Other1", "Other2", "Other3",
                           "Other4"),
                       col = myColors)
linetypes = c("dashed", "solid")[c(1,2,2,1,2,2,1,rep(2,5))]
linetype_catalog = tibble(name = c("Pro", "Syn1", "Syn2", "Total_Syn", "Pico1",
                                   "Pico2", "Total_Pico", "Bead", "Other1", "Other2",
                                   "Other3", "Other4"),
                          linetype = factor(linetypes))

linewidths = c(1.3, 1)[c(1,2,2,1,2,2,1,rep(2,5))]
linewidth_catalog = tibble(name = c("Pro", "Syn1", "Syn2", "Total_Syn", "Pico1",
                                   "Pico2", "Total_Pico", "Bead", "Other1", "Other2",
                                   "Other3", "Other4"),
                          linewidth = linewidths)
TT = length(pro)
time_catalog = tibble(time=1:TT, date = rownames(res_flowmix$X) %>% lubridate::as_datetime())
lat_catalog = tibble(time=1:TT, lat = datobj$lat) %>%
  mutate(transition = ifelse(lat > 30 & lat < 32, TRUE, FALSE)) %>% 
  mutate(lat = (lat-min(lat))) %>% mutate(lat = lat/max(lat))

long_probmat = probmat %>%
  ## add_column(time=datobj$time) %>%
  ## mutate(time = lubridate::as_datetime(time)) %>%
  mutate(time=row_number()) %>% 
  pivot_longer(-time)

## Create the plot
g1 = 
  long_probmat %>%
  left_join(time_catalog, by = "time") %>% 
  left_join(lat_catalog, by = "time") %>% 
  left_join(linetype_catalog, by = "name") %>% 
  left_join(linewidth_catalog, by = "name") %>% 
  ggplot() +
  theme_minimal() +
  geom_line(aes(x=date, y=value, group = name, col = name, linetype = linetype, linewidth = linewidth),
            lineend="round", alpha = .9)  +
  geom_line(aes(x=date, y=lat),
            data =.  %>% group_by(date) %>% summarize(date=unique(date),lat=unique(lat)), col = 'grey35', linetype='dashed') +
  annotate("text", x =  time_catalog %>% .[75,] %>% pull(date), y = 1.04, label = "Ship turns around\n at 42.4°N", lineheight = .8, col = 'grey35')+
  annotate("text", x =  time_catalog %>% .[24,] %>% pull(date), y = .82, label = "Latitude",
           angle = 44, col= 'grey35')+
  scale_color_manual(values=myColors) +
  scale_linewidth(range = c(.5,1)) +
  ylim(c(0,1.05)) +
  geom_vline(xintercept=167, linetype = "dashed", col = 'grey50') +
  ylab("Relative biomass") +
  xlab("") +
  scale_x_datetime(date_breaks = "days")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = time_catalog %>% .[167,] %>% pull(date), col = 'red', linetype=2) +
  ## geom_vline(xintercept = time_catalog %>% .[15,] %>% pull(date), col = 'red', linetype=2) +
  ## geom_vline(xintercept = lat_catalog %>% left_join(time_catalog) %>% subset(transition==TRUE) %>% pull(date), col = 'red', linetype=2) +
  ## geom_vline(xintercept = lat_catalog %>% left_join(time_catalog) %>% subset(transition==TRUE) %>% pull(date), col = 'red', linetype=2) +
  ggtitle("Automatic gating of Seaflow Flow Cytometry data (Gradients 2)")  +
  theme(legend.position="none") +
  annotate("text", x= time_catalog %>% .[167+23,] %>% pull(date), y=0.95, label = "North Pacific \nTransition Zone \n(NPTZ) starts", lineheight = .8)  +
  annotate("text", x= time_catalog %>% .[195,] %>% pull(date),
           y = as.numeric(probmat[195, "Total_Syn"])+.1,
           label = "Synechococcus",
           col = color_catalog %>% subset(name=="Total_Syn") %>% pull(col)) +
  annotate("text", x= time_catalog %>% .[250,] %>% pull(date),
           y = as.numeric(probmat[250, "Pro"])+.1, label = "Prochlorococcus",
           col = color_catalog %>% subset(name=="Pro") %>% pull(col)) +
  annotate("text", x= time_catalog %>% .[50,] %>% pull(date), y = as.numeric(probmat[50, "Total_Pico"])+.15,
           label = "PicoEukaryotes",
           col = color_catalog %>% subset(name=="Total_Pico") %>% pull(col))   +
  annotate("text", x= time_catalog %>% .[285,] %>% pull(date), 
           y = 0.8, label = "C", fontface = "bold", size = 10)
g1
ggsave(g1, file = file.path(plotdir, "rel-biomass.pdf"), width = 9, height = 4.5)


###### PCA Line Plots

# PC1

lat_scaled <- (lat - min(lat)) / (max(lat) - min(lat)) * (10) - 1.5

cov_pc_df <- cbind(X, X_pc)
cov_pc_df <- as.data.frame(cov_pc_df)
cov_pc_df$time <- lubridate::as_datetime(rownames(cov_pc_df))
cov_pc_df$lat <- lat_scaled

colnames(cov_pc_df)[c(26, 36, 37)] <- c("sdns", "nitrate", "phosphate")

cov_pc_df_long <- melt(cov_pc_df, id.vars = "time", variable.name = "Variable", value.name = "Value")
head(cov_pc_df_long)

pc1_covs <- c("sss", "sst", "PP", "Si", "NO3", "CHL", "PHYC", "PO4", "O2", "PC1")
pc1_covs_lat <- c(pc1_covs, "lat")

ggplot(cov_pc_df_long[cov_pc_df_long$Variable %in% pc1_covs_lat,]) + 
  geom_path(data = cov_pc_df_long[cov_pc_df_long$Variable %in% pc1_covs,], 
          aes(time, Value, group = Variable, color = Variable, 
              linewidth = Variable)) + 
  scale_x_datetime(date_breaks = "2 days") + 
  guides(x = guide_axis(angle = 15)) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(9, "Set1"), "black")) + 
  scale_linewidth_manual(values = c(rep(1, 9), 2)) + 
  geom_path(data = cov_pc_df_long[cov_pc_df_long$Variable == "lat",], 
            aes(time, Value), color = "black", linetype = "dashed", linewidth = 1) + 
  annotate("text", x = cov_pc_df$time[130], y = 8, label = "latitude (for comparison;\nnot used as a covariate)", 
            fontface = "bold", size = 5) + 
  geom_vline(xintercept = cov_pc_df$time[203], linetype = "dotdash", linewidth = 1) + 
  geom_vline(xintercept = cov_pc_df$time[221], linetype = "dotdash", linewidth = 1) + 
  annotate("text", x = cov_pc_df$time[265], y = -3.5, 
          label = "latitudes 34°N-31.9°N;\nconditions change\nabruptly", 
          fontface = "bold", size = 5) + 
  annotate("text", x = cov_pc_df$time[140], y = -5.5, 
            label = "1st Principal Component:\n\"negative latitude,\"\nlearned from\n environmental conditions", 
            fontface = "bold", size = 5) + 
  theme(text = element_text(size = 22), legend.key.size = unit(1, "cm"))

###### PC2

pc2_covs <- c("Fe", "sla", "wind_stress", "wind_speed", "vgos", "vgosa", 
              "sdns", "northward_wind", "nitrate", "phosphate", "PC2")

pdf("paper_plots/PC2_line.pdf", 14, 7.9)
ggplot(cov_pc_df_long[cov_pc_df_long$Variable %in% pc2_covs,]) + 
  geom_path(data = cov_pc_df_long[cov_pc_df_long$Variable %in% pc2_covs,], 
          aes(time, Value, group = Variable, color = Variable, 
              linewidth = Variable)) + 
  scale_x_datetime(date_breaks = "2 days") + 
  guides(x = guide_axis(angle = 15)) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(10, "Set3"), "black")) + 
  scale_linewidth_manual(values = c(rep(1, 10), 2)) + 
  geom_vline(xintercept = cov_pc_df$time[203], linetype = "dotdash", linewidth = 1) + 
  geom_vline(xintercept = cov_pc_df$time[221], linetype = "dotdash", linewidth = 1) + 
  # annotate("text", x = cov_pc_df$time[140], y = -5.5, 
  #           label = "1st Principal Component:\n\"negative latitude,\"\nlearned from\n environmental conditions", 
  #           fontface = "bold", size = 5) + 
  theme(text = element_text(size = 22), legend.key.size = unit(1, "cm"))
graphics.off()

###### PC3

pc3_covs <- c("p1", "p2", "ugos", "ugosa", "PC3")

pdf("paper_plots/PC3_line.pdf", 14, 7.9)
ggplot(cov_pc_df_long[cov_pc_df_long$Variable %in% pc3_covs,]) + 
  geom_path(data = cov_pc_df_long[cov_pc_df_long$Variable %in% pc3_covs,], 
          aes(time, Value, group = Variable, color = Variable, 
              linewidth = Variable)) + 
  scale_x_datetime(date_breaks = "2 days") + 
  guides(x = guide_axis(angle = 15)) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(4, "Set3"), "black")) + 
  scale_linewidth_manual(values = c(rep(1, 4), 2)) + 
  geom_vline(xintercept = cov_pc_df$time[203], linetype = "dotdash", linewidth = 1) + 
  geom_vline(xintercept = cov_pc_df$time[221], linetype = "dotdash", linewidth = 1) + 
  # annotate("text", x = cov_pc_df$time[140], y = -5.5, 
  #           label = "1st Principal Component:\n\"negative latitude,\"\nlearned from\n environmental conditions", 
  #           fontface = "bold", size = 5) + 
  theme(text = element_text(size = 22), legend.key.size = unit(1, "cm"))
graphics.off()

###### PC4

pc4_covs <- c("p3", "p4", "par", "PC4")

pdf("paper_plots/PC4_line.pdf", 14, 7.9)
ggplot(cov_pc_df_long[cov_pc_df_long$Variable %in% pc4_covs,]) + 
  geom_path(data = cov_pc_df_long[cov_pc_df_long$Variable %in% pc4_covs,], 
          aes(time, Value, group = Variable, color = Variable, 
              linewidth = Variable)) + 
  scale_x_datetime(date_breaks = "2 days") + 
  guides(x = guide_axis(angle = 15)) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(3, "Set3"), "black")) + 
  scale_linewidth_manual(values = c(rep(1, 3), 2)) + 
  geom_vline(xintercept = cov_pc_df$time[203], linetype = "dotdash", linewidth = 1) + 
  geom_vline(xintercept = cov_pc_df$time[221], linetype = "dotdash", linewidth = 1) + 
  # annotate("text", x = cov_pc_df$time[140], y = -5.5, 
  #           label = "1st Principal Component:\n\"negative latitude,\"\nlearned from\n environmental conditions", 
  #           fontface = "bold", size = 5) + 
  theme(text = element_text(size = 22), legend.key.size = unit(1, "cm"))
graphics.off()

################################

# look at predictions at low latitude, transition zone, high latitude

ice_X <- tcrossprod(rep(1, 30 * 3), colMeans(X_pc))
PC2_seq <- seq(min(X_pc[,2]), max(X_pc[,2]), length.out = 30)
ice_X[,2] <- rep(PC2_seq, times = 3)
PC1_levels <- c(-4, 1, 6)
ice_X[,1] <- rep(PC1_levels, each = 30)

set.seed(0)
ice_X_nl <- X_hidden(ice_X, 9, n.h, 0.5)

lin_ice <- predict(linear_best, newx = ice_X, logits = FALSE)
nl_ice <- predict(half_best, newx = ice_X_nl, logits = FALSE)

lin_diam <- lin_ice$mn[,1,c(10, 1, 5, 7)]
nl_diam <- nl_ice$mn[,1,c(8, 3, 7, 5)]
colnames(lin_diam) <- c("pro", "syn", "pico1", "pico2")
colnames(nl_diam) <- c("pro", "syn", "pico1", "pico2")

lin_diam <- as.data.frame(lin_diam)
nl_diam <- as.data.frame(nl_diam)
lin_diam$PC2 <- nl_diam$PC2 <- ice_X[,2]
lin_diam$PC1 <- nl_diam$PC1 <- factor(ice_X[,1], levels = c(-4, 1, 6), labels = c("High Latitude", 
                                                                                  "Transition Zone", 
                                                                                  "Low Latitude"))

lin_diam$Response <- nl_diam$Response <- "Log_Diam"

lin_diam_melt <- melt(lin_diam, id.vars = c("PC1", "PC2", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

nl_diam_melt <- melt(nl_diam, id.vars = c("PC1", "PC2", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

lin_diam_melt$model <- "linear"
nl_diam_melt$model <- "nonlinear"

lin_chl <- lin_ice$mn[,2,c(10, 1, 5, 7)]
nl_chl <- nl_ice$mn[,2,c(8, 3, 7, 5)]
colnames(lin_chl) <- c("pro", "syn", "pico1", "pico2")
colnames(nl_chl) <- c("pro", "syn", "pico1", "pico2")

lin_chl <- as.data.frame(lin_chl)
nl_chl <- as.data.frame(nl_chl)
lin_chl$PC2 <- nl_chl$PC2 <- ice_X[,2]
lin_chl$PC1 <- nl_chl$PC1 <- factor(ice_X[,1], levels = c(-4, 1, 6), labels = c("High Latitude", 
                                                                                  "Transition Zone", 
                                                                                  "Low Latitude"))

lin_chl$Response <- nl_chl$Response <- "Log_chl"

lin_chl_melt <- melt(lin_chl, id.vars = c("PC1", "PC2", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

nl_chl_melt <- melt(nl_chl, id.vars = c("PC1", "PC2", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

lin_chl_melt$model <- "linear"
nl_chl_melt$model <- "nonlinear"

lin_pe <- lin_ice$mn[,3,c(10, 1, 5, 7)]
nl_pe <- nl_ice$mn[,3,c(8, 3, 7, 5)]
colnames(lin_pe) <- c("pro", "syn", "pico1", "pico2")
colnames(nl_pe) <- c("pro", "syn", "pico1", "pico2")

lin_pe <- as.data.frame(lin_pe)
nl_pe <- as.data.frame(nl_pe)
lin_pe$PC2 <- nl_pe$PC2 <- ice_X[,2]
lin_pe$PC1 <- nl_pe$PC1 <- factor(ice_X[,1], levels = c(-4, 1, 6), labels = c("High Latitude", 
                                                                                  "Transition Zone", 
                                                                                  "Low Latitude"))

lin_pe$Response <- nl_pe$Response <- "Log_pe"

lin_pe_melt <- melt(lin_pe, id.vars = c("PC1", "PC2", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

nl_pe_melt <- melt(nl_pe, id.vars = c("PC1", "PC2", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

lin_pe_melt$model <- "linear"
nl_pe_melt$model <- "nonlinear"

lin_prob <- lin_ice$prob[,c(10, 1, 5, 7)]
nl_prob <- nl_ice$prob[,c(8, 3, 7, 5)]
colnames(lin_prob) <- c("pro", "syn", "pico1", "pico2")
colnames(nl_prob) <- c("pro", "syn", "pico1", "pico2")

lin_prob <- as.data.frame(lin_prob)
nl_prob <- as.data.frame(nl_prob)
lin_prob$PC2 <- nl_prob$PC2 <- ice_X[,2]
lin_prob$PC1 <- nl_prob$PC1 <- factor(ice_X[,1], levels = c(-4, 1, 6), labels = c("High Latitude", 
                                                                                  "Transition Zone", 
                                                                                  "Low Latitude"))

lin_prob$Response <- nl_prob$Response <- "Relative_Abundance"

lin_prob_melt <- melt(lin_prob, id.vars = c("PC1", "PC2", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

nl_prob_melt <- melt(nl_prob, id.vars = c("PC1", "PC2", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

lin_prob_melt$model <- "linear"
nl_prob_melt$model <- "nonlinear"

prob_melt <- rbind(lin_prob_melt, nl_prob_melt)

pred_df <- rbind(lin_prob_melt, nl_prob_melt, lin_diam_melt, nl_diam_melt, lin_chl_melt, nl_chl_melt, lin_pe_melt, nl_pe_melt)
pred_df$Response <- factor(pred_df$Response, levels = c("Relative_Abundance", 
                                                        "Log_Diam", 
                                                        "Log_chl", 
                                                        "Log_pe"))

ggplot(pred_df) + 
  geom_line(aes(x = PC2, y = Prediction, linetype = model, color = PC1)) + 
  facet_wrap(Population ~ Response, nrow = 4, scales = "free") + 
  scale_linetype_manual(values = c("linear" = "dashed", "nonlinear" = "solid"))
 
diam_melt <- rbind(lin_diam_melt, nl_diam_melt)

ggplot(diam_melt[diam_melt$Population == "pico1",]) + 
  geom_line(aes(x = PC2, y = Relative_Abundance, linetype = model, color = PC1), linewidth = 1) + 
  facet_wrap(~ Population) + 
  theme(text = element_text(size = 22), legend.key.size = unit(1, "cm")) + 
  ylab("")

    # scale_color_discrete(breaks = NULL) + 
  # scale_linetype_manual(values = c("dotdash", "solid")) + 
  # facet_wrap(~ Population, scales = "free", nrow = 2, as.table = FALSE)# + 
  # geom_vline(xintercept = -X_pc[203,1], linetype = "dotted") + 
  # geom_vline(xintercept = -X_pc[221,1], linetype = "dotted") + 
  # geom_text(data = positions, aes(-PC1, Relative_Abundance, label = label))

###### now vary PC1

ice_X <- tcrossprod(rep(1, 30), colMeans(X_pc))
ice_X[,1] <- PC1_seq <- seq(min(X_pc[,1]), max(X_pc[,1]), length.out = 30)

set.seed(0)
ice_X_nl <- X_hidden(ice_X, 9, n.h, 0.5)

lin_ice <- predict(linear_best, newx = ice_X, logits = FALSE)
nl_ice <- predict(half_best, newx = ice_X_nl, logits = FALSE)

lin_diam <- lin_ice$mn[,1,c(10, 1, 5, 7)]
nl_diam <- nl_ice$mn[,1,c(8, 3, 7, 5)]
colnames(lin_diam) <- c("pro", "syn", "pico1", "pico2")
colnames(nl_diam) <- c("pro", "syn", "pico1", "pico2")
lin_diam <- as.data.frame(lin_diam)
nl_diam <- as.data.frame(nl_diam)
lin_diam$PC1 <- nl_diam$PC1 <- PC1_seq

lin_diam_melt <- melt(lin_diam, id.vars = "PC1", 
                                  variable.name = "Population", 
                                  value.name = "Log_Diameter")
nl_diam_melt <- melt(nl_diam, id.vars = "PC1", 
                                  variable.name = "Population", 
                                  value.name = "Log_Diameter")

lin_diam_melt$model <- "linear"
nl_diam_melt$model <- "nonlinear"

diam_melt <- rbind(lin_diam_melt, nl_diam_melt)

positions <- data.frame(Population = rep(c("pro", "syn", "pico1", "pico2"), each = 2), 
                        PC1 = rep(c(5, 0) - 1.75, times = 4), 
                        Log_Diameter = c(rep(0.8, 2), rep(0.15, 2), rep(0.2, 2), 0.15, 0.125),
                        label = rep(c("31.9°N", "34°N"), times = 4))


# POSTER/JSM  PLOTS
ggplot(diam_melt[diam_melt$Population == "pico1",]) + 
  geom_line(aes(x = -PC1, y = Log_Diameter, linetype = model, color = Population), linewidth = 1) + 
  scale_color_discrete(breaks = NULL) + 
  scale_linetype_manual(values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
  facet_wrap(~ Population) + 
  geom_vline(xintercept = -X_pc[203,1], linetype = "dotted") + 
  geom_vline(xintercept = -X_pc[221,1], linetype = "dotted") +
  annotate("text", x = -6, y = 5.7, label = "31.9°N", size = 8) +
  annotate("text", x = 1.25, y = 5.7, label = "34°N", size = 8) + 
  facet_wrap(~ Population) + 
  # annotate("text", x = 4.5, y = 5.55, label = "Picoeukaryotes", size = 8) + 
  theme(text = element_text(size = 22))
  # geom_text(data = positions, aes(-PC1, Relative_Abundance, label = label))

# I don't know what this is
ggplot(prob_melt[prob_melt$Population == "pro",]) + 
  geom_line(aes(x = PC2, y = Relative_Abundance, linetype = model, color = Population), linewidth = 1) + 
  scale_color_discrete(breaks = NULL) + 
  scale_linetype_manual(values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
  # facet_wrap(~ Population, scales = "free", nrow = 2, as.table = FALSE) + 
  geom_vline(xintercept = -X_pc[203,1], linetype = "dotted") + 
  geom_vline(xintercept = -X_pc[221,1], linetype = "dotted") +
  annotate("text", x = -6, y = 5.7, label = "31.9°N", size = 8) +
  annotate("text", x = 1.25, y = 5.7, label = "34°N", size = 8) + 
  facet_wrap(~ Population) + 
  # annotate("text", x = 4.5, y = 5.55, label = "Picoeukaryotes", size = 8) + 
  theme(text = element_text(size = 22))


old_prob_melt <- prob_melt 
prob_melt <- prob_melt %>% 
  rename(Relative_Abundance = Prediction) %>% 
  select(-Response)

# JSM

png("pro_pc2_prob.png", 10.9, 9.1, "in", res = 300)
ggplot(prob_melt[prob_melt$Population == "pro",]) + 
  geom_line(aes(x = PC2, y = Relative_Abundance, linetype = model, color = PC1), linewidth = 1) + 
  scale_linetype_manual(values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
  facet_wrap(~ Population, scales = "free", nrow = 2, as.table = FALSE) + 
  facet_wrap(~ Population) + 
  theme(text = element_text(size = 22))
graphics.off()

ggplot(diam_melt[diam_melt$Population == "syn",]) + 
  geom_line(aes(x = -PC1, y = Log_Diameter, linetype = model, color = Population), linewidth = 1) + 
  scale_color_discrete(breaks = NULL) + 
  scale_linetype_manual(values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
  facet_wrap(~ Population) + 
  geom_vline(xintercept = -X_pc[203,1], linetype = "dotted") + 
  geom_vline(xintercept = -X_pc[221,1], linetype = "dotted") +
  annotate("text", x = -6, y = 5.7, label = "31.9°N", size = 8) +
  annotate("text", x = 1.25, y = 5.7, label = "34°N", size = 8) + 
  facet_wrap(~ Population) + 
  # annotate("text", x = 4.5, y = 5.55, label = "Picoeukaryotes", size = 8) + 
  theme(text = element_text(size = 22))

lin_chl_melt$model <- "linear"
nl_chl_melt$model <- "nonlinear"
chl_melt <- rbind(lin_chl_melt, nl_chl_melt)
chl_melt <- chl_melt %>% 
  rename(Log_chl = Prediction) %>% 
  select(-Response)

# JSM
png("syn_pc1_chl.png", 8.9, 7.9, "in", res = 300)
ggplot(chl_melt[chl_melt$Population == "syn",]) + 
  geom_line(aes(x = -PC1, y = Log_chl, linetype = model, color = Population), linewidth = 1) + 
  scale_color_discrete(breaks = NULL) + 
  scale_linetype_manual(values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
  facet_wrap(~ Population) + 
  geom_vline(xintercept = -X_pc[203,1], linetype = "dotted") + 
  geom_vline(xintercept = -X_pc[221,1], linetype = "dotted") +
  # annotate("text", x = -6, y = 5.7, label = "31.9°N", size = 8) +
  # annotate("text", x = 1.25, y = 5.7, label = "34°N", size = 8) + 
  facet_wrap(~ Population) + 
  # annotate("text", x = 4.5, y = 5.55, label = "Picoeukaryotes", size = 8) + 
  theme(text = element_text(size = 22))
graphics.off()

# TO REPRODUCE THE ABOVE PLOT: 
# Run lines 1406-1413, 1455-1473, 2579-2582

# Need to make PC1.svg and PC1_PCX.csv for X = 2, 3, 4
# (it's not clear where these files came from)

##################################

# PAPER PLOTS

# First, exploration. We want ICEs for every response, every population
#    -PC1
#     PC2, PC3, and PC4 for three different levels of PC1 (low, transition, high)

# Create fine grids of X values
PC1_seq <- seq(min(X_pc[,1]), max(X_pc[,1]), length.out = 30)
PC2_seq <- seq(min(X_pc[,2]), max(X_pc[,2]), length.out = 30)
PC3_seq <- seq(min(X_pc[,3]), max(X_pc[,3]), length.out = 30)
PC4_seq <- seq(min(X_pc[,4]), max(X_pc[,4]), length.out = 30)

# Start with all PCs fixed at their mean
PC1_ice_X <- tcrossprod(rep(1, 30), colMeans(X_pc))
# vary PC1 
PC1_ice_X[,1] <- PC1_seq

# Start with all PCs fixed at their mean
other_ice_X <- tcrossprod(rep(1, 30 * 3), colMeans(X_pc))

# check the ICEs of PC2, PC3, and PC4 with PC1 fixed at values corresponding 
# to low, medium, and high latitude
PC1_levels <- c(-4, 1, 6)
other_ice_X[,1] <- rep(PC1_levels, each = 30)

PC2_ice_X <- PC3_ice_X <- PC4_ice_X <- other_ice_X

PC2_ice_X[,2] <- rep(PC2_seq, times = 3)
PC3_ice_X[,3] <- rep(PC3_seq, times = 3)
PC4_ice_X[,4] <- rep(PC4_seq, times = 3)

set.seed(0)
PC1_ice_X_nl <- X_hidden(PC1_ice_X, 9, n.h, 0.5)
set.seed(0)
PC2_ice_X_nl <- X_hidden(PC2_ice_X, 9, n.h, 0.5)
set.seed(0)
PC3_ice_X_nl <- X_hidden(PC3_ice_X, 9, n.h, 0.5)
set.seed(0)
PC4_ice_X_nl <- X_hidden(PC4_ice_X, 9, n.h, 0.5)

PC1_lin_ice <- predict(linear_best, newx = PC1_ice_X, logits = FALSE)
PC1_nl_ice <- predict(half_best, newx = PC1_ice_X_nl, logits = FALSE)
PC2_lin_ice <- predict(linear_best, newx = PC2_ice_X, logits = FALSE)
PC2_nl_ice <- predict(half_best, newx = PC2_ice_X_nl, logits = FALSE)
PC3_lin_ice <- predict(linear_best, newx = PC3_ice_X, logits = FALSE)
PC3_nl_ice <- predict(half_best, newx = PC3_ice_X_nl, logits = FALSE)
PC4_lin_ice <- predict(linear_best, newx = PC4_ice_X, logits = FALSE)
PC4_nl_ice <- predict(half_best, newx = PC4_ice_X_nl, logits = FALSE)

# TODO: functionize attaching & formatting the predictions 
# (because this is where it starts to get really tedious)

PC1_lin_probs <- PC1_lin_ice$prob[,c(10, 1, 5, 7)]
PC1_nl_probs <- PC1_nl_ice$prob[,c(8, 3, 7, 5)]
colnames(PC1_lin_probs) <- c("pro", "syn", "pico1", "pico2")
colnames(PC1_nl_probs) <- c("pro", "syn", "pico1", "pico2")
PC1_lin_probs <- as.data.frame(PC1_lin_probs)
PC1_nl_probs <- as.data.frame(PC1_nl_probs)
PC1_lin_probs$PC1 <- PC1_nl_probs$PC1 <- PC1_seq
PC1_lin_probs$Response <- PC1_nl_probs$Response <- "Relative_Abundance"
PC1_lin_probs$model <- "linear" 
PC1_nl_probs$model <- "nonlinear"

PC1_lin_diam <- PC1_lin_ice$mn[,1,c(10, 1, 5, 7)]
PC1_nl_diam <- PC1_nl_ice$mn[,1,c(8, 3, 7, 5)]
colnames(PC1_lin_diam) <- c("pro", "syn", "pico1", "pico2")
colnames(PC1_nl_diam) <- c("pro", "syn", "pico1", "pico2")
PC1_lin_diam <- as.data.frame(PC1_lin_diam)
PC1_nl_diam <- as.data.frame(PC1_nl_diam)
PC1_lin_diam$PC1 <- PC1_nl_diam$PC1 <- PC1_seq
PC1_lin_diam$Response <- PC1_nl_diam$Response <- "Log_Diameter"
PC1_lin_diam$model <- "linear" 
PC1_nl_diam$model <- "nonlinear"

PC1_lin_chl <- PC1_lin_ice$mn[,2,c(10, 1, 5, 7)]
PC1_nl_chl <- PC1_nl_ice$mn[,2,c(8, 3, 7, 5)]
colnames(PC1_lin_chl) <- c("pro", "syn", "pico1", "pico2")
colnames(PC1_nl_chl) <- c("pro", "syn", "pico1", "pico2")
PC1_lin_chl <- as.data.frame(PC1_lin_chl)
PC1_nl_chl <- as.data.frame(PC1_nl_chl)
PC1_lin_chl$PC1 <- PC1_nl_chl$PC1 <- PC1_seq
PC1_lin_chl$Response <- PC1_nl_chl$Response <- "Log_chl"
PC1_lin_chl$model <- "linear" 
PC1_nl_chl$model <- "nonlinear"

PC1_lin_pe <- PC1_lin_ice$mn[,3,c(10, 1, 5, 7)]
PC1_nl_pe <- PC1_nl_ice$mn[,3,c(8, 3, 7, 5)]
colnames(PC1_lin_pe) <- c("pro", "syn", "pico1", "pico2")
colnames(PC1_nl_pe) <- c("pro", "syn", "pico1", "pico2")
PC1_lin_pe <- as.data.frame(PC1_lin_pe)
PC1_nl_pe <- as.data.frame(PC1_nl_pe)
PC1_lin_pe$PC1 <- PC1_nl_pe$PC1 <- PC1_seq
PC1_lin_pe$Response <- PC1_nl_pe$Response <- "Log_pe"
PC1_lin_pe$model <- "linear" 
PC1_nl_pe$model <- "nonlinear"

PC1_pred <- rbind(PC1_lin_probs, PC1_nl_probs, 
                  PC1_lin_diam, PC1_nl_diam, 
                  PC1_lin_chl, PC1_nl_chl, 
                  PC1_lin_pe, PC1_nl_pe)

PC1_pred_melt <- melt(PC1_pred, id.vars = c("PC1", "model", "Response"), 
                                  variable.name = "Population", 
                                  value.name = "Prediction")

PC1_pred_melt$Response <- factor(PC1_pred_melt$Response, 
  levels = c("Relative_Abundance", "Log_Diameter", "Log_chl", "Log_pe"))

PC1_plot <- ggplot(PC1_pred_melt) + 
  geom_line(aes(x = -PC1, y = Prediction, linetype = model, color = model), linewidth = 1) + 
  scale_linetype_manual(values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
  facet_wrap(Population ~ Response, scales = "free") + 
  geom_vline(xintercept = -X_pc[203,1], linetype = "dotted") + 
  geom_vline(xintercept = -X_pc[221,1], linetype = "dotted") +
  # annotate("text", x = -6, y = 5.7, label = "31.9°N", size = 8) +
  # annotate("text", x = 1.25, y = 5.7, label = "34°N", size = 8) + 
  theme(text = element_text(size = 10))

# code to output a list of these plots
PC1_ice_list <- lapply(c("pro", "syn", "pico1", "pico2"), function(pop) {
  lapply(c("Relative_Abundance", "Log_Diameter", "Log_chl", "Log_pe"), function(resp) {
    ggplot(filter(PC1_pred_melt, Population == pop, Response == resp)) + 
  geom_line(aes(x = -PC1, y = Prediction, linetype = model, color = model), linewidth = 1) + 
  scale_linetype_manual(name = "Model", labels = c("Linear", "Nonlinear"), values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
  scale_color_discrete(name = "Model", labels = c("Linear", "Nonlinear")) +
  # facet_wrap(Population ~ Response, scales = "free") + 
  geom_vline(xintercept = -X_pc[203,1], linetype = "dotted") + 
  geom_vline(xintercept = -X_pc[221,1], linetype = "dotted") +
  # annotate("text", x = -6, y = 5.7, label = "31.9°N", size = 8) +
  # annotate("text", x = 1.25, y = 5.7, label = "34°N", size = 8) + 
  theme(text = element_text(size = 22)) + 
  labs(title = paste(pop, resp))
  })
}) %>% unlist(recursive = FALSE)

# FRANCOIS PLOTS

# PC1

# syn chl - "roughly agrees with shape and scale"
pdf("paper_plots/PC1_syn_chl.pdf", 10.5, 9.36)
PC1_ice_list[[7]] + 
  labs(title  = "syn chl")
graphics.off()

# pro diam and pe
pro1 <- PC1_ice_list[[2]] + 
  labs(title  = "pro Log Diameter", x = "") 

leg <- ggpubr::get_legend(pro1)
leg <- ggpubr::as_ggplot(leg)

pro1 <- pro1 + theme(legend.position = "")

pro2 <- PC1_ice_list[[4]] + 
  labs(title  = "pro pe", x = "", y = "") + 
  theme(legend.position = "")

# pro_title <- ggpubr::text_grob("Pro", size = 26, hjust = 2.1)
pro_x_lab <- ggpubr::text_grob("-PC1", size = 22, hjust = 1.25, vjust = -0.5)
pdf("paper_plots/PC1_pro_diam_pe.pdf", 21, 9.36)
grid.arrange(pro1, pro2, leg,  pro_x_lab, nrow = 2, widths = c(1, 1, 0.25), 
              layout_matrix = matrix(c(1, 2, 3, 
                                  4, 4, 4), byrow = TRUE, nrow = 2), heights = c(1, 0.05))
graphics.off()

# pico1 diam & chl
pico1 <- PC1_ice_list[[10]] + 
  labs(title  = "Pico1 Log Diameter", x = "") + 
  theme(legend.position = "")

pico2 <- PC1_ice_list[[11]] + 
  labs(title  = "pico1 pro chl", x = "", y = "") + 
  theme(legend.position = "")

pdf("paper_plots/PC1_pico1_diam_chl.pdf", 21, 9.36)
grid.arrange(pico1, pico2, leg,  pro_x_lab, nrow = 2, widths = c(1, 1, 0.25), 
              layout_matrix = matrix(c(1, 2, 3, 
                                  4, 4, 4), byrow = TRUE, nrow = 2), heights = c(1, 0.05))
graphics.off()

# pro pe, pico1 diam, pico2 pi
pro2 <- PC1_ice_list[[4]] + 
  labs(title  = "Pro Phycoerythrin", x = "", y = "") +
  labs(y = "Prediction")
  
leg <- ggpubr::get_legend(pro2)
leg <- ggpubr::as_ggplot(leg)

pro2 <- pro2 + 
  theme(title = element_text(size = 20), 
        axis.title.y = element_text(size = 18), 
        legend.text = element_text(size = 18), legend.title = element_text(size = 20), 
        legend.position = "none") + 
  scale_x_continuous(breaks = seq(-8, 8, by = 4)) + 
  scale_y_continuous(breaks = seq(0.31, 0.39, by = 0.02), limits = c(0.31, 0.39))


pico1 <- pico1 + 
  labs(y = "") + 
  theme(legend.position = "none", title = element_text(size = 20), 
        axis.title.y = element_text(size = 18), 
        legend.text = element_text(size = 18), legend.title = element_text(size = 20)) + 
  scale_x_continuous(breaks = seq(-8, 8, by = 4)) + 
  scale_y_continuous(breaks = seq(5.4, 5.8, 0.1), limits = c(5.4, 5.8))

pico2 <- PC1_ice_list[[13]] + 
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.04), limits = c(0, 0.2)) + 
  scale_x_continuous(breaks = seq(-8, 8, by = 4)) + 
  labs(title = "Pico2 Relative Abundance", y = "", x = "") + 
  theme(legend.position = "none", title = element_text(size = 20), 
        axis.title.y = element_text(size = 18), 
        legend.text = element_text(size = 18), legend.title = element_text(size = 20))

pro_x_lab <- ggpubr::text_grob("-PC1 (Proxy for Latitude)", size = 18, hjust = 0.375, vjust = -0.5)
pdf("paper_plots/PC1_paper_plot.pdf", 20.18, 6.65)
grid.arrange(pro2, pico1, pico2, leg, pro_x_lab, nrow = 2, 
            layout_matrix = matrix(c(1, 2, 3, 4, 
                                     NA, 5, NA, NA), 2, byrow = TRUE), 
            widths = c(1, 1, 1, 0.25), heights = c(1, 0.05))
graphics.off()

# pico2 relative abundance, "optimal latitude for abundance"
pdf("paper_plots/PC1_pico2_prob.pdf", 10.5, 9.36)
PC1_ice_list[[13]] + 
  geom_vline(xintercept = -X_pc[32, 1], linetype = "dashed") + 
  annotate("text", x = 1.9, y = 0.145, label = "38.8°N", size = 7) + 
  # geom_hline(yintercept = 0.167, linetype = "dashed") + 
  geom_vline(xintercept = 4.25, linetype = "dashed") + 
  annotate("text", x = 5.1, y = 0.145, label = "40°N", size = 7) + 
  # geom_hline(yintercept = 0.162, linetype = "dashed") + 
  scale_y_continuous(breaks = seq(0, 0.18, by = 0.02), limits = c(0, 0.18)) + 
  labs(title = "pico2 Relative Abundance")
graphics.off()

# other PCs

plot_ice <- function(PC_lin_ice, PC_nl_ice, PC_num, PC_seq, return_plot_list = FALSE) {
  PC_lin_probs <- PC_lin_ice$prob[,c(10, 1, 5, 7)]
  PC_nl_probs <- PC_nl_ice$prob[,c(8, 3, 7, 5)]
  colnames(PC_lin_probs) <- c("pro", "syn", "pico1", "pico2")
  colnames(PC_nl_probs) <- c("pro", "syn", "pico1", "pico2")
  PC_lin_probs <- as.data.frame(PC_lin_probs)
  PC_nl_probs <- as.data.frame(PC_nl_probs)
  PC_lin_probs$PC1 <- PC_nl_probs$PC1 <- factor(other_ice_X[,1], levels = c(-4, 1, 6), 
                                                labels = c("Subpolar", 
                                                           "Transition Zone", 
                                                           "Subtropical"))
  PC_lin_probs[,PC_num] <- PC_nl_probs[,PC_num] <- PC_seq
  PC_lin_probs$Response <- PC_nl_probs$Response <- "Relative_Abundance"
  PC_lin_probs$model <- "linear" 
  PC_nl_probs$model <- "nonlinear"

  PC_lin_diam <- PC_lin_ice$mn[,1,c(10, 1, 5, 7)]
  PC_nl_diam <- PC_nl_ice$mn[,1,c(8, 3, 7, 5)]
  colnames(PC_lin_diam) <- c("pro", "syn", "pico1", "pico2")
  colnames(PC_nl_diam) <- c("pro", "syn", "pico1", "pico2")
  PC_lin_diam <- as.data.frame(PC_lin_diam)
  PC_nl_diam <- as.data.frame(PC_nl_diam)
  PC_lin_diam$PC1 <- PC_nl_diam$PC1 <- factor(other_ice_X[,1], levels = c(-4, 1, 6), labels = c("Subpolar", 
                                                                              "Transition Zone", 
                                                                              "Subtropical"))
  PC_lin_diam[,PC_num] <- PC_nl_diam[,PC_num] <- PC_seq
  PC_lin_diam$Response <- PC_nl_diam$Response <- "Log_Diameter"
  PC_lin_diam$model <- "linear" 
  PC_nl_diam$model <- "nonlinear"

  PC_lin_chl <- PC_lin_ice$mn[,2,c(10, 1, 5, 7)]
  PC_nl_chl <- PC_nl_ice$mn[,2,c(8, 3, 7, 5)]
  colnames(PC_lin_chl) <- c("pro", "syn", "pico1", "pico2")
  colnames(PC_nl_chl) <- c("pro", "syn", "pico1", "pico2")
  PC_lin_chl <- as.data.frame(PC_lin_chl)
  PC_nl_chl <- as.data.frame(PC_nl_chl)
  PC_lin_chl$PC1 <- PC_nl_chl$PC1 <- factor(other_ice_X[,1], levels = c(-4, 1, 6), labels = c("Subpolar", 
                                                                              "Transition Zone", 
                                                                              "Subtropical"))
  PC_lin_chl[,PC_num] <- PC_nl_chl[,PC_num] <- PC_seq
  PC_lin_chl$Response <- PC_nl_chl$Response <- "Log_chl"
  PC_lin_chl$model <- "linear" 
  PC_nl_chl$model <- "nonlinear"

  PC_lin_pe <- PC_lin_ice$mn[,3,c(10, 1, 5, 7)]
  PC_nl_pe <- PC_nl_ice$mn[,3,c(8, 3, 7, 5)]
  colnames(PC_lin_pe) <- c("pro", "syn", "pico1", "pico2")
  colnames(PC_nl_pe) <- c("pro", "syn", "pico1", "pico2")
  PC_lin_pe <- as.data.frame(PC_lin_pe)
  PC_nl_pe <- as.data.frame(PC_nl_pe)
  PC_lin_pe$PC1 <- PC_nl_pe$PC1 <- factor(other_ice_X[,1], levels = c(-4, 1, 6), labels = c("Subpolar", 
                                                                              "Transition Zone", 
                                                                              "Subtropical"))
  PC_lin_pe[,PC_num] <- PC_nl_pe[,PC_num] <- PC_seq
  PC_lin_pe$Response <- PC_nl_pe$Response <- "Log_pe"
  PC_lin_pe$model <- "linear" 
  PC_nl_pe$model <- "nonlinear"

  PC_pred <- rbind(PC_lin_probs, PC_nl_probs, 
                    PC_lin_diam, PC_nl_diam, 
                    PC_lin_chl, PC_nl_chl, 
                    PC_lin_pe, PC_nl_pe)

  PC_pred_melt <- melt(PC_pred, id.vars = c("PC1", PC_num, "model", "Response"), 
                                    variable.name = "Population", 
                                    value.name = "Prediction")

  PC_pred_melt$Response <- factor(PC_pred_melt$Response, 
    levels = c("Relative_Abundance", "Log_Diameter", "Log_chl", "Log_pe"))

  subsampled_pts <- PC_pred_melt %>% 
      group_split(Population, Response, model, PC1) %>% 
      lapply(function(cur_df) {
        cur_seq <- seq(1, nrow(cur_df), length.out = 6) %>% round()
        cur_seq <- cur_seq[2:5]
        cur_df[cur_seq,]
      }) %>% bind_rows()

  if(return_plot_list) {
    lapply(c("pro", "syn", "pico1", "pico2"), function(pop) {
      lapply(c("Relative_Abundance", "Log_Diameter", "Log_chl", "Log_pe"), function(resp) {
        ggplot(filter(PC_pred_melt, Response == resp, Population == pop)) + 
          geom_line(aes(x = .data[[PC_num]], y = Prediction, linetype = model, color = PC1), linewidth = 1) + 
          geom_point(data = filter(subsampled_pts, Response == resp, Population == pop), aes(x = .data[[PC_num]], y = Prediction, group = model, color = PC1, shape = PC1), size = 4) + 
          scale_linetype_manual(name = "Model", labels = c("Linear", "Nonlinear"), values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
          # facet_wrap(Population ~ Response, scales = "free")
          labs(title = pop, y = resp) + 
          theme(text = element_text(size = 20))
      })
    }) %>% unlist(recursive = FALSE)
  } else { # put four points on lines to distinguish if plots are printed in greyscale
    ggplot(PC_pred_melt) + 
          geom_line(aes(x = .data[[PC_num]], y = Prediction, linetype = model, color = PC1), linewidth = 1) + 
          geom_point(data = subsampled_pts, aes(x = .data[[PC_num]], y = Prediction, group = model, color = PC1, shape = PC1), size = 2.5) + 
          scale_linetype_manual(name = "Model", labels = c("Linear", "Nonlinear"), values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
          facet_wrap(Population ~ Response, scales = "free") + 
          theme(text = element_text(size = 10))
  }
}

# Previous (EDA)
# pdf("paper_plots/PC1_plot.pdf", 10, 8.5)
# PC1_plot
# graphics.off()
# pdf("paper_plots/PC2_plot.pdf", 10, 8.5)
# plot_ice(PC2_lin_ice, PC2_nl_ice, "PC2", PC2_ice_X[,2])
# graphics.off()
# pdf("paper_plots/PC3_plot.pdf", 10, 8.5)
# plot_ice(PC3_lin_ice, PC3_nl_ice, "PC3", PC3_ice_X[,3])
# graphics.off()
# pdf("paper_plots/PC4_plot.pdf", 10, 8.5)
# plot_ice(PC4_lin_ice, PC4_nl_ice, "PC4", PC4_ice_X[,4])
# graphics.off()

# PC2

PC2_plot <- plot_ice(PC2_lin_ice, PC2_nl_ice, "PC2", PC2_seq)
PC2_plist <- plot_ice(PC2_lin_ice, PC2_nl_ice, "PC2", PC2_seq, return_plot_list = TRUE)

# pro abundance
PC2_1 <- PC2_plist[[1]] + 
  labs(title = "Pro Relative Abundance", y = "Prediction") + 
  theme(legend.position = "none", title = element_text(size = 20), 
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 16), 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18)) + 
  labs(x = "", y = "Prediction")

# syn pe
PC2_2 <- PC2_plist[[8]] + 
  labs(title = "Syn Phycoerythrin", y = "Prediction") + 
  theme(legend.position = "none", title = element_text(size = 20), 
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 16), 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18)) + 
  labs(x = "", y = "")

# pico2 abundance
PC2_3 <- PC2_plist[[13]] + 
  labs(title = "Pico2 Relative Abundance", y = "Prediction") + 
  labs(x = "", y = "") + 
  theme(legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18), 
        legend.key.size = unit(1, "cm"))

PC2_leg <- ggpubr::get_legend(PC2_3) %>% 
  ggpubr::as_ggplot()

PC2_3 <- PC2_3 + 
  theme(legend.position = "none", title = element_text(size = 20), 
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 16), 
        )

# pdf("paper_plots/PC2_pro_pi.pdf", 10.5, 9.36)
# pdf("paper_plots/PC2_syn_pe.pdf", 10.5, 9.36)
# pdf("paper_plots/PC2_pico2_pi.pdf", 10.5, 9.36)
PC2_xlab <- ggpubr::text_grob("PC2", hjust = 2.5, vjust = -0.5, size = 18)

pdf("paper_plots/PC2_selected.pdf", 21.3, 6.7)
grid.arrange(PC2_1, PC2_2, PC2_3, PC2_leg, PC2_xlab, 
             layout_matrix = matrix(c(1, 2, 3, 4, 
                                       5, 5, 5, 5), nrow = 2, byrow = TRUE), 
             widths = c(1, 1, 1, 0.4), heights = c(1, 0.05))
graphics.off()
##########

# Start with all PCs fixed at their mean
other_ice_X_3_4 <- tcrossprod(rep(1, 296 * 3), colMeans(X_pc))

# check the ICEs of PC2, PC3, and PC4 with PC1 fixed at values corresponding 
# to low, medium, and high latitude
PC1_levels <- c(-4, 1, 6)
other_ice_X_3_4[,1] <- rep(PC1_levels, each = 296)

PC3_4_ice_X <- other_ice_X_3_4

PC3_4_ice_X[,3] <- rep(X_pc[,3], times = 3)
PC3_4_ice_X[,4] <- rep(X_pc[,4], times = 3)

# par should be taken from original dataset and repeated 3 times

set.seed(0)
PC3_4_ice_X_nl <- X_hidden(PC3_4_ice_X, 9, n.h, 0.5)

PC3_4_lin_ice <- predict(linear_best, newx = PC3_4_ice_X, logits = FALSE)
PC3_4_nl_ice <- predict(half_best, newx = PC3_4_ice_X_nl, logits = FALSE)


plot_ice_2 <- function(PC_lin_ice, PC_nl_ice, PC_num, PC_seq, other_ice_X, return_plot_list = FALSE) {
  PC_lin_probs <- PC_lin_ice$prob[,c(10, 1, 5, 7)]
  PC_nl_probs <- PC_nl_ice$prob[,c(8, 3, 7, 5)]
  colnames(PC_lin_probs) <- c("pro", "syn", "pico1", "pico2")
  colnames(PC_nl_probs) <- c("pro", "syn", "pico1", "pico2")
  PC_lin_probs <- as.data.frame(PC_lin_probs)
  PC_nl_probs <- as.data.frame(PC_nl_probs)
  PC_lin_probs$PC1 <- PC_nl_probs$PC1 <- factor(other_ice_X[,1], levels = c(-4, 1, 6), 
                                                labels = c("Subpolar", 
                                                           "Transition Zone", 
                                                           "Subtropical"))
  PC_lin_probs[,PC_num] <- PC_nl_probs[,PC_num] <- PC_seq
  PC_lin_probs$Response <- PC_nl_probs$Response <- "Relative_Abundance"
  PC_lin_probs$model <- "linear" 
  PC_nl_probs$model <- "nonlinear"

  PC_lin_diam <- PC_lin_ice$mn[,1,c(10, 1, 5, 7)]
  PC_nl_diam <- PC_nl_ice$mn[,1,c(8, 3, 7, 5)]
  colnames(PC_lin_diam) <- c("pro", "syn", "pico1", "pico2")
  colnames(PC_nl_diam) <- c("pro", "syn", "pico1", "pico2")
  PC_lin_diam <- as.data.frame(PC_lin_diam)
  PC_nl_diam <- as.data.frame(PC_nl_diam)
  PC_lin_diam$PC1 <- PC_nl_diam$PC1 <- factor(other_ice_X[,1], levels = c(-4, 1, 6), labels = c("Subpolar", 
                                                                              "Transition Zone", 
                                                                              "Subtropical"))
  PC_lin_diam[,PC_num] <- PC_nl_diam[,PC_num] <- PC_seq
  PC_lin_diam$Response <- PC_nl_diam$Response <- "Log_Diameter"
  PC_lin_diam$model <- "linear" 
  PC_nl_diam$model <- "nonlinear"

  PC_lin_chl <- PC_lin_ice$mn[,2,c(10, 1, 5, 7)]
  PC_nl_chl <- PC_nl_ice$mn[,2,c(8, 3, 7, 5)]
  colnames(PC_lin_chl) <- c("pro", "syn", "pico1", "pico2")
  colnames(PC_nl_chl) <- c("pro", "syn", "pico1", "pico2")
  PC_lin_chl <- as.data.frame(PC_lin_chl)
  PC_nl_chl <- as.data.frame(PC_nl_chl)
  PC_lin_chl$PC1 <- PC_nl_chl$PC1 <- factor(other_ice_X[,1], levels = c(-4, 1, 6), labels = c("Subpolar", 
                                                                              "Transition Zone", 
                                                                              "Subtropical"))
  PC_lin_chl[,PC_num] <- PC_nl_chl[,PC_num] <- PC_seq
  PC_lin_chl$Response <- PC_nl_chl$Response <- "Log_chl"
  PC_lin_chl$model <- "linear" 
  PC_nl_chl$model <- "nonlinear"

  PC_lin_pe <- PC_lin_ice$mn[,3,c(10, 1, 5, 7)]
  PC_nl_pe <- PC_nl_ice$mn[,3,c(8, 3, 7, 5)]
  colnames(PC_lin_pe) <- c("pro", "syn", "pico1", "pico2")
  colnames(PC_nl_pe) <- c("pro", "syn", "pico1", "pico2")
  PC_lin_pe <- as.data.frame(PC_lin_pe)
  PC_nl_pe <- as.data.frame(PC_nl_pe)
  PC_lin_pe$PC1 <- PC_nl_pe$PC1 <- factor(other_ice_X[,1], levels = c(-4, 1, 6), labels = c("Subpolar", 
                                                                              "Transition Zone", 
                                                                              "Subtropical"))
  PC_lin_pe[,PC_num] <- PC_nl_pe[,PC_num] <- PC_seq
  PC_lin_pe$Response <- PC_nl_pe$Response <- "Log_pe"
  PC_lin_pe$model <- "linear" 
  PC_nl_pe$model <- "nonlinear"

  PC_pred <- rbind(PC_lin_probs, PC_nl_probs, 
                    PC_lin_diam, PC_nl_diam, 
                    PC_lin_chl, PC_nl_chl, 
                    PC_lin_pe, PC_nl_pe)

  PC_pred_melt <- melt(PC_pred, id.vars = c("PC1", PC_num, "model", "Response"), 
                                    variable.name = "Population", 
                                    value.name = "Prediction")

  PC_pred_melt$Response <- factor(PC_pred_melt$Response, 
    levels = c("Relative_Abundance", "Log_Diameter", "Log_chl", "Log_pe"))

  if(return_plot_list) {
    lapply(c("pro", "syn", "pico1", "pico2"), function(pop) {
      lapply(c("Relative_Abundance", "Log_Diameter", "Log_chl", "Log_pe"), function(resp) {
        if(resp == "Relative_Abundance") {
          ggplot(filter(PC_pred_melt, Response == resp, Population == pop)) + 
            geom_smooth(aes(x = .data[[PC_num]], y = Prediction, linetype = model, color = PC1), se = FALSE, linewidth = 1) + 
            scale_linetype_manual(name = "Model", labels = c("Linear", "Nonlinear"), values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
            labs(title = paste0(pop, resp)) + 
            theme(text = element_text(size = 22))
        } else {
          ggplot(filter(PC_pred_melt, model == "nonlinear", Response == resp, Population == pop)) + 
            geom_smooth(aes(x = .data[[PC_num]], y = Prediction, linetype = model, color = PC1), se = FALSE, linewidth = 1) + 
            geom_smooth(data = filter(PC_pred_melt, model == "linear", Response == resp, Population == pop), aes(x = .data[[PC_num]], y = Prediction, linetype = model, color = PC1), se = FALSE, method = "lm", linewidth = 1) + 
            scale_linetype_manual(values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
            labs(title = paste0(pop, resp)) + 
            theme(text = element_text(size = 22))
        }
      })
    }) %>% unlist(recursive = FALSE)
  } else {
    # using lm to straighten the line model predictions, 
    # except for the relative abundance plots
    # TODO: Look up stat_smooth + after_stat, but if it doesn't work out, just use lm() 
    # and gather fitted values in a df
    p1 <- ggplot(filter(PC_pred_melt, model == "nonlinear", Response == "Relative_Abundance")) + 
            geom_smooth(aes(x = .data[[PC_num]], y = Prediction, linetype = model, color = PC1), se = FALSE, linewidth = 0.4) + 
            geom_smooth(data = filter(PC_pred_melt, model == "linear", Response == "Relative_Abundance"), aes(x = .data[[PC_num]], y = Prediction, linetype = model, color = PC1), se = FALSE, linewidth = 0.4) + 
            scale_linetype_manual(values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
            facet_wrap(Population ~ Response, nrow = 4, scales = "free") + 
            theme(text = element_text(size = 10)) + 
            labs(x = "")

    leg <- ggpubr::get_legend(p1) %>% 
      ggpubr::as_ggplot()
    
    p1 <- p1 + theme(legend.position = "none")

    p2 <- ggplot(filter(PC_pred_melt, model == "nonlinear", Response != "Relative_Abundance")) + 
            geom_smooth(aes(x = .data[[PC_num]], y = Prediction, linetype = model, color = PC1), se = FALSE, linewidth = 0.4) + 
            geom_smooth(data = filter(PC_pred_melt, model == "linear", Response != "Relative_Abundance"), aes(x = .data[[PC_num]], y = Prediction, linetype = model, color = PC1), se = FALSE, linewidth = 0.4, method = "lm") + 
            scale_linetype_manual(values = c("linear" = "dotdash", "nonlinear" = "solid")) + 
            facet_wrap(Population ~ Response, nrow = 4, scales = "free") + 
            theme(text = element_text(size = 10), legend.position = "none") + 
            labs(y = "", x = "")
    
    x_lab <- ggpubr::text_grob("PC4", vjust = -1, hjust = 2)
    
    grid.arrange(p1, p2, leg, x_lab, layout_matrix = matrix(c(1, 2, 3, 
                                                              4,4, 4), 
                                                            byrow = TRUE, nrow = 2), 
                widths = c(1, 2.75, 0.5), heights = c(1, 0.05))
  }
}

par_seq <- rep(X[,"par"], times = 3)

plot_ice_2(PC3_4_lin_ice, PC3_4_nl_ice, "PC4", PC3_4_ice_X[,4], other_ice_X_3_4)
plot_ice_2(PC3_4_lin_ice, PC3_4_nl_ice, "par", par_seq, other_ice_X_3_4)
plist_2 <- plot_ice_2(PC3_4_lin_ice, PC3_4_nl_ice, "PC4", PC3_4_ice_X[,4], other_ice_X_3_4, return_plot_list = TRUE)

# Example datetime
dt <- as.POSIXct(time, format = "%Y-%m-%dT%T", )

# Extract hour and minute
hm <- format(dt, "%H:%M")
hm <- as.POSIXct(hm, format = "%H:%M")

ord <- order(hm)
hm_sort <- hm[ord]
values <- linear_best$mn[ord, 1, 10]
values <- half_best$mn[ord, 1, 8]

# TODO: work on this
# as an initial exercise, color the points by PC1

# PC_lin_ice$prob[,c(10, 1, 5, 7)]
# PC_nl_probs <- PC_nl_ice$prob[,c(8, 3, 7, 5)]

diam_time <- data.frame(Time_of_Day = rep(hm_sort, 8), 
                        Log_Diameter = c(linear_best$mn[ord, 1, 10], half_best$mn[ord,1,8], 
                                         linear_best$mn[ord, 1, 1], half_best$mn[ord,1,3], 
                                         linear_best$mn[ord, 1, 5], half_best$mn[ord,1,7], 
                                         linear_best$mn[ord, 1, 7], half_best$mn[ord,1,5]),
                        Model = rep(rep(c("Linear", "Nonlinear"), each = 296), 4), 
                        Cluster = rep(c("pro", "syn", "pico1", "pico2"), each = 296 * 2),
                        PC1 = rep(X_pc[,1], 8), 
                        PC2 = rep(X_pc[,2], 8), 
                        PC3 = rep(X_pc[,3], 8), 
                        PC4 = rep(X_pc[,4], 8))

pdf("paper_plots/diam_time_of_day.pdf", 19.4, 9.2)
ggplot(diam_time) + 
  geom_point(aes(Time_of_Day, Log_Diameter)) + 
  facet_wrap(Cluster ~ Model, scales = "free", nrow = 2)
graphics.off()

# TODO: figure out what explains the group structure

# Plot with formatted time labels
par(mfrow = c(1, 2))
plot(hm_sort, values, type = "p", xaxt = "n", xlab = "Time of Day", ylab = "Value")
axis(1, at = hm_sort, labels = format(hm_sort, "%H:%M"), cex.axis = 0.8)

par <- X[ord,"par"]
plot(hm_sort, par, type = "p", xaxt = "n", xlab = "Time of Day", ylab = "Value")
axis(1, at = hm_sort, labels = format(hm_sort, "%H:%M"), cex.axis = 0.8)

# plot(X[,"par"], linear_best$prob[,1])
# plot(X[,"par"], half_best$prob[,3])

# interesting ones
# pdf("paper_plots/PC4_pro_pi.pdf", 10.5, 9.36)
# pdf("paper_plots/PC4_pico2_pi.pdf", 10.5, 9.36)
# pdf("paper_plots/PC4_pro_diam.pdf", 10.5, 9.36)
# pdf("paper_plots/PC4_syn_diam.pdf", 10.5, 9.36)
# pdf("paper_plots/PC4_pico1_diam.pdf", 10.5, 9.36)
# pdf("paper_plots/PC4_pico2_diam.pdf", 10.5, 9.36)
# pdf("paper_plots/PC4_syn_chl.pdf", 10.5, 9.36)
# pdf("paper_plots/PC4_pico1_chl.pdf", 10.5, 9.36)
# pdf("paper_plots/PC4_pico1_pe.pdf", 10.5, 9.36)

# pro relative abundance
PC4_pro_1 <- plist_2[[1]] + 
  theme(legend.position = "none") + 
  labs(title = "pro Relative Abundance", x = "", y = "Prediction")

# pro diam
PC4_pro_2 <- plist_2[[2]] + 
  labs(title = "pro Log Diameter", x = "", y = "")

PC4_leg <- ggpubr::get_legend(PC4_pro_2) %>% 
  ggpubr::as_ggplot()

PC4_pro_2 <- PC4_pro_2 + 
  theme(legend.position = "none")

PC4_xlab_2 <- ggpubr::text_grob("PC4", hjust = 2.25, vjust = -0.5, size = 18)
pdf("paper_plots/PC4_pro.pdf", 15.75, 7)
grid.arrange(PC4_pro_1, PC4_pro_2, PC4_leg, PC4_xlab_2, 
             layout_matrix = matrix(c(1, 2, 3, 
                                      4, 4, 4), nrow = 2, byrow = TRUE), 
              widths = c(1, 1, 0.4), heights = c(1, 0.05))
graphics.off()

# syn diam
PC4_syn1 <- plist_2[[6]] + 
  theme(legend.position = "none") + 
  labs(title = "syn Log Diameter", x = "", y = "Prediction")

# syn chl 
PC4_syn2 <- plist_2[[7]] + 
  theme(legend.position = "none") + 
  labs(title = "syn chl", x = "", y = "")

pdf("paper_plots/PC4_syn.pdf", 15.75, 7)
grid.arrange(PC4_syn1, PC4_syn2, PC4_leg, PC4_xlab_2, 
             layout_matrix = matrix(c(1, 2, 3, 
                                      4, 4, 4), nrow = 2, byrow = TRUE), 
              widths = c(1, 1, 0.4), heights = c(1, 0.05))
graphics.off()

# pico1 diam
PC4_pico1_1 <- plist_2[[10]] + 
  theme(legend.position = "none") + 
  labs(title = "pico1 Log Diameter", x = "", y = "Prediction")

# pico1 chl 
PC4_pico1_2 <- plist_2[[11]] + 
  theme(legend.position = "none") + 
  labs(title = "pico1 chl", x = "", y = "")

# pico1 pe
PC4_pico1_3 <- plist_2[[12]] + 
  theme(legend.position = "none") + 
  labs(title = "pico1 pe", x = "", y = "")

PC4_xlab_3 <- ggpubr::text_grob("PC4", hjust = 2.25, vjust = -0.5, size = 18)
pdf("paper_plots/PC4_pico1.pdf", 18, 5.75)
grid.arrange(PC4_pico1_1, PC4_pico1_2, PC4_pico1_3, PC4_leg, 
            PC4_xlab_3, 
             layout_matrix = matrix(c(1, 2, 3, 4, 
                                      5, 5, 5, 5), nrow = 2, byrow = TRUE), 
              widths = c(1, 1, 1, 0.6), heights = c(1, 0.05))
graphics.off()

# pico2 relative abundance
PC4_pico2_1 <- plist_2[[13]] + 
  theme(legend.position = "none") + 
  labs(title = "pico2 Relative Abundance", x = "", y = "Prediction")

# pico2 diam 
PC4_pico2_2 <- plist_2[[14]] + 
  theme(legend.position = "none") + 
  labs(title = "pico2 Log Diameter", x = "", y = "")

pdf("paper_plots/PC4_pico2.pdf", 15.75, 7)
grid.arrange(PC4_pico2_1, PC4_pico2_2, PC4_leg, PC4_xlab_2, 
             layout_matrix = matrix(c(1, 2, 3, 
                                      4, 4, 4), nrow = 2, byrow = TRUE), 
              widths = c(1, 1, 0.4), heights = c(1, 0.05))
graphics.off()


############

# TODO: compute 2d KDEs and make predictions
my_kde2d <- function (x, y, h, gx, gy, lims = c(range(x), range(y))) {
    nx <- length(x)
    if (length(y) != nx) 
        stop("data vectors must be the same length")
    if (any(!is.finite(x)) || any(!is.finite(y))) 
        stop("missing or infinite values in the data are not allowed")
    if (any(!is.finite(lims))) 
        stop("only finite values are allowed in 'lims'")
    # n <- rep(n, length.out = 2L)
    # gx <- seq.int(lims[1L], lims[2L], length.out = n[1L])
    # gy <- seq.int(lims[3L], lims[4L], length.out = n[2L])
    h <- if (missing(h)) 
        c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))
    else rep(h, length.out = 2L)
    if (any(h <= 0)) 
        stop("bandwidths must be strictly positive")
    h <- h/4
    ax <- outer(gx, x, "-")/h[1L]
    ay <- outer(gy, y, "-")/h[2L]
    z <- tcrossprod(matrix(dnorm(ax), , nx), matrix(dnorm(ay), 
        , nx))/(nx * h[1L] * h[2L])
    list(x = gx, y = gy, z = z)
}


# Start with all PCs fixed at their mean
ice_X_2d <- tcrossprod(rep(1, 30**2 * 3), colMeans(X_pc))

# low, medium, and high latitude
PC1_levels <- c(-4, 1, 6)
ice_X_2d[,1] <- rep(PC1_levels, each = 30**2)

PC3_4_grid <- expand.grid(PC3_seq, PC4_seq)

ice_X_2d[,c(3, 4)] <- kronecker(rep(1, 3), as.matrix(PC3_4_grid))

set.seed(0)
ice_X_2d_nl <- X_hidden(ice_X_2d, 9, n.h, 0.5)

lin_ice_2d <- predict(linear_best, newx = ice_X_2d, logits = FALSE)
nl_ice_2d <- predict(half_best, newx = ice_X_2d_nl, logits = FALSE)

kde <- MASS::kde2d(PC3_seq, PC4_seq, n = 30)

######

# TODO: 2d ICE

PC_lin_ice <- lin_ice_2d 
PC_nl_ice <- nl_ice_2d

PC_lin_probs <- PC_lin_ice$prob[,c(10, 1, 5, 7)]
PC_nl_probs <- PC_nl_ice$prob[,c(8, 3, 7, 5)]
colnames(PC_lin_probs) <- c("pro", "syn", "pico1", "pico2")
colnames(PC_nl_probs) <- c("pro", "syn", "pico1", "pico2")
PC_lin_probs <- as.data.frame(PC_lin_probs)
PC_nl_probs <- as.data.frame(PC_nl_probs)
PC_lin_probs$PC1 <- PC_nl_probs$PC1 <- factor(ice_X_2d[,1], levels = c(-4, 1, 6), 
                                              labels = c("Subpolar", 
                                                         "Transition Zone", 
                                                         "Subtropical"))

PC_lin_probs$PC3 <- PC_nl_probs$PC3 <- ice_X_2d[,3]
PC_lin_probs$PC4 <- PC_nl_probs$PC4 <- ice_X_2d[,4]

PC_lin_probs$Response <- PC_nl_probs$Response <- "Relative_Abundance"
PC_lin_probs$model <- "linear" 
PC_nl_probs$model <- "nonlinear"

PC_pred <- rbind(PC_lin_probs, PC_nl_probs)

PC_pred_melt <- melt(PC_pred, id.vars = c("PC1", "PC3", "PC4", "model", "Response"), 
                                    variable.name = "Population", 
                                    value.name = "Prediction")

  PC_pred_melt$Response <- factor(PC_pred_melt$Response, 
    levels = c("Relative_Abundance", "Log_Diameter", "Log_chl", "Log_pe"))

ggplot(filter(PC_pred_melt)) + 
        geom_contour_filled(aes(x = PC3, y = PC4, z = Prediction)) + 
        facet_wrap(~ model + PC1 + Population, scales = "free") + 
        theme(text = element_text(size = 10)) 
head(PC_pred_melt)
