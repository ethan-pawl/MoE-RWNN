library(myUtils)
library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(flowmix)
library(tibble)
library(ggpubr)
library(corrplot)

simdata_dir <- file.path("~", 
                         "00_Cyto", 
                         "data",
                         "simdata", 
                         "simdata_04_01")

plot_dir <- "plots"
if(!dir.exists(plot_dir)) dir.create(plot_dir)

cvres_files <- list.files("results", recursive = TRUE)

datobj <- 

###################################

# Helper Functions

get_simdata <- function(fname) {
    simdata_fname <- fname %>% 
        strsplit("-") %>%
        unlist() %>% 
        .[1:3] %>% 
        paste(collapse = "-") %>% 
        paste0("simdata-", ., "-1.Rdata") %>%
        file.path(simdata_dir, .)
    
    load(simdata_fname)
    res
}

get_oos <- function(fname) {
    simdata_fname <- fname %>% 
        strsplit("-") %>%
        unlist() %>% 
        .[1:3] %>% 
        paste(collapse = "-") %>% 
        paste0("simdata-", ., "-0.Rdata") %>%
        file.path(simdata_dir, .)
    
    load(simdata_fname)
    res
}

get_title <- function(cv, plot_type) {
    params <- cv$destin %>% 
        strsplit("/") %>% 
        unlist() %>%
        .[2] %>% 
        strsplit("-") %>% 
        unlist()
    
    data_config <- switch(paste(params[1], params[2], sep = "-"), 
                          `1-1` = "Linear Data", 
                          `2-1` = "Interaction Mean", 
                          `3-1` = "Quadratic Mean", 
                          `4-1` = "Logistic Mean", 
                          `1-2` = "Interaction Logit", 
                          `1-3` = "Quadratic Logit", 
                          `1-4` = "Logistic Logit")

    model_fit <- switch(params[4], 
                        `l` = "Linear", 
                        `n` = "Nonlinear")

    label <- switch(plot_type, 
                    mn = "Cluster Means", 
                    prob = "Cluster Probabilities")

    paste0(label, ": ", data_config, ", Signal Size", params[3], 
           ", ", model_fit, " Model Fit")
}

library(RColorBrewer)
display.brewer.pal(n = 8, name = "Set1")


# Extract colors from the "Paired" palette
brew_colors <- brewer.pal(n = 8, name = "Set1")

# Assign specific colors to each category (manually choosing 4 for example)
named_colors <- c("1" = brew_colors[1],
                  "2" = brew_colors[2],
                  "1, Truth" = brew_colors[4],
                  "2, Truth" = brew_colors[5])

# Labels for the legend
labels <- c("1" = "1, Estimated",
            "2" = "2, Estimated",
            "1, Truth" = "1, Truth",
            "2, Truth" = "2, Truth")

plot_sim_means <- function(simdat, cv, plot_band = TRUE, plot_model = TRUE) {
    if(plot_model) {
        flowtrend::plot_1d(simdat$ybin_list, simdat$countslist, cv$bestres, 
                           bin = TRUE, plot_band = plot_band) + 
            geom_line(data = data.frame(x = rep(1:296, 2), 
                                        y = c(simdat$mean_spec, simdat$pico_mu), 
                                        prob = c(simdat$prob_spec, 1 - simdat$prob_spec),
                                        clust = rep(c("1, Truth", "2, Truth"), 
                                                    each = 296)), 
                      mapping = aes(x = x, y = y, color = clust), lineend = "round", linewidth = 1, 
                      linetype = "twodash") + 
            scale_color_manual(name = "Cluster", values = named_colors, labels = labels) + 
            # scale_color_manual(name = "Cluster", 
            #                    values = c("1" = "red", 
            #                               "2" = "blue", 
            #                               "1, Truth" = "purple", 
            #                               "2, Truth" = "#E69F00"), 
            #                    labels = c("1" = "1, Estimated", 
            #                               "2" = "2, Estimated", 
            #                               "1, Truth" = "1, Truth", 
            #                               "2, Truth" = "2, Truth")) + 
            scale_linewidth_continuous(name = "Cluster Probability", range = c(0.25, 4)) + 
            labs(title = get_title(cv, 1))
    } else {
        flowtrend::plot_1d(simdat$ybin_list, simdat$countslist, 
                           bin = TRUE) + 
            geom_line(data = data.frame(x = rep(1:296, 2), 
                                        y = c(simdat$mean_spec, simdat$pico_mu), 
                                        prob = c(simdat$prob_spec, 1 - simdat$prob_spec),
                                        clust = rep(c("1, Truth", "2, Truth"), 
                                                    each = 296)), 
                      mapping = aes(x = x, y = y, color = clust, linewidth = prob), lineend = "round") + 
            scale_color_manual(name = "Cluster", 
                               values = c("1, Truth" = "purple", 
                                          "2, Truth" = "#E69F00")) + 
            scale_linewidth_continuous(name = "Cluster Probability", 
                                       range = c(0.25, 4)) + 
            labs(title = get_title(cv, 1))
    } 
}

plot_sim_probs <- function(simdat, cv, plot_band = TRUE) {
    flowtrend::plot_prob(cv$bestres) + 
        geom_line(data = data.frame(x = rep(1:296, 2), 
                                    y = c(simdat$prob_spec, 
                                          1 - simdat$prob_spec), 
                                    clust = rep(c("1, Truth", "2, Truth"), 
                                                each = 296)), 
                  mapping = aes(x = x, y = y, color = clust)) + 
        scale_color_manual(values = c("1" = "red", 
                                      "2" = "blue", 
                                      "1, Truth" = "purple", 
                                      "2, Truth" = "#E69F00"), 
                           labels = c("1" = "1, Estimated", 
                               "2" = "2, Estimated", 
                               "1, Truth" = "1, Truth", 
                               "2, Truth" = "2, Truth")) + 
        labs(title = get_title(cv, 2))
}

#####################

# Gather results from every simulation into one data frame

results_list <- lapply(cvres_files, function(cv_file) {
    cv <- readRDS(file.path("results", cv_file))
    is_data <- get_simdata(cv_file)
    oos_data <- get_oos(cv_file)

    # fix label switching
    clust_1_pos <- cv$bestres$mn[,1,] %>% 
        apply(2, min) %>%
        which.min()

    if(clust_1_pos == 2) {
        cv$bestres$alpha <- cv$bestres$alpha[2:1,1:10]
        cv$bestres$beta <- cv$bestres$beta[2:1]
        cv$bestres$mn <- cv$bestres$mn[,,2:1, drop = FALSE]
        cv$bestres$prob <- cv$bestres$prob[,2:1, drop = FALSE]
        cv$bestres$sigma <- cv$bestres$sigma[2:1, 1, 1, drop = FALSE]
    }

    # plot clustering over time
    short_name <- substring(cv_file, 1, 7)

    pdf(file.path(plot_dir, paste0(short_name, "_means.pdf")), 9.7, 8.1)
    print(plot_sim_means(is_data, cv))
    graphics.off()

    pdf(file.path(plot_dir, paste0(short_name, "_probs.pdf")), 9.7, 8.1)
    print(plot_sim_probs(is_data, cv))
    graphics.off()

    # log likelihoods
    # I need to collect all the NLLs in a data frame
    oos_nll <- objective_newdat(cv$bestres, oos_data$ybin_list, oos_data$countslist)
    oracle <- list(mn = array(c(is_data$mean_spec, is_data$pico_mu), dim = c(296, 1, 2)), 
                   prob = array(c(is_data$prob_spec, 1 - is_data$prob_spec), dim = c(296, 2)), 
                   sigma = array(is_data$clust_sig^2, dim = c(2, 1, 1)))

    oracle_oos_nll <- objective_newdat(oracle, oos_data$ybin_list, oos_data$countslist)               

    nll_oos_oracle_diff <- oos_nll - oracle_oos_nll

    params <- cv$destin %>% 
            strsplit("/") %>% 
            unlist() %>%
            .[2] %>% 
            strsplit("-") %>% 
            unlist()

    data_config <- switch(paste(params[1], params[2], sep = "-"), 
                          `1-1` = "Linear Data", 
                          `2-1` = "Interaction Mean", 
                          `3-1` = "Quadratic Mean", 
                          `4-1` = "Logistic Mean", 
                          `1-2` = "Interaction Logit", 
                          `1-3` = "Quadratic Logit", 
                          `1-4` = "Logistic Logit")

    model_fit <- switch(params[4], 
                        `l` = "Linear", 
                        `n` = "Nonlinear")

    results <- data.frame(scenario = data_config, model_fit = model_fit, 
                          signal_size = params[3], nll_oos = oos_nll, 
                          nll_oos_oracle = oracle_oos_nll, 
                          nll_oos_oracle_diff = nll_oos_oracle_diff)

    return(list(results = results, cvobj = cv, is_data = is_data, oos_data = oos_data))
}) 

save(results_list, file = "results_list.Rdata")
load("results_list.Rdata")

results <- lapply(results_list, function(res) {
    res$results
})%>% bind_rows()

scenarios <- c("Linear Data", "Interaction Mean", 
               "Quadratic Mean", "Logistic Mean", 
               "Interaction Logit", "Quadratic Logit", 
               "Logistic Logit")

results$signal_size <- factor(results$signal_size, levels = 1:10)
results$scenario <- factor(results$scenario, 
                            levels = scenarios)
save(results, file = "results.Rdata")

##########################

load("results.Rdata")

p1 <- ggplot(filter(results, scenario == "Linear Data")) + 
    geom_line(aes(signal_size, nll_oos_oracle_diff, group = model_fit, col = model_fit, 
                  linetype = model_fit), linewidth = 1) + 
    scale_linetype_manual(values = c("Linear" = "dashed", "Nonlinear" = "solid")) + 
    theme(legend.title = element_text(size = 14), legend.text = element_text(size = 14),
        legend.key.size = unit(1.5, "lines")) + 
    labs(title = "Linear Data", y = "NLL (Above Oracle)", x = "", color = "Fitted Model", 
         linetype = "Fitted Model")

leg <- get_legend(p1)
leg <- as_ggplot(leg)

p1 <- p1 + theme(legend.position = "none")

p_list <- lapply(scenarios, function(cur_scenario) {
    if(cur_scenario == "Linear Data") {
        return(p1)
    } else {
        if(cur_scenario == "Interaction Logit") {
            y_lab <- "NLL (Above Oracle)"
        } else {
            y_lab <- ""
        }
        p <- ggplot(filter(results, scenario == cur_scenario)) + 
                geom_line(aes(signal_size, nll_oos_oracle_diff, group = model_fit, col = model_fit, 
                              linetype = model_fit), linewidth = 1) + 
                scale_linetype_manual(values = c("Linear" = "dashed", "Nonlinear" = "solid")) + 
                theme(legend.position = "none", 
                      plot.title = element_text(size = 20, margin = margin(t = 5, b = 10)),
                      axis.title = element_text(size = 18),  # Axis labels
                      axis.text = element_text(size = 14), 
                      axis.title.y.left = element_text(margin = margin(r = 10))) + 
                labs(title = cur_scenario, y = y_lab, x = "")
    }
})

p_list[[length(p_list) + 1]] <- leg
p_list[[9]] <- text_grob("Signal Size", size = 14, y = 0.9)
layout_mat <- matrix(c(1, 2, 3, 4, 8, 
                       NA, 5, 6, 7, 8, 
                       9, 9, 9, 9, 9), nrow = 3, byrow = TRUE)

pdf("plots/00_nll.pdf", 17.5, 7.9)
grid.arrange(grobs = p_list, layout_matrix = layout_mat, widths = c(1, 1, 1, 1, 0.4), 
            heights = c(1, 1, 0.1))
graphics.off()

ggplot(new_results) + 
    geom_line(aes(signal_size, nll_oos_oracle_diff, group = model_fit, col = model_fit)) + 
    geom_point(aes(signal_size, nll_oos_oracle_diff, group = model_fit, col = model_fit)) + 
    facet_wrap(~ scenario, nrow = 2, scales = "free") + 
    scale_x_continuous(breaks = seq(0, 1, 0.2))

# replace signal size index with actual signal size
load("signals.Rdata", verbose = TRUE)

# helper to replace signal size index with actual signal size
get_signals <- function(scenario) {
    switch(as.character(scenario), 
       `Linear Data` = 1, 
       `Interaction Mean` = 2, 
       `Quadratic Mean` = 3, 
       `Logistic Mean` = 4, 
       `Interaction Logit` = 1, 
       `Quadratic Logit` = 1, 
       `Logistic Logit` = 1)
}

new_results <- results %>%
    group_split(scenario) %>% 
    lapply(function(cur_df) {
        cur_scenario <- get_signals(cur_df$scenario[1])
        cur_df$signal_size <- signals[[cur_scenario]][cur_df$signal_size] # gets indices
        cur_df
    }) %>% 
    bind_rows()

p1 <- ggplot(filter(new_results, scenario == "Linear Data")) + 
    geom_line(aes(signal_size, nll_oos_oracle_diff, group = model_fit, col = model_fit, 
                  linetype = model_fit), linewidth = 1) + 
    scale_linetype_manual(values = c("Linear" = "dashed", "Nonlinear" = "solid")) + 
    scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) + 
    theme(legend.title = element_text(size = 14), legend.text = element_text(size = 14),
        legend.key.size = unit(1.5, "lines"), 
        legend.margin = margin(0, 0, 20, 60), 
        plot.title = element_text(size = 16, margin = margin(t = 5, b = 10)),
        axis.title = element_text(size = 16),  # Axis labels
        axis.text = element_text(size = 14), 
        axis.title.y.left = element_text(margin = margin(r = 10))) + 
    labs(title = "Linear Data", y = "NLL (Above Oracle)", x = "", color = "Fitted Model", 
         linetype = "Fitted Model")

leg <- get_legend(p1)
leg <- as_ggplot(leg)

p1 <- p1 + theme(legend.position = "none")

p_list <- lapply(scenarios, function(cur_scenario) {
    if(cur_scenario == "Linear Data") {
        return(p1)
    } else {
        if(cur_scenario == "Interaction Logit") {
            y_lab <- "NLL (Above Oracle)"
        } else {
            y_lab <- ""
        }
        p <- ggplot(filter(new_results, scenario == cur_scenario)) + 
                geom_line(aes(signal_size, nll_oos_oracle_diff, group = model_fit, col = model_fit, 
                              linetype = model_fit), linewidth = 1) + 
                scale_linetype_manual(values = c("Linear" = "dashed", "Nonlinear" = "solid")) + 
                scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) + 
                theme(legend.position = "none", 
                      plot.title = element_text(size = 16, margin = margin(t = 5, b = 10)),
                      axis.title = element_text(size = 16),  # Axis labels
                      axis.text = element_text(size = 14), ) + 
                labs(title = cur_scenario, y = y_lab, x = "")
    }
})

p_list[[length(p_list) + 1]] <- leg
p_list[[9]] <- text_grob("Signal Size", size = 14, y = 0.9)
layout_mat <- matrix(c(1, 2, 3, 4, 
                       8, 5, 6, 7, 
                       9, 9, 9, 9), nrow = 3, byrow = TRUE)

pdf("plots/00_nll_cont_signal.pdf", 16, 7.9)
grid.arrange(grobs = p_list, layout_matrix = layout_mat, widths = c(1, 1, 1, 1, 0.1), 
            heights = c(1, 1, 0.1))
graphics.off()

#################

# Showing the data

# panels:

# one of each mean (4)
# complete overlap, medium, reasonable separation (3) of linear model
# estimated mean for linear and nonlinear (2) (choose one good case)

rownames(results) <- as.character(1:nrow(results))

results %>% filter(scenario == "Linear Data", 
                   model_fit == "Linear", 
                   signal_size == 10)

results %>% filter(scenario == "Interaction Mean", 
                   model_fit == "Linear", 
                   signal_size == 10)

results %>% filter(scenario == "Quadratic Mean", 
                   model_fit == "Linear", 
                   signal_size == 10)

results %>% filter(scenario == "Logistic Mean", 
                   model_fit == "Linear", 
                   signal_size == 10)

# Compare mean functions
p1 <- plot_sim_means(results_list[[3]]$is_data, results_list[[3]]$cvobj, plot_model = FALSE) + 
    theme_gray() +
    theme(legend.position = "none") + 
    labs(x = "", title = "")

p2 <- plot_sim_means(results_list[[83]]$is_data, results_list[[83]]$cvobj, plot_model = FALSE) + 
    theme(legend.position = "none") + 
    labs(x = "", y = "", title = "")

p3 <- plot_sim_means(results_list[[103]]$is_data, results_list[[103]]$cvobj, plot_model = FALSE) + 
    theme(legend.position = "none") + 
    labs(x = "", y = "", title = "")

p4 <- plot_sim_means(results_list[[123]]$is_data, results_list[[123]]$cvobj, plot_model = FALSE) + 
    theme(legend.position = "none") + 
    labs(x = "", y = "", title = "")

x_lab <- text_grob("Time", size = 18, y = 0.9)
pdf("plots/00_means_comp.pdf", 22.7, 6.5)
grid.arrange(p1, p2, p3, p4, x_lab, layout_matrix = matrix(c(1, 2, 3, 4, 
                                                       5, 5, 5, 5), byrow = TRUE, nrow = 2), 
            heights = c(1, 0.025))
graphics.off()

# Compare Signal Sizes
results %>% filter(scenario == "Linear Data", 
                   model_fit == "Linear", 
                   signal_size == 1)

results %>% filter(scenario == "Linear Data", 
                   model_fit == "Linear", 
                   signal_size == 5)

p5 <- plot_sim_means(results_list[[1]]$is_data, results_list[[1]]$cvobj, plot_model = FALSE) + 
    theme_gray() +
    theme(legend.position = "none") + 
    labs(x = "", y = "", title = "")

p6 <- plot_sim_means(results_list[[11]]$is_data, results_list[[11]]$cvobj, plot_model = FALSE) + 
    theme(legend.position = "none") + 
    labs(x = "", y = "", title = "")

grid.arrange(p6, p5, nrow = 1)

# Compare Estimated Models
results %>% filter(scenario == "Interaction Mean", 
                   model_fit == "Linear", 
                   signal_size == 7)

results %>% filter(scenario == "Interaction Mean", 
                   model_fit == "Nonlinear", 
                   signal_size == 7)

p7 <- plot_sim_means(results_list[[95]]$is_data, results_list[[95]]$cvobj, plot_model = TRUE) + 
    theme_gray() +
    # theme(legend.position = "none") + 
    guides(fill = "none") + 
    labs(x = "", y = "Data (Interaction Mean)", title = "Linear Model Fit") + 
    theme(legend.title = element_text(size = 18),   # Title size
          legend.text = element_text(size = 16),    # Item text size
          legend.key.size = unit(2, "lines"))

leg <- get_legend(p7) %>% 
    as_ggplot()

p7 <- p7 + 
    theme(legend.position = "none", 
          plot.title = element_text(size = 20, margin = margin(t = 5, b = 10)),
          axis.title = element_text(size = 18),  # Axis labels
          axis.text = element_text(size = 14), 
          axis.title.y.left = element_text(margin = margin(r = 10))) # , axis.text.y = element_text(margin = margin(r = 15)))

p8 <- plot_sim_means(results_list[[96]]$is_data, results_list[[96]]$cvobj, plot_model = TRUE) + 
    theme_gray() + 
    theme(legend.position = "none", 
          plot.title = element_text(size = 20, margin = margin(t = 5, b = 10)),
          axis.title = element_text(size = 18),  # Axis labels
          axis.text = element_text(size = 14)) + 
    labs(x = "", y = "", title = "Nonlinear Model Fit")

p9 <- plot_sim_probs(results_list[[95]]$is_data, results_list[[95]]$cvobj) + 
        theme(legend.position = "none") + 
        labs(x = "", y = "Cluster Probability", title = "")

p10 <- plot_sim_probs(results_list[[96]]$is_data, results_list[[96]]$cvobj) + 
        theme(legend.position = "none") + 
        labs(x = "", y = "", title = "")

pdf("plots/00_inter_mean_model_comp.pdf", 10.7, 9.4)
grid.arrange(p7, p8, x_lab, p9, p10, x_lab, leg, layout_matrix = matrix(c(1, 2, 7,
                                                                         3, 3, 7,
                                                                         4, 5, 7,
                                                                         6, 6, 7), byrow = TRUE, 
                                                           nrow = 4), 
            heights = c(1, 0.02, 1, 0.1), widths = c(1, 1, 0.3))
graphics.off()

pdf("plots/00_mega_rough.pdf", 8.1, 12.3)
grid.arrange(p1, p2, p3, p4, 
             p5, p6, 
             p7, p8, 
             p9, p10, 
             x_lab, leg, layout_matrix = matrix(c(1, 2, 3, 4, 12,
                                                  5, 5, 6, 6, 12,
                                                  7, 7, 8, 8, 12,
                                                  9, 9, 10, 10, 12,
                                                  11, 11, 11, 11, 12), byrow = TRUE, nrow = 5), 
             widths = c(1, 1, 1, 1, 0.5), heights = c(3, 4, 4, 4, 0.2))
graphics.off()

x_lab <- text_grob("Time", size = 18, x = 0.52, y = 0.9)
pdf("plots/00_sim_data.pdf", 14, 6.5)
grid.arrange(p7, p8,
             x_lab, leg, layout_matrix = matrix(c(1, 2, 4,
                                                  3, 3, 4), byrow = TRUE, nrow = 2), 
             widths = c(1, 1, 0.45), heights = c(1, 0.1))
graphics.off()

grid.arrange(p7, p8,
             x_lab, leg, layout_matrix = matrix(c(1, 4,
                                                  2, 4, 
                                                  3, 4), byrow = TRUE, nrow = 3), 
             widths = c(1, 0.25), heights = c(1, 1, 0.1))

# TODO: continue from here. Stitch all the plots together into one big plot
# TODO: then make additional distinguishing marks by mapping cluster to linetype
# TODO: finally, upload to Overleaf, then alter text size accordingly

##################

cvscore_heatmap <- function(sim) {
    res <- readRDS(file.path("results", paste0(sim, "_summary.RDS")))
    print(res$min.inds)
    corrplot(res$cvscore.mat, is.corr = FALSE, method = "color", 
             col = COL2('RdBu', 200))

}

min_inds <- function(sim) {
    res <- readRDS(file.path("results", paste0(sim, "_summary.RDS")))
    res$min.inds
}

cvscore_heatmap("1-1-2-l")
cvscore_heatmap("1-1-3-l")
cvscore_heatmap("1-1-4-l")

cvscore_heatmap("1-1-8-l")
cvscore_heatmap("1-1-9-l")
cvscore_heatmap("1-1-10-l")

cvscore_heatmap("2-1-1-l")
cvscore_heatmap("2-1-2-l")
cvscore_heatmap("2-1-3-l")
cvscore_heatmap("2-1-7-l")
cvscore_heatmap("2-1-8-l")

cvscore_heatmap("3-1-1-l")
cvscore_heatmap("3-1-2-l")
cvscore_heatmap("3-1-3-l")

cvscore_heatmap("3-1-1-n")
cvscore_heatmap("3-1-2-n")
cvscore_heatmap("3-1-3-n")

cvscore_heatmap("4-1-1-l")
cvscore_heatmap("4-1-2-l")
cvscore_heatmap("4-1-3-l")
cvscore_heatmap("4-1-9-l")

cvscore_heatmap("4-1-1-n")
cvscore_heatmap("4-1-2-n")
cvscore_heatmap("4-1-3-n")

cvscore_heatmap("1-2-1-l")
cvscore_heatmap("1-2-2-l")
cvscore_heatmap("1-2-3-l")
cvscore_heatmap("1-2-7-l")
cvscore_heatmap("1-2-8-l")

cvscore_heatmap("1-3-2-l")
cvscore_heatmap("1-3-2-n")

cvscore_heatmap("1-4-1-n")
cvscore_heatmap("1-4-2-n")

# need to decrease prob regularization
# do this for a few of the worst cases to define a large-scale solution

# first: 1-1-2-l

res_dir <- dir_create("results_retry")

library(parallel)
library(parallelly)
make_seedtab <- function(imean, iprob, iint, modelFit) {
    RNGkind("L'Ecuyer-CMRG")
    nalpha <- 10
    nbeta <- 10
    nfold <- 5
    nrep <- 10
    nrows <- nalpha * nbeta * (nfold + 1) * nrep

    seed_destin <- dir_create("seedtabs_retry")

    seedfile <- file.path(seed_destin,
                          paste0(imean, "-", iprob, "-", iint, "-", modelFit, "_seedtab.csv"))
    
    if(file.exists(seedfile)) {
        print("seedtab already exists.")
    } else {
        ## Generate the random number /states/ (7 integers each)
        set.seed(NULL)
        s <- list(.Random.seed)
        for (ii in 2:nrows){
            s[[ii]] <- nextRNGStream(s[[ii-1]])
        }
        s = s %>% do.call(rbind, .) %>% as_tibble()
        colnames(s) = paste0("seed", 1:7)
        # ifold == 0 corresponds to the refit step
        tab = expand.grid(ialpha = 1:nalpha, 
                            ibeta = 1:nbeta, 
                            ifold = 0:nfold, 
                            irep = 1:nrep) %>% as_tibble()
        seedtab = tab %>% bind_cols(s)
        ## Write the table to file
        write.csv(seedtab, file = seedfile, row.names = FALSE)
        print(paste0("Wrote table containing seeds to", seedfile))
    }
    RNGkind("default")
}

RhpcBLASctl::blas_set_num_threads(1)
sim_dir <- dir_create(file.path(res_dir, "1-1-2-l"))

imean <- 1
iprob <- 1
iint <- 2
modelFit <- "l"

make_seedtab(imean, iprob, iint, modelFit)

load(file.path("~", 
               "00_Cyto", 
               "data", 
               "simdata", 
               "simdata_04_01",
               paste0("simdata-", imean, "-", iprob, "-", iint, "-", "1", ".Rdata")))

ylist <- res$ybin_list
countslist <- res$countslist
maxdev <- diff(range(res$mean_spec)) / 2

numclust <- 2
nfold <- 5 
nrep <- 10
blocksize <- 20
seedtab <- read.csv(file.path("seedtabs_retry",
                              paste0(imean, "-", iprob, "-", iint, "-", modelFit, "_seedtab.csv")))

destin <- file.path("results_retry", 
                    paste0(imean, "-", iprob, "-", iint, "-", modelFit))

if(!dir.exists(destin)) {
    dir.create(destin, recursive = TRUE)
}

load(file.path("~", 
               "00_Cyto", 
               "data",
               "X_data", 
                switch(modelFit, 
                       l = "X_pc.Rdata", 
                       n = "X_nl_70.Rdata")))
# either way, loads in an object named X

prob_lambdas <- c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4)
mean_lambdas <- c(0.00196, 0.00528, 0.0142, 0.0383, 0.103)

save(prob_lambdas, file = file.path(destin, "prob_lambdas.Rdata"))
save(mean_lambdas, file = file.path(destin, "mean_lambdas.Rdata"))

folds <- make_cv_folds(ylist, nfold, blocksize)

cv.flowmix(ylist, 
           countslist, 
           X, 
           destin, 
           mean_lambdas, 
           prob_lambdas,
           NULL, 
           maxdev, 
           numclust, 
           nfold, 
           nrep, 
           TRUE, 
           FALSE, 
           TRUE, 
           availableCores() - 1, 
           blocksize, 
           folds, 
           seedtab)

prob_lambdas <- c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4)
mean_lambda <- 0.0383

cvreslist <- lapply(1:length(prob_lambdas), function(i) {
    cur_destin <- dir_create(file.path(destin, i))
    prob_lambda <- prob_lambdas[i]

    cvres <- cv.flowmix(ylist, 
               countslist, 
               X, 
               cur_destin, 
               mean_lambda, 
               prob_lambda,
               NULL, 
               maxdev, 
               numclust, 
               nfold, 
               nrep, 
               TRUE, 
               FALSE, 
               TRUE, 
               availableCores() - 1, 
               blocksize, 
               folds, 
               seedtab)
    
    refitres <- cv.flowmix(ylist, 
               countslist, 
               X, 
               cur_destin, 
               mean_lambda, 
               prob_lambda,
               NULL, 
               maxdev, 
               numclust, 
               nfold, 
               nrep, 
               TRUE, 
               TRUE, 
               FALSE, 
               availableCores() - 1, 
               blocksize, 
               folds, 
               seedtab)

    cvres <- cv_summary(destin = cur_destin, 
           save = TRUE, 
           filename = paste0(imean, "-", iprob, "-", iint, "-", modelFit, "_summary.RDS"))
    
    print(cvres$cvscore.mat)
})

res1 <- readRDS("results/1-1-2-l_summary.RDS")

# TODO: keep going here
cv.flowmix(ylist, 
           countslist, 
           X, 
           destin, 
           mean_lambdas, 
           prob_lambdas,
           NULL, 
           maxdev, 
           numclust, 
           nfold, 
           nrep, 
           TRUE, 
           TRUE, 
           FALSE, 
           availableCores() - 1, 
           blocksize, 
           folds, 
           seedtab)

# TODO: may need to try a finer grid around (higher and lower than) 1e-4
#       - for one scenario, try different things until NLL is lowered
#       - TODO:  try a very fine grid between 
#                  prob_lambdas[2] and prob_lambdas[1]
#       - linear grid may be best

plot(1:3, res1$cvscore.mat[1:3,7])

cv_summary(destin = destin, 
           save = TRUE, 
           filename = paste0(imean, "-", iprob, "-", iint, "-", modelFit, "_summary.RDS"))
