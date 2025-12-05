library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(flowmix)
library(tibble)
library(ggpubr)
library(RColorBrewer)
library(ggtext)
library(grid)

# TODO: alter this to work correctly with the new setup

# Collect results
simdata_dir <- file.path("1_simulation", "simdata")
cvres_files <- list.files(file.path("1_simulation", "results"), 
                          recursive = TRUE, 
                          pattern = "*.RDS", 
                          full.names = TRUE)

###################################

# Helper Functions

# start with simulation results file name, 
# load training simdata and return
get_simdata <- function(fname) {
    simdata_fname <- fname %>% 
        strsplit("-") %>%
        unlist() %>% 
        .[1:3] %>% 
        paste(collapse = "-") %>% 
        paste0("simdata-", ., "-1.Rdata") %>%
        file.path(simdata_dir, .)
    
    load(simdata_fname)
    simdata
}

# start with simulation results file name, 
# load test simdata and return
get_oos <- function(fname) {
    simdata_fname <- fname %>% 
        strsplit("-") %>%
        unlist() %>% 
        .[1:3] %>% 
        paste(collapse = "-") %>% 
        paste0("simdata-", ., "-0.Rdata") %>%
        file.path(simdata_dir, .)
    
    load(simdata_fname)
    simdata
}

# Convert data configuration indices to 
# a plot title in words 
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

results_list_file <- file.path("1_simulation", "results", "00_results_list.Rdata")

if(!file.exists(results_list_file)) {
    # Gather results from every simulation into one big list 
    # and one big data frame
    results_list <- lapply(cvres_files, function(cv_file) {
        # Load CV results
        cv <- readRDS(file.path("results", cv_file))
        is_data <- get_simdata(cv_file)
        oos_data <- get_oos(cv_file)

        # Fix label switching between cluster 1 and cluster 2
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

        # Calculate out-of-sample negative log-likelihoods (NLLs) in a data frame
        oos_nll <- objective_newdat(cv$bestres, oos_data$ybin_list, oos_data$countslist)

        # Create oracle model (the model which generates the data)
        oracle <- list(mn = array(c(is_data$mean_spec, is_data$pico_mu), dim = c(296, 1, 2)), 
                       prob = array(c(is_data$prob_spec, 1 - is_data$prob_spec), dim = c(296, 2)), 
                       sigma = array(is_data$clust_sig^2, dim = c(2, 1, 1)))

        # Get oracle NLL
        oracle_oos_nll <- objective_newdat(oracle, oos_data$ybin_list, oos_data$countslist)               

        # Get model NLLs, less oracle NLL
        nll_oos_oracle_diff <- oos_nll - oracle_oos_nll

        # Extract data configuration indices
        params <- cv$destin %>% 
                strsplit("/") %>% 
                unlist() %>%
                .[2] %>% 
                strsplit("-") %>% 
                unlist()

        # Convert indices to descriptive data configuration names
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

        # Collect model evaluation results in one row of a data frame
        results <- data.frame(scenario = data_config, model_fit = model_fit, 
                              signal_size = params[3], nll_oos = oos_nll, 
                              nll_oos_oracle = oracle_oos_nll, 
                              nll_oos_oracle_diff = nll_oos_oracle_diff)

        # return eval results, estimated model, training data, and test data
        return(list(results = results, cvobj = cv, is_data = is_data, oos_data = oos_data))
    }) 

    save(results_list, file = results_list_file)
} else {
    # Load the results list if this step has already been done
    load(results_list_file)
}

results_file <- file.path("1_simulation", "results", "00_results.Rdata")

scenarios <- c("Linear Data", "Interaction Mean", 
               "Quadratic Mean", "Logistic Mean", 
               "Interaction Logit", "Quadratic Logit", 
               "Logistic Logit")

if(!file.exists(results_file)) {
    # Collect model evaluation results into one data frame 
    results <- lapply(results_list, function(res) {
        res$results
    })%>% bind_rows()

    results$signal_size <- factor(results$signal_size, levels = 1:10)
    results$scenario <- factor(results$scenario, 
                                levels = scenarios)
    save(results, file = results_file)
} else {
    # Load the results data frame if this step has already been done
    load(results_file)
}

##########################

# Figures

# Figure 2

load(file.path("1_simulation", "aux_simdata", "signals.Rdata"))

# helper to map each scenario to its 
# corresponding vector of signal sizes
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

# Replace signal size index with actual signal size 
# (distance between clusters)
new_results <- results %>%
    group_split(scenario) %>% 
    lapply(function(cur_df) {
        cur_scenario <- get_signals(cur_df$scenario[1]) # get vector
        cur_df$signal_size <- signals[[cur_scenario]][cur_df$signal_size] # index vector to get exact signal size
        cur_df
    }) %>% 
    bind_rows()

# Plot each panel individually

# Plot first panel to extract legend
p1 <- ggplot(filter(new_results, scenario == "Linear Data")) + 
    geom_line(aes(signal_size, nll_oos_oracle_diff, group = model_fit, col = model_fit, 
                  linetype = model_fit), linewidth = 1) + 
    scale_linetype_manual(values = c("Linear" = "dotdash", "Nonlinear" = "solid")) + 
    scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) + 
    theme(legend.title = element_text(size = 18, hjust = 0.4), legend.text = element_text(size = 16),
        legend.key.size = unit(2.5, "lines"), 
        legend.margin = margin(0, 0, 20, 60), 
        plot.title = element_textbox_simple(
          size = 20,
          margin = margin(t = 5, b = 10),
          face = "bold",
          fill = "grey85",
          box.color = NA,
          padding = margin(3, 6, 0.5, 6), 
          halign = 0.5
        ),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), 
        axis.title.y.left = element_text(margin = margin(r = 10))) + 
    labs(title = "Linear Means\nand Probabilities", y = "NLPL (Above Oracle)", x = "", 
         color = "Model", 
         linetype = "Model") + 
    guides(color = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2))

leg <- get_legend(p1)
leg <- as_ggplot(leg)

p1 <- p1 + 
    theme(legend.position = "none")

# Create the remaining panels
p_list <- lapply(scenarios, function(cur_scenario) {
    if(cur_scenario == "Linear Data") {
        return(p1)
    } else {
        p <- ggplot(filter(new_results, scenario == cur_scenario)) + 
                geom_line(aes(signal_size, nll_oos_oracle_diff, group = model_fit, col = model_fit, 
                              linetype = model_fit), linewidth = 1) + 
                scale_linetype_manual(values = c("Linear" = "dotdash", "Nonlinear" = "solid")) + 
                scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) + 
                theme(legend.position = "none", 
                      plot.title = element_text(size = 16, margin = margin(t = 5, b = 10)),
                      axis.title = element_text(size = 16),
                      axis.text = element_text(size = 14)) + 
                labs(title = strsplit(cur_scenario, " ")[[1]][1], y = "", x = "")
    }
})

# Attach legend and x-axis label
p_list[[length(p_list) + 1]] <- leg
p_list[[9]] <- text_grob(expression("Signal Size (" * Delta * ")"), size = 14, x = 0.53, y = 0.9)
layout_mat <- matrix(c(1, 2, 3, 4, 
                       8, 5, 6, 7, 
                       9, 9, 9, 9), nrow = 3, byrow = TRUE)

p_list[[10]] <- grobTree(
    rectGrob(gp = gpar(fill = "grey85", col = NA), height = unit(0.75, "cm"), 
             y = 0.275, width = unit(0.92, "npc"), x = unit(0.534, "npc")),
    text_grob("Nonlinear Means", size = 20, face = "bold", just = "top", 
              x = 0.5075)
)

p_list[[11]] <- grobTree(
    rectGrob(gp = gpar(fill = "grey85", col = NA), height = unit(0.75, "cm"), 
             y = 0.275, width = unit(0.91, "npc"), x = unit(0.5375, "npc")),
    text_grob("Nonlinear Probabilities", size = 20, face = "bold", just = "top", 
        x = 0.5325, y = 0.8)
)

layout_mat <- matrix(c(10, 10, 10, 10,
                       1, 2, 3, 4, 
                       11, 11, 11, 11,
                       8, 5, 6, 7, 
                       9, 9, 9, 9), nrow = 5, byrow = TRUE)

p_list[[12]] <- nullGrob()
layout_mat <- matrix(c(12, 10, 10, 10,
                       12, 2, 3, 4, 
                       1, 2, 3, 4, 
                       1, 11, 11, 11,
                       1, 12, 12, 12,
                       1, 5, 6, 7,
                       8, 5, 6, 7, 
                       9, 9, 9, 9), ncol = 4, byrow = TRUE)

pdf(file.path("plots", "Figure03.pdf"), 16, 7.9)
grid.arrange(grobs = p_list, layout_matrix = layout_mat, 
             heights = c(0.25, 0.89, 1.11, 0.1, 0.05, 1.11, 0.89, 0.15), 
             widths = c(1.05, 1, 1, 1))
graphics.off()

#################

# Figure 2

brew_colors <- brewer.pal(n = 8, name = "Set1")

named_colors <- c("1" = brew_colors[5],
                  "2" = brew_colors[2],
                  "Truth" = "black", 
                  "True 95% CI" = "black")

named_linetypes <- c("1" = "solid", 
                     "2" = "solid", 
                    "Truth" = "twodash", 
                    "True 95% CI" = "dotted")

# Labels for the legend
labels <- c("1" = "Cluster 1 Estimated Mean",
            "2" = "Cluster 2 Estimated Mean",
            "Truth" = "True Mean", 
            "True 95% CI" = "True 95% Probability Region")

plot_sim_means <- function(simdat, cv, plot_band = TRUE, plot_model = TRUE) {
    if(plot_model) {
        flowtrend::plot_1d(simdat$ybin_list, simdat$countslist, cv$bestres, 
                           bin = TRUE, plot_band = plot_band) + 
            geom_line(data = data.frame(x = rep(1:296, 2), 
                                        y = c(simdat$mean_spec, simdat$pico_mu), 
                                        prob = c(simdat$prob_spec, 1 - simdat$prob_spec),
                                        clust = "Truth",
                                        clust_group = rep(c("1, Truth", "2, Truth"), 
                                                    each = 296)), 
                      mapping = aes(x = x, y = y, group = clust_group, 
                                    color = clust, linetype = clust), 
                      lineend = "round", linewidth = 1) + 
            scale_color_manual(name = "", 
                               values = named_colors, 
                               labels = labels) + 
            scale_linetype_manual(name = "", 
                                  values = named_linetypes, 
                                  labels = labels) + 
            scale_linewidth_continuous(name = "Cluster Probability", 
                                       range = c(1, 4)) + 
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

true_1_upper <- results_list[[91]]$is_data$mean_spec[,1] + 1.96 * 0.2
true_1_lower <- results_list[[91]]$is_data$mean_spec[,1] - 1.96 * 0.2

true_2_upper <- results_list[[91]]$is_data$pico_mu + 1.96 * 0.2
true_2_lower <- results_list[[91]]$is_data$pico_mu - 1.96 * 0.2

ci_df <- data.frame(x = rep(1:296, 4), 
                    y = c(true_1_upper, true_1_lower, true_2_upper, true_2_lower), 
                    Cluster = rep(c("1, Truth", "2_Truth"), each = 296 * 2), 
                    type = rep(rep(c("upper", "lower"), each = 296), 2), 
                    clust = "True 95% CI")

p11 <- plot_sim_means(results_list[[91]]$is_data, results_list[[91]]$cvobj, plot_model = TRUE) + 
    geom_line(data = ci_df, 
             mapping = aes(x = x, y = y,
                           group = interaction(Cluster, type), 
                           linetype = clust, 
                           color = clust), 
             lineend = "round", 
             linewidth = 1) + 
    theme_gray() + 
    theme(plot.title = element_text(size = 20, margin = margin(t = 5, b = 10)), 
          axis.title = element_text(size = 18), 
          axis.text = element_text(size = 14),
          legend.position = "right",
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18, margin = margin(r = 50)), 
          legend.key.size = unit(2, "lines"), 
          legend.margin = margin(0, 0, 0, 0)) + 
    guides(color = guide_legend(override.aes = list(linetype = c("solid",  
                                                                 "solid", 
                                                                 "dotted", 
                                                                 "twodash")), 
                                nrow = 2, ncol = 2),
           linetype = "none",
           linewidth = "none",
           fill = "none") + 
    scale_y_continuous(limits = c(0.25, 2)) + 
    labs(x = "", y = "Data", title = "Estimated Linear Model")

leg <- get_legend(p11) %>% 
    as_ggplot()

p11 <- p11 + 
    theme(legend.position = "none")

p12 <- plot_sim_means(results_list[[92]]$is_data, results_list[[92]]$cvobj, plot_model = TRUE) + 
    geom_line(data = ci_df, 
             mapping = aes(x = x, y = y,
                           group = interaction(Cluster, type), 
                           linetype = clust, 
                           color = clust), 
             lineend = "round", 
             linewidth = 1)  + 
    scale_y_continuous(limits = c(0.25, 2)) + 
    theme_gray() + 
    theme(legend.position = "none", 
          plot.title = element_text(size = 20, margin = margin(t = 5, b = 10)),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14)) + 
    labs(x = "", y = "", title = "Estimated Nonlinear Model")

x_lab <- text_grob("Time (t)", size = 18, x = 0.53, y = 0.9)

pdf(file.path("plots", "Figure02.pdf"), 14.6, 8.5)
grid.arrange(p11, p12,
             x_lab, leg, nullGrob(), 
             layout_matrix = matrix(c(1, 1, 2, 2,
                                      3, 3, 3, 3,
                                      5, 4, 4, 5), byrow = TRUE, nrow = 3), 
             widths = c(0.27, 0.73, 0.95, 0.05), heights = c(1, 0.01, 0.2))
graphics.off()
