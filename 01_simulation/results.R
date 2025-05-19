library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(flowmix)
library(tibble)
library(ggpubr)
library(RColorBrewer)

# Collect results
simdata_dir <- file.path("01_simulation", "simdata")
cvres_files <- list.files(file.path("01_simulation", "results"), recursive = TRUE)

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

# Use RColorBrewer palette in plotting
brew_colors <- brewer.pal(n = 8, name = "Set1")

# Assign specific colors to each category
named_colors <- c("1" = brew_colors[1],
                  "2" = brew_colors[2],
                  "1, Truth" = brew_colors[4],
                  "2, Truth" = brew_colors[5])

# Labels for the legend
labels <- c("1" = "1, Estimated",
            "2" = "2, Estimated",
            "1, Truth" = "1, Truth",
            "2, Truth" = "2, Truth")

# Plot training simdata and overlay estimated parameters
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

results_list_file <- file.path("01_simulation", "results", "00_results_list.Rdata")

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

results_file <- file.path("01_simulation", "results", "00_results.Rdata")

if(!file.exists(results_file)) {
    # Collect model evaluation results into one data frame 
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
    save(results, file = results_file)
} else {
    # Load the results data frame if this step has already been done
    load(results_file)
}

##########################

# Figures

# Figure 2

load(file.path("01_simulation", "aux_simdata", "signals.Rdata"))

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

# Create the remaining panels
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

# Attach legend and x-axis label
p_list[[length(p_list) + 1]] <- leg
p_list[[9]] <- text_grob("Signal Size", size = 14, y = 0.9)
layout_mat <- matrix(c(1, 2, 3, 4, 
                       8, 5, 6, 7, 
                       9, 9, 9, 9), nrow = 3, byrow = TRUE)

# Figure 2
pdf(file.path("plots", "nll_cont_signal.pdf"), 16, 7.9)
grid.arrange(grobs = p_list, layout_matrix = layout_mat, widths = c(1, 1, 1, 1, 0.1), 
            heights = c(1, 1, 0.1))
graphics.off()

#################

# Figure 1

p1 <- plot_sim_means(results_list[[95]]$is_data, results_list[[95]]$cvobj, plot_model = TRUE) + 
    theme_gray() +
    # theme(legend.position = "none") + 
    guides(fill = "none") + 
    labs(x = "", y = "Data (Interaction Mean)", title = "Linear Model Fit") + 
    theme(legend.title = element_text(size = 18),   # Title size
          legend.text = element_text(size = 16),    # Item text size
          legend.key.size = unit(2, "lines"))

leg <- get_legend(p1) %>% 
    as_ggplot()

p1 <- p1 + 
    theme(legend.position = "none", 
          plot.title = element_text(size = 20, margin = margin(t = 5, b = 10)),
          axis.title = element_text(size = 18),  # Axis labels
          axis.text = element_text(size = 14), 
          axis.title.y.left = element_text(margin = margin(r = 10))) # , axis.text.y = element_text(margin = margin(r = 15)))

p2 <- plot_sim_means(results_list[[96]]$is_data, results_list[[96]]$cvobj, plot_model = TRUE) + 
    theme_gray() + 
    theme(legend.position = "none", 
          plot.title = element_text(size = 20, margin = margin(t = 5, b = 10)),
          axis.title = element_text(size = 18),  # Axis labels
          axis.text = element_text(size = 14)) + 
    labs(x = "", y = "", title = "Nonlinear Model Fit")

x_lab <- text_grob("Time", size = 18, x = 0.52, y = 0.9)
pdf("plots/00_sim_data.pdf", 14, 6.5)
grid.arrange(p1, p2,
             x_lab, leg, layout_matrix = matrix(c(1, 2, 4,
                                                  3, 3, 4), byrow = TRUE, nrow = 2), 
             widths = c(1, 1, 0.45), heights = c(1, 0.1))
graphics.off()
