isim <- as.integer(commandArgs(trailingOnly = TRUE)[1])

library(dotty)
library(flowmix)
library(RhpcBLASctl)

blas_set_num_threads(1)
omp_set_num_threads(1)

res_dir <- file.path("1_simulation", "NLPLs")
if(!dir.exists(res_dir)) dir.create(res_dir)


hps_and_seeds <- readRDS(
  file.path("1_simulation", "hyperparameters_and_seeds.RDS")
)

.[prob_lambda, mean_lambda, seed] <- hps_and_seeds[[isim]]

sim_name <- names(hps_and_seeds)[[isim]] 
sim_name_split <- sim_name |> 
  strsplit("-")
sim_name_split <- sim_name_split[[1]]
.[imean, iprob, iint, nh, NNseed] <- sim_name_split

simdata <- readRDS(
  file.path(
    "1_simulation", "simdata", 
    paste(
      "simdata_imean", imean, "iprob", iprob, "iint", iint, 
      "dataSeed_1.RDS", sep = "_"
    )
  )
)

ybin_list <- simdata$ybin_list 
countslist <- simdata$countslist 
numclust <- 2
maxdev <- diff(range(simdata$mean_spec)) / 2

X <- readRDS(
  file.path(
    "data", 
    "X_variations",
    paste0("X_pc_9_nh_", nh, "_seed_", NNseed, "_ofold_NA_ifold_NA.RDS")
  )
)

res <- flowmix_once(
  ylist = ybin_list, 
  X = X, 
  countslist = countslist, 
  numclust = numclust, 
  prob_lambda = prob_lambda, 
  mean_lambda = mean_lambda, 
  verbose = TRUE, 
  maxdev = maxdev, 
  seed = seed
)

# Fix potential label switching b/w clusters 1 and 2
clust_1_pos <- res$mn[,1,] |>
  apply(2, min) |>
  which.min()

if(clust_1_pos == 2) {
  res$alpha <- res$alpha[2:1,1:10]
  res$beta <- res$beta[2:1]
  res$mn <- res$mn[,,2:1, drop = FALSE]
  res$prob <- res$prob[,2:1, drop = FALSE]
  res$sigma <- res$sigma[2:1, 1, 1, drop = FALSE]
}

if(paste(sim_name_split[1:3], collapse = "-") == "2-1-5") {
  res_dir <- file.path("1_simulation", "results", "2-1-5")
  if(!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)
  saveRDS(res, file.path(res_dir, paste0(sim_name, "_fit_model.RDS")))
}

oos_data <- readRDS(
  file.path(
    "1_simulation", "simdata", 
    paste(
      "simdata_imean", imean, "iprob", iprob, "iint", iint, 
      "dataSeed_0.RDS", sep = "_"
    )
  )
)

# Calculate out-of-sample negative log-pseudolikelihoods (NLPLs) in a data frame
oos_nll <- objective_newdat(res, oos_data$ybin_list, oos_data$countslist)

# Create oracle model (the model which generates the data)
TT <- length(ybin_list)
oracle <- list(
  mn = array(c(simdata$mean_spec, simdata$pico_mu), dim = c(TT, 1, numclust)), 
  prob = array(c(simdata$prob_spec, 1 - simdata$prob_spec), dim = c(TT, numclust)),
  sigma = array(simdata$clust_sig^2, dim = c(numclust, 1, 1))
)

# Get oracle NLPL
oracle_oos_nll <- objective_newdat(oracle, oos_data$ybin_list, oos_data$countslist)               

# Get model NLPL, less oracle NLPL
nll_oos_oracle_diff <- oos_nll - oracle_oos_nll

model <- if(nh == "NA") "Linear" else "Nonlinear"
scenario <- switch(paste0(imean, iprob), 
  "11" = "Linear", 
  "12" = "Interaction in Logit", 
  "13" = "Quadratic Logit", 
  "14" = "Logistic Logit", 
  "21" = "Interaction in Mean", 
  "31" = "Quadratic Mean", 
  "41" = "Logistic Mean"
)

nll_df <- data.frame(
  `NLPL (Above Oracle)` = nll_oos_oracle_diff, 
  Model = model, 
  `Signal Size` = as.integer(iint), 
  Scenario =  scenario, 
  check.names = FALSE
)

saveRDS(nll_df, file.path(res_dir, paste0(sim_name, ".RDS")))
