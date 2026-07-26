cvscore_dir <- file.path("5_competitors", "50_sim", "501_gam", "cvscores")

cv_gridsize <- 10
nfold <- 5

cvscore_arr <- array(NA, c(cv_gridsize, nfold))

for(i in 1:cv_gridsize) {
  for(k in 1:nfold) {
    fname <- file.path(cvscore_dir, paste(i, k, "cvscore.Rdata", sep = "-"))
    if(file.exists(fname)) {
      load(fname)
      # Loads igamma, ifold, cvscore
      cvscore_arr[igamma,ifold] <- cvscore
    }
  }
}

saveRDS(cvscore_arr, file = file.path("5_competitors", "50_sim", "501_gam", "cv-summary.RDS"))