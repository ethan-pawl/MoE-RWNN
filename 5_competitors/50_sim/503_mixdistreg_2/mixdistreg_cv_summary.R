cvscore_dir <- file.path("5_competitors", "50_sim", "503_mixdistreg", "cvscores")

l1_grid <- c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100)
cv_gridsize <- length(l1_grid)
nfold <- 5

cvscore_arr <- array(NA, c(cv_gridsize, cv_gridsize, nfold))
dimnames(cvscore_arr) <- list(
  ialpha = l1_grid, 
  ibeta = l1_grid, 
  ifold = 1:nfold
)

for(i in 1:cv_gridsize) {
  for(j in 1:cv_gridsize) {
    for(k in 1:nfold) {
      fname <- file.path(cvscore_dir, paste(i, j, k, "cvscore.Rdata", sep = "-"))
      if(file.exists(fname)) {
        load(fname)
      # Loads ialpha, ibeta, cvscore
        cvscore_arr[ialpha,ibeta,ifold] <- cvscore
      }
    }
  }
}

cvscore_mat <- apply(cvscore_arr, c(1, 2), function(x) mean(x, na.rm = TRUE))

png("corrplot.png")
corrplot::corrplot(log(cvscore_mat + abs(min(cvscore_mat)) + 1), method = "color", is.corr = FALSE, xlab = l1_grid, ylab = l1_grid)
graphics.off()

which.min(cvscore_mat) |> arrayInd(c(10, 10), dimnames(cvscore_mat), useNames = TRUE)
#       ialpha ibeta
# 1e-05      3     4
rownames(cvscore_mat)[3] # 1e-5
colnames(cvscore_mat)[4] # 1e-4
