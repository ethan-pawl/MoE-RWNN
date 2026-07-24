cvscore_dir <- file.path("5_competitors", "50_sim", "503_mixdistreg", "cvscores")

l1_grid <- c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100)

cvscore_mat <- matrix(NA, 10, 10)
dimnames(cvscore_mat) <- list(
  ialpha = l1_grid, 
  ibeta = l1_grid
)

for(i in 1:10) {
  for(j in 1:10) {
    load(file.path(cvscore_dir, paste0(i, "-", j, "-cvscore.Rdata")))

    # Loads ialpha, ibeta, cvscore
    cvscore_mat[ialpha,ibeta] <- cvscore
  }
}

corrplot::corrplot(log(cvscore_mat + abs(min(cvscore_mat)) + 1), method = "color", is.corr = FALSE, xlab = l1_grid, ylab = l1_grid)
which.min(cvscore_mat) |> arrayInd(c(10, 10), dimnames(cvscore_mat), useNames = TRUE)
#       ialpha ibeta
# 1e-05      3     5
rownames(cvscore_mat)[3] # 1e-5
colnames(cvscore_mat)[5] # 1e-3
