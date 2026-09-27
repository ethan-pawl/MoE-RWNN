num_PCs <- c("9", "18", "27", "37")
NNseeds <- c(4, 1, 2, 1)

cvscores <- vapply(1:4, function(i) {
  fname <- paste("nl", num_PCs[i], "70", NNseeds[i], sep = "_")
  destin <- file.path("2_application", "3_PCs", "results", fname)
  
  cvres <- flowmix::cv_aggregate(destin)
  return(cvres$cvscore.mat[1])

}, numeric(1))

mean(cvscores)
sd(cvscores) / 2