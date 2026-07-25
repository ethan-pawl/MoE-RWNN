cvres_dir <- file.path("5_competitors", "50_sim", "503_mixdistreg", "cvscores")

job_grid <- expand.grid( 
  list( 
    ialpha = 1:10, 
    ibeta = 1:10, 
    ifold = 1:5
  )
)

jobs_failed <- logical(nrow(job_grid))
for(ijob in 1:nrow(job_grid)) {
  ialpha <- job_grid[ijob,"ialpha"]
  ibeta <- job_grid[ijob,"ibeta"]
  ifold <- job_grid[ijob,"ifold"]
  fname <- paste(ialpha, ibeta, ifold, "cvscore.Rdata", sep = "-")
  if(!file.exists(file.path(cvres_dir, fname))) jobs_failed[ijob] <- TRUE
}

which(jobs_failed) |> 
  paste(collapse = ",") |> 
  cat("\n")