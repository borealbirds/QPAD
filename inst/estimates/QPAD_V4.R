ln <- system.file(paste0("estimates/QPAD_v", version, ".rds"), 
                  package="QPAD")
bamcoefs <- readRDS(ln)
.BAMCOEFS <- list2env(bamcoefs)
rm(bamcoefs)
rm(ln)
