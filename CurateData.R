library(data.table)
library(pbapply)
fullResults <- list.files("Results/", full.names = TRUE)
fullResults_names <- list.files("Results/")

results <- pblapply(fullResults, fread)
names(results) <- fullResults_names
results <- rbindlist(results, idcol="Analysis")

results[, c("model", "NA", "AlnSize", "TreeSize") := tstrsplit(Analysis, "_", fixed=TRUE)]

results$model <- sub("model", "", results$model)
results$AlnSize <- sub("length", "", results$AlnSize)
results$TreeSize <- sub("size", "", results$TreeSize )
results$TreeSize  <- sub(".csv", "", results$TreeSize , fixed = TRUE)

results$TreeSize <- as.numeric(results$TreeSize)
results$AlnSize <- as.numeric(results$AlnSize)
results$internalNodes <- results$TreeSize - 1

fwrite(results, "Results/complete.csv")


##
library(dplyr)
sumResults <- results %>%
  group_by(model, AlnSize, TreeSize) %>%
  summarise_at(.vars = c("x"), .funs = list(mean = mean, sd = sd, min = min, max = max))

fwrite(sumResults, "Results/sumResults.csv")






