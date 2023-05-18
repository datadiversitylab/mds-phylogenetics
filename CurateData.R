library(data.table)
library(pbapply)
library(ggplot2)

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
sumResultsCounts <- results %>% 
  group_by(model, AlnSize, TreeSize, x) %>%
  summarise(count = n()) %>%
  mutate(prop = count/sum(count))

fwrite(sumResultsCounts, "Results/sumResults.csv")

sumResultsCounts <- fread("Results/sumResults.csv")

p <- ggplot(data = sumResultsCounts) +
  geom_point(aes(x = TreeSize, y = prop, color = ifelse(x == 0, 0, 1))) +
  ylim(0, 1) +
  facet_wrap(.~model + AlnSize, ncol = 4) +
  theme(legend.position = "none")

pdf("Fig1.pdf", 10, 5)
print(p)
dev.off()




