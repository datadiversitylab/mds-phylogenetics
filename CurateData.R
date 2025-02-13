library(data.table)
library(pbapply)
library(ggplot2)

fullResults <- list.files("Results/", full.names = TRUE)
fullResults_names <- list.files("Results/")

results <- pblapply(fullResults, fread)
names(results) <- fullResults_names
results <- rbindlist(results, idcol="Analysis", fill = TRUE)

results[, c("model", "NA", "AlnSize", "TreeSize") := tstrsplit(Analysis, "_", fixed=TRUE)]

results$model <- sub("model", "", results$model)
results$AlnSize <- sub("length", "", results$AlnSize)
results$TreeSize <- sub("size", "", results$TreeSize )
results$TreeSize  <- sub(".csv", "", results$TreeSize , fixed = TRUE)

results$TreeSize <- as.numeric(results$TreeSize)
results$AlnSize <- as.numeric(results$AlnSize)
results$internalNodes <- results$TreeSize - 1

##Estimating number of int branches
#TR <- rmtopology(1, 10)[[1]]
#(length(which(TR$edge[,2] > Ntip(TR))) *2)
results$internalBranches <- ifelse(results$TreeSize == 4,
                                   1, ifelse(results$TreeSize == 5,
                                             2, 7))

fwrite(results, "Results/complete.csv")

head(results)

##
library(dplyr)
sumResultsCounts <- results %>%
  group_by(model, AlnSize, TreeSize) %>%
  mutate(prop = x/(2*internalBranches)) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop = sd(prop, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(sumResultsCounts, "Results/sumResults.csv")

sumResultsCounts <- fread("Results/sumResults.csv")

p <- ggplot(data = sumResultsCounts) +
  geom_point(aes(x = TreeSize, y = mean_prop, color = model, shape = factor(AlnSize))) +
  ylim(0, 1) +
  #facet_wrap(.~model + AlnSize, ncol = 4) +
  theme(legend.position = "none")


pdf("Fig1.pdf", 10, 5)
print(p)
dev.off()




