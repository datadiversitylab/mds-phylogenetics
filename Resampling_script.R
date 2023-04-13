
#import packages
library(ape)
library(phangorn)
library(MASS)
library(smacof)

#simulate data & create distant matrix
set.seed(1)

#Function to create a number of random trees with size `size`. The same 
#function fits an MDS and the inverse-MDS for each of the trees. This function
#also compares the frequency in which inverse-MDS trees match the true tree.

singleTreeSize <- function(size, replicates){
  
  dist <- sapply(1:replicates, function(y){
    #This is the true tree
    tree <- rtree(size)
    
    #Let's now simulate the sequence alignment
    data <- simSeq(tree, l = 500, type = "DNA", bf = c(.1,.2,.3,.4), Q = 1:6, ancestral = FALSE)
    
    #We need to estimate the distance matrix for the estimated alignment
    dm <- dist.ml(data)
    dm[is.na(dm)] <- 0
    
    #performs multidimensional scaling (MDS) & plot results 
    fit_reg <- mds(dm, ndim = 2)
    mds_tree <- as.phylo(hclust(dm))
    
    #performs inverse MDS
    fit_in <- inverseMDS(fit_reg$conf)
    fit_in_trees <- lapply(fit_in, function(x) as.phylo(hclust(x)))
    
    TopDistance <- sapply(fit_in_trees, function(x) dist.topo(unroot(tree), unroot(x)))
    FreqSame <- length(which(TopDistance == 0))/length(TopDistance)
    
    list("FreqEquivalent" =  FreqSame)
  })
  dist <- cbind.data.frame(Tree = paste0("Tree_", 1:replicates), "FreqEquivalent" = unlist(dist))
  return(dist)
}
res <- singleTreeSize(size = 7, replicates = 100)
hist(res$FreqEquivalent, main = "Frequency of the i-MDS \ntrees matching the true tree")




