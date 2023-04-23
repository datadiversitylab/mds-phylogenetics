
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

singleTreeSize <- function(size, 
                           replicates, 
                           aln_prefix = "alignment",
                           aln_model = "JC", 
                           aln_number = 2,
                           aln_length =  1000,
                           aln_gamma4C = NULL,
                           aln_I = NULL,
                           aln_indel = NULL #"0.03,0.1"
){
  
  
  #This is the true tree
  tree <- rtree(size)
  
  #Let's now simulate the sequence alignment
  system(paste0("iqtree2 --alisim ", aln_prefix, " -m ", aln_model, 
                if(!is.null(aln_I)){paste0("+I{",aln_I,"}")}, 
                if(!is.null(aln_gamma4C)){paste0("+G4{",aln_gamma4C,"}")}, 
                " -t tree.nwk -seed 123 --num-alignments ", aln_number,
                " --length " , aln_length,
                if(!is.null(aln_indel)){paste0(" --indel ",aln_indel)},
                " --no-unaligned"
  ))
  
  tfiles <- list.files(pattern = aln_prefix)
  
  dist <- sapply(tfiles, function(y){
    data <- as.phyDat(read.dna(y))
    
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
  
  return(dist)
}
res <- singleTreeSize(size = 7, replicates = 100)
hist(res$FreqEquivalent, main = "Frequency of the i-MDS \ntrees matching the true tree")




