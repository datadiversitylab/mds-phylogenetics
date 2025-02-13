library(ape)
library(phytools)
library(combinat)
library(phangorn)

#Tree building
size = 20
replicates = 100
aln_prefix = "alignment"
aln_model = "GTR"
aln_number = 100
aln_length =  10000
aln_gamma4C = NULL
aln_I = NULL
aln_indel = NULL

## Create true tree (50 species)

tree <- rtree(size)
write.tree(tree, "tree.nwk")

## Generate alignment with parameters

file.remove(list.files(pattern = "alignment_"))

system(paste0("iqtree2 --alisim ", aln_prefix, " -m ", aln_model,
              if(!is.null(aln_I)){paste0("+I{",aln_I,"}")},
              if(!is.null(aln_gamma4C)){paste0("+G4{",aln_gamma4C,"}")},
              " -t tree.nwk -seed 123 --num-alignments ", aln_number,
              " --length " , aln_length,
              if(!is.null(aln_indel)){paste0(" --indel ",aln_indel)},
              " --no-unaligned"
))

## Partition the alignment to run on 4 species each time
dist_trees <- function(alns = alns, sfraction = 0.5, export = TRUE){

  trees_p_aln <- pblapply(alns, function(y){

    data <- as.phyDat(read.dna(y)) # full dataset
    testComb <- combn(names(data), 4) # find all combinations (4 terminals)
    ## Create trees for each alignment (MDS etc)
    subsetCols <- sample(1:ncol(testComb), round(sfraction*ncol(testComb)) )

    trees_4 <- lapply(subsetCols, function(x){
      tryCatch({
        subdata <- data[testComb[,x]]
        dm <- dist.ml(subdata)
        dm[is.na(dm)] <- 0

        fit_reg <- mds(dm, ndim = 2)
        mds_tree <- as.phylo(hclust(dm))
        mds_tree
      }, error=function(e){})
    })

    trees_4[sapply(trees_4, is.null)] <- NULL
    trees_4
    class(trees_4) <- "multiPhylo"
    if(export == TRUE){
    write.tree(trees_4, paste0(y, "tree_dist.tre"))
    }
    trees_4
  })

  return(trees_p_aln)
}

## Reconstruct phylogeny from tree distribution
alns <- list.files(pattern = "alignment_")
trees_d <- dist_trees(alns = alns[1:10], sfraction = 1, export = TRUE)

treesR <- list.files(pattern = ".phytree_dist.tre")

lapply(treesR, function(y){
  system(paste0("astral -i ", y, " -o ", y,".out.tre"))
})

##Read in reconstructed trees and compare them to the true topology

treesRA <- list.files(pattern = ".out.tre")
recTrees <- lapply(treesRA, read.tree)
tree <- read.tree("tree.nwk")
treeDist <- sapply(recTrees, function(x){
dist.topo(unroot(tree), unroot(x))/((length(which(tree$edge[,2] > Ntip(tree))) *2)) #Based on the number of internal branches
})



