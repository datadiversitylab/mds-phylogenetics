
#import packages
library(ape)
library(phangorn)
library(MASS)

#packageurl <- "https://cran.r-project.org/src/contrib/Archive/smacof/smacof_1.10-8.tar.gz" 
#install.packages(packageurl, repos=NULL, type="source")
library(smacof) #1.10-8


#simulate data & create distant matrix
set.seed(1)

find_mode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

#Function to create a number of random trees with size `size`. The same 
#function fits an MDS and the inverse-MDS for each of the trees. This function
#also compares the frequency in which inverse-MDS trees match the true tree.

singleTreeSize <- function(size, 
                           replicates, 
                           aln_prefix = "alignment",
                           aln_model = "JC", 
                           aln_number = 100,
                           aln_length =  1000,
                           aln_gamma4C = NULL,
                           aln_I = NULL,
                           aln_indel = NULL #"0.03,0.1"
){
  
  compReps <- lapply(1:replicates, function(x){
  
  #This is the true tree
  tree <- rtree(size)
  write.tree(tree, "tree.nwk")

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
    tryCatch({
    data <- as.phyDat(read.dna(y))
    
    #We need to estimate the distance matrix for the estimated alignment
    dm <- dist.ml(data)
    dm[is.na(dm)] <- 0
    
    #performs multidimensional scaling (MDS) & plot results 
    fit_reg <- mds(dm, ndim = 2)
    mds_tree <- as.phylo(hclust(dm))
    
    #performs inverse MDS
    fit_in <- inverseMDS(fit_reg$conf)
    fit_in_trees <- lapply(fit_in$dissmat, function(x) as.phylo(hclust(x)))
    
    TopDistance <- sapply(fit_in_trees, function(x) dist.topo(unroot(tree), unroot(x)))
    
    list("TopDistance" =  TopDistance)
    }, error=function(e){})
  })
  
  #Clean up the folder
  file.remove(tfiles)
  file.remove("tree.nwk")
  file.remove("tree.nwk.log")
  
  #basic stats
  target <- unlist(dist)
  return(target)
  })
  return(unlist(compReps))
}


params <- expand.grid(model = c("JC", "GTR"),
                      aln_length = c(100, 500, 1000, 10000),
                      size = c(4, 5, 10, 100, 1000)
            )
params$model <- as.character(params$model)

lapply(1:nrow(params), function(x){
  
  partialres <- singleTreeSize(size = params$size[x], 
                 replicates = 1000, 
                 aln_model = params$model[x], 
                 aln_length =  params$aln_length[x])
  
  nameFile <- paste0(names(params), params[x,], collapse = "_")
  write.csv(partialres, paste0(nameFile, ".csv"))
  
})





