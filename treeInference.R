library(ape)
library(phytools)
library(combinat)
library(phangorn)
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/smacof/smacof_1.10-8.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
library(smacof) #1.10-8

#Tree building
size = 10
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
generate_fixed_element_combinations <- function(vector, fixed_element, group_size, num_combinations) {
  # Ensure fixed_element is in the vector
  if (!(fixed_element %in% vector)) {
    stop("The specified fixed element is not in the given vector.")
  }

  # Remove the fixed element from the available choices
  remaining_elements <- setdiff(vector, fixed_element)

  # Check if there are enough elements to form valid groups
  if (length(remaining_elements) < (group_size - 1)) {
    stop("Not enough elements to form the required group size.")
  }

  # Generate num_combinations unique random groups
  combinations <- replicate(num_combinations,
                            c(fixed_element, sample(remaining_elements, group_size - 1, replace = FALSE)),
                            simplify = FALSE)

  # Convert to a matrix for easy viewing
  combination_matrix <- data.frame(do.call(cbind, combinations))

  return(combination_matrix)
}

# For each alignment
# For each target species (100 in this case), find 10 unique quartets
# 20k quartets per species per alignment

dist_trees <- function(alns = alns, group_size = 4, num_combinations = 10, export = TRUE){
  trees_p_aln <- pblapply(alns, function(y){
    data <- as.phyDat(read.dna(y)) # full dataset
    t_taxa <- names(data)
    t_all_t_taxa <-  pblapply(t_taxa, function(z){
      subsetCols <- generate_fixed_element_combinations(vector = t_taxa, fixed_element = z, group_size = group_size, num_combinations = num_combinations)
      trees_ttaxa <- lapply(subsetCols, function(x){
        tryCatch({
        subdata <- data[x,]
        dm <- dist.ml(subdata)
        dm[is.na(dm)] <- 0

        fit_reg <- mds(dm, ndim = 2) #where's the inverse MDS?
        mds_tree <- as.phylo(hclust(dm))
        mds_tree
        }, error=function(e){})
      })
      trees_ttaxa[sapply(trees_ttaxa, is.null)] <- NULL
      trees_ttaxa
    })
    fullTrees <- do.call(c, unlist(t_all_t_taxa, recursive = FALSE))
    if(export == TRUE){
      write.tree(fullTrees, paste0(y, "tree_dist.tre"))
    }
    fullTrees
  })
  return(trees_p_aln)
}

## Reconstruct phylogeny from tree distribution (for each alignment)
alns <- list.files(pattern = "alignment_")
trees_d <- dist_trees(alns = alns[1:10], group_size = 4, num_combinations = 10)

treesR <- list.files(pattern = ".phytree_dist.tre")

lapply(treesR, function(y){
  system(paste0("astral -i ", y, " -o ", y,".out.tre"))
})

##Read in reconstructed trees and compare them to the true topology

library(treeio)

treesRA <- list.files(pattern = ".out.tre")
recTrees <- lapply(treesRA, read.astral)
recTrees <- lapply(recTrees, function(x){
  stree <- x@phylo
  stree$edge.length <- NULL
  stree
  })

class(recTrees) <- "multiPhylo"

write.tree(recTrees, "dist.astral.tree")

tree <- read.tree("tree.nwk")
tree$edge.length <- NULL
treeDist <- lapply(recTrees, function(x){
  disttopo = dist.topo(unroot(tree), unroot(x))
  cbind(raw = disttopo, stn =disttopo/((length(which(tree$edge[,2] > Ntip(tree))) *2))) #Based on the number of internal branches
})

treeDist <-  do.call(rbind, treeDist)

library(ggtree)
library(ggpubr)

p1 <- ggdensitree(recTrees, alpha=.3, colour='steelblue',
            tip.order = tree$tip.label) +
  geom_tiplab(size=3) + hexpand(.35)

p2 <- ggtree(tree, layout="slanted", ladderize = FALSE)+
  geom_tiplab(size=3) + hexpand(.35)

p3 <- ggarrange(p1, p2)

ggsave("comparisonTree.pdf", p3, width = 100, height = 100, units = "mm")

