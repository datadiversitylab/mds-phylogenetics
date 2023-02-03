#Load the necessary packages 
library(ape) 
library(phangorn)

set.seed(1)

tree <- rtree(5)
data <- simSeq(tree, l = 10, type="DNA", bf=c(.1,.2,.3,.4), Q=1:6,ancestral=TRUE)
dm  <- dist.ml(data)
treeUPGMA  <- upgma(dm)
treeNJ  <- NJ(dm)

treeNJ #Target Tree
tree #True tree
data #Alignment 

# https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html



