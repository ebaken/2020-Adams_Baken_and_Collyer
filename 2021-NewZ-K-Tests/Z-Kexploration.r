#does total length of tree make a difference? 
  #rescale keeps total length the same. That's not it.
  
#something about missing tree size in Z or Z_12 (sqrt(N) often in Z-denominator)

library(geomorph)
library(phytools)
library(geiger)
n=64
tree <- pbtree(n=n)
vcv.phylo(tree)[1,1]  #length of ultrametric tree changes each time

vcv.phylo(rescale(tree,model = "lambda", 0.5))[1,1]

n=256
tree <- pbtree(n=n)
tree.1 <-rescale(tree, model = "depth",1)
y <- fastBM(tree)

#rescaling tree does not affect K, p, or Z_K
physignal(y,tree)$Z/n
physignal(y,tree.1)

## Perhaps the issue is that for Z_K, denominator should account for sqrt(N) 
   #much like t-statistic
treesizes <- 2^(5:10) 
ntrees <- 10

results <- lapply(1:length(treesizes), function(i){ 
  tree <- pbtree(n=treesizes[i], scale = 1, nsim = ntrees)
  out <- lapply(1:ntrees, function(j) {
    y <- fastBM(tree[[j]])
    K <- physignal(y,tree[[j]])
  })
})

 z.k <- simplify2array(lapply(1:length(treesizes), function(j) {
   unlist(simplify2array(results[[j]])[4,])
 }))

plot(t(replicate(treesizes,n = 10)),z.k)
plot(t(replicate(treesizes,n = 10)),z.k/sqrt(treesizes))
