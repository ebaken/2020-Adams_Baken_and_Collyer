# Simulation analyses

# general phy-signal analyses

#### pure-birth tree analyses

nsim <- 100
treesizes <- 2^(5:10)
labmdas <- seq(0, 1, 0.05)
location <- getwd() # should be same as this file location
destination <- paste(getwd(), "/results", sep = "")

pb.args <- list(treesize = treesizes[1], tree.type = "pb",
                nsim = nsim, lambdas = lambdas, poly.per = 0,
                beta = 1, mod.type = "mean")

pb.sets <- lapply(1:length(treesizes), function(j){
  pb.args$treesize = treesizes[j]
  R <- do.call(effects.simulations.PP, pb.args)

  cat(j, "tree of", length(treesizes), "trees, completed\n")
  R
})

names(pb.sets) <- treesizes

# save files
setwd(destination)
lapply(1:length(pb.sets), function(j){
  R <- pb.sets[[j]]
  write.csv(R$Result, file = R$file.name, row.names = FALSE)
})

save(pb.sets, file = "pb.sets.Rdata")
setwd(location)


# plot lambdas
par(mfrow = c(2, length(treesizes)/2))
for(i in 1:length(treesizes)) {
  R <- pb.sets[[i]]
  lambda.plot(R, pch = 19, 
              xlab = latex2exp::TeX("Input Phylogenetic Signal($\\lambda_{in}$)"),
              ylab = latex2exp::TeX("Estimated Phylogenetic Signal ($\\lambda_{est}$)"),
              ylim = c(0, 1.1))
  legend("topleft", paste("N =", treesizes[i]), bty = "n", cex = 1.5)
}



