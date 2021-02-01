#Tests of Z_K and effect of N and P (as done with Z_CR, Z_PLS)
library(geiger)
library(phytools)

#1: Dependency on p
nsims <- 100
p <- c(1,seq(5,25,by=5))
n = treesizes <- 2^(6) 
lambdas <- c(0,0.2, 0.4, 0.6, 0.8, 1)
kapp.max <- 400


tree <- lapply(1:nsims, function(j) {
  kapp <- 1000
  while(kapp > kapp.max) {
    tree <- pbtree(n=n, scale = 1)
    Cov <- vcv.phylo(tree)
    kapp <- kappa(Cov)
  }
  tree
})
#    tree <- pbtree(n=treesizes, nsim = nsims,scale = 1)
results <- lapply(1:length(lambdas), function(k){
  K.p <- simplify2array(lapply(1:length(p), function(i) {
    K.mult <- unlist(lapply(1:nsims, function(j){
      tree_scaled <- rescale(tree[[j]], model = "lambda",  lambda = lambdas[k])
      sigma <- diag(1,p[i])
      e <- sim.char(phy = tree_scaled, par = sigma, nsim = 1)[,,1]
      kappa_est <- physignal(e, phy = tree[[j]], print.progress = F)$Z
    }))
  }))
})

output <- simplify2array(lapply(1:length(lambdas), function(j) {
  mns <- apply(results[[j]],2,mean)
  sds <- apply(results[[j]],2,sd)
  simplify2array(list(mns,sds))
}))
dimnames(output)[[1]] <- p
dimnames(output)[[2]] <- c("mean", "sd")
dimnames(output)[[3]] <- lambdas

save(output, file = "PResults.Rdata")

################################################
#2: Dependency on N
nsims <- 100
p <- 1 #c(1,seq(5,25,by=5))
treesizes <- 2^(4:8) 
lambdas <- c(0,0.2, 0.4, 0.6, 0.8, 1)

resultsN <- lapply(1:length(lambdas), function(k){
  K.p <- simplify2array(lapply(1:length(treesizes), function(i) {
    tree <- lapply(1:nsims, function(j) {
      kapp <- 1000
      while(kapp > kapp.max) {
        tree <- pbtree(n=treesizes[i], scale = 1)
        Cov <- vcv.phylo(tree)
        kapp <- kappa(Cov)
      }
      tree
    })
    K.mult <- unlist(lapply(1:nsims, function(j){
      tree_scaled <- rescale(tree[[j]], model = "lambda",  lambda = lambdas[k])
      sigma <- diag(1,p)
      e <- sim.char(phy = tree_scaled, par = sigma, nsim = 1)[,,1]
      kappa_est <- physignal(e, phy = tree[[j]], print.progress = F)$Z
    }))
  }))
})

outputN <- simplify2array(lapply(1:length(lambdas), function(j) {
  mns <- apply(resultsN[[j]],2,mean)
  sds <- apply(results[[j]],2,sd)
  simplify2array(list(mns,sds))
}))
dimnames(outputN)[[1]] <- treesizes
dimnames(outputN)[[2]] <- c("mean", "sd")
dimnames(outputN)[[3]] <- lambdas

save(outputN, file = "NResults.Rdata")

############## Plots
load("NResults.Rdata")
    outputN <- output
load("PResults.Rdata")


pdf("Z-by-N and p.pdf")
par (mfrow = c(2,2))
p <- c(1,seq(5,25,by=5))
n = treesizes <- 2^(6) 
lambdas <- c(0,0.2, 0.4, 0.6, 0.8, 1)

plot(p,output[,,1][,1],type="l",lwd=2,ylim = c(-3,3),ylab = "Z",
     main = "Null Pattern (Lambda = 0)")
lines(p, output[,,1][,1]-(1.96*output[,,1][,2]),lty=2)
lines(p, output[,,1][,1]+(1.96*output[,,1][,2]),lty=2)


plot(p,output[,,1][,1],type="l",ylim=c(-1,20), ylab = "Z", col=1,
     main = "Patterns across different Lambdas")
lapply(2:length(lambdas), function(j){
  lines(p,output[,,j][,1], col=j)    
})
legend("topright", legend = paste("lambda = ",lambdas), cex = 0.5,lwd=2, lty = 1, col = 1:length(lambdas))

### by N
p <- 1 #c(1,seq(5,25,by=5))
treesizes <- 2^(4:8) 
lambdas <- c(0,0.2, 0.4, 0.6, 0.8, 1)

plot(treesizes,outputN[,,1][,1],type="l",ylim=c(-1,20), ylab = "Z", col=1,
     main = "Patterns across different Lambdas")
lapply(2:length(lambdas), function(j){
  lines(treesizes,outputN[,,j][,1], col=j)    
})
legend("topleft", legend = paste("lambda = ",lambdas), cex = 0.5,lwd=2, lty = 1, col = 1:length(lambdas))


par (mfrow = c(1,1))
dev.off()

