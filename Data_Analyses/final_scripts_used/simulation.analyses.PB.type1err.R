# simulation.pb.type1error
library(readxl)
library(latex2exp)
source("simulation.source.R")

# read in simulation data

nsim <- 50 
treesizes <- 2^(5:10) 
lambdas <- seq(0, 1, 0.05) 
sim.sets_32 <- read.csv("results/pb.32.50.sim.0.poly.mean.type.csv")
sim.sets_64 <- read.csv("results/pb.64.50.sim.0.poly.mean.type.csv")
sim.sets_128 <- read.csv("results/pb.128.50.sim.0.poly.mean.type.csv")
sim.sets_256 <- read.csv("results/pb.256.50.sim.0.poly.mean.type.csv")
sim.sets_512 <- read.csv("results/pb.512.50.sim.0.poly.mean.type.csv")
sim.sets_1024 <- read.csv("results/pb.1024.50.sim.0.poly.mean.type.csv")

sim.sets <- list(sim.sets_32, sim.sets_64, sim.sets_128, sim.sets_256, sim.sets_512, sim.sets_1024)
names(sim.sets) <- treesizes

#### compare false discovery rates between Z_Ks calculated from data on *different* random trees ####

type1error.zk <- matrix(NA, nrow = 21, ncol = 6)
colnames(type1error.zk) <-  treesizes
row.names(type1error.zk) <- lambdas

type1error.zl <- matrix(NA, nrow = 21, ncol = 6)
colnames(type1error.zl) <-  treesizes
row.names(type1error.zl) <- lambdas

dist.zk <- matrix(NA, nrow= 50, ncol = 50)
dist.zl <- matrix(NA, nrow= 50, ncol = 50)

for (i in 1:length(treesizes)) { # for each tree size
  for (k in 1:length(lambdas)){ # for each lambda input
    data.i <- sim.sets[[i]][which(as.character(sim.sets[[i]]$lambda.in) == as.character(lambdas[k])),] # extract the 50 sim of that tree size at that lambda input
    for (l in 1:nsim) {
      dist.zk[,l] <- (data.i$kappa.z-data.i$kappa.z[l])/sqrt(2) 
      dist.zl[,l] <- (data.i$lambda.z-data.i$lambda.z[l])/sqrt(2)    
    }
    ps.zk <- pnorm(dist.zk[lower.tri(dist.zk)], lower.tail = F) # calculate p value from Z_12 K stat produced above
    ps.zl <- pnorm(dist.zl[lower.tri(dist.zl)], lower.tail = F) # calculate p value from Z_12 lambda stat produced above
    type1error.zk[k,i] <- length(which(ps.zk < 0.05))/(length(ps.zk)) # rate of p < 0.05 
    type1error.zl[k,i] <- length(which(ps.zl < 0.05))/(length(ps.zl)) # rate of p < 0.05 
  }
}

write.csv(type1error.zk, "results/pb.error.difftrees.zk.csv")
write.csv(type1error.zl, "results/pb.error.difftrees.zl.csv")


# error and false discovery rates over original Fig 3
# same code as what generated original fig 3 except added points, abline, and axis lines, plus edit to legends to add color

type1error.zk <- read.csv("results/pb.error.difftrees.zk.csv")
type1error.zl <- read.csv("results/pb.error.difftrees.zl.csv")


zk.fd.sum <- rowMeans(type1error.zk[,-1])
zk.fd.sum[1] # 0.04680272 <- type 1 error for Z_K on diff trees
mean(zk.fd.sum[-c(1)]) # 0.1185102 <- mean of false discovery rate for Z_K on diff trees

zl.fd.sum <- rowMeans(type1error.zl[,-1]) # same stats for Z_lambda
zl.fd.sum[1] # 0
mean(zl.fd.sum[-c(1)]) # 0.1862857

# Fig 3 with these data overlain
par(mfrow = c(1,2), 
    mar = c(5,5,1,2))

LZ <- sapply(1:length(sim.sets), function(j) by(sim.sets[[j]]$lambda.z, sim.sets[[j]]$lambda.in, mean))

plot(lambdas, LZ[,1], ylim = c(-3, 50), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{\\lambda}$"), type = "l", cex.lab=1.25, cex.axis= 1.25)
for(i in 2:ncol(LZ)) points(lambdas, LZ[,i], type = "l", lwd = i)

points(x = lambdas, y = (type1error.zl[,2])*50, type = "l", lwd = 2, col = "purple") # 32s
points(x = lambdas, y = (type1error.zl[,3])*50, type = "l", lwd = 2, col = "blue") # 64s
points(x = lambdas, y = (type1error.zl[,4])*50, type = "l", lwd = 2, col = "darkgreen") # 128s
points(x = lambdas, y = (type1error.zl[,5])*50, type = "l", lwd = 2, col = "green") # 256s
points(x = lambdas, y = (type1error.zl[,6])*50, type = "l", lwd = 2, col = "yellow") # 512s
points(x = lambdas, y = (type1error.zl[,7])*50, type = "l", lwd = 2, col = "orange") # 1024s
abline(h = 2.5, col = "red") # acceptable type 1 error rates
axis(side = 4, seq(0,50, 5), seq(0,1, .1))

legend("topleft", as.character(treesizes), lwd = 1:length(treesizes), col =  c("purple", "blue", "darkgreen", "green", "yellow", "orange"))

KZ <- sapply(1:length(sim.sets), function(j) by(sim.sets[[j]]$kappa.z, sim.sets[[j]]$lambda.in, mean))
plot(lambdas, KZ[,1], ylim = c(-3, 50), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{K}$"),type = "l", cex.lab=1.25, cex.axis= 1.25)
for(i in 2:ncol(KZ)) points(lambdas, KZ[,i], type = "l", lwd = i)

points(x = lambdas, y = (type1error.zk[,2])*50, type = "l", lwd = 2, col = "purple") # 32s
points(x = lambdas, y = (type1error.zk[,3])*50, type = "l", lwd = 2, col = "blue") # 64s
points(x = lambdas, y = (type1error.zk[,4])*50, type = "l", lwd = 2, col = "darkgreen") # 128s
points(x = lambdas, y = (type1error.zk[,5])*50, type = "l", lwd = 2, col = "green") # 256s
points(x = lambdas, y = (type1error.zk[,6])*50, type = "l", lwd = 2, col = "yellow") # 512s
points(x = lambdas, y = (type1error.zk[,7])*50, type = "l", lwd = 2, col = "orange") # 1024s
abline(h = 2.5, col = "red") # acceptable type 1 error rates
axis(side = 4, seq(0,50, 5), seq(0,1, .1))

legend("topleft", as.character(treesizes), lwd = 1:length(treesizes), col =  c("purple", "blue", "darkgreen", "green", "yellow", "orange"))

par(mfcol = c(1,1))


#### compare false discovery rates between Z_Ks calculated from data on *the same* random trees ####
# this code ran for ~5 days, and I mostly edited the get.stats function from Data_Analyses/updated_simulations/simulation.source.R

nsim.internal <- 50 # number of datasets generated per tree

allresults.zl <- allresults.zk <- array(NA, dim = c(nsim.internal, length(lambdas), length(treesizes)))
dimnames(allresults.zl)[[2]] <- dimnames(allresults.zk)[[2]] <- lambdas
dimnames(allresults.zl)[[3]] <- dimnames(allresults.zk)[[3]] <- treesizes

get.stats.type1 <- function(tree, x, beta, lambda, type, nsim.internal) {
  res_0 <- rep(NA, 2)
  tree_scaled <- tree$tree_scaled
  tree <- tree$tree
  Sig <- matrix(0, 2, 2)
  diag(Sig) <- 1
  resid_all <- if(type == "reg") sim.char(tree_scaled, par = Sig)[,,1] else 
    sim.char(tree_scaled, par = 1, nsim = nsim.internal)
  for(i in 1:nsim.internal){ # this internal loop could be made more efficient if we want
    resid <- as.matrix(resid_all[,,i])
    if(type == "reg") x <- resid[,2]
    mu <- x * beta
    y <- if(type != "mean") mu + resid else resid
    y <- as.matrix(y)
    x <- as.matrix(x)
    resid <- as.matrix(resid)
    rownames(y) <- rownames(x) <- rownames(resid)
    
    df <- data.frame(x = x, y = y, Species = rownames(y))
    df <- comparative.data(tree, df, "Species", vcv = F, vcv.dim = 2, force.root = TRUE)
    
    pgls_est  <- if(type == "mean") try(pgls(y ~ 1, data = df, lambda='ML'), silent = T) else
      try(pgls(y ~ x, data = df, lambda='ML'), silent = T)
    kappa_est <- physignal(cbind(y,0), phy = tree, iter = 999, print.progress = F)
    if(!inherits(pgls_est, "try-error")) {
      summ <- summary(pgls_est)
      Pval <- unlist(summ$coefficients[,4])
      names(Pval) <- paste("Pval", names(Pval), sep = ".")
      #res <- c(unlist(pgls_est$model), Pval, unlist(pgls_est$param.CI)[c(1,6,7)],
      res <- c(lambda.z = lambda.Z(pgls_est), kappa_est["Z"])
      #unlist(kappa_est[c("phy.signal", "Z")]))
      names(res)[length(res)] <- "kappa.z"
    } else res <- pgls_est
    
    res_0 <- rbind(res_0, res)
  }
  res_0[-1,]
}

for (i in 1:length(treesizes)) {
  for (j in 1:length(lambdas)) {
    tree <- make.tree(n = treesizes[i], type = "pb", poly.per = 0, lambda = lambdas[j])
    stats <- get.stats.type1(tree, x = matrix(1, treesizes[i]), 
                             beta = 1, lambda = lambdas[j], type = "mean", 
                             nsim.internal = nsim.internal)
    allresults.zl[,j,i] <- unlist(stats[,1])
    allresults.zk[,j,i] <- unlist(stats[,2])
  }
} 



type1error.zk <- matrix(NA, nrow = length(treesizes), ncol = length(lambdas) )
row.names(type1error.zk) <-  treesizes
colnames(type1error.zk) <- lambdas

type1error.zl <- matrix(NA,  nrow = length(treesizes), ncol = length(lambdas) )
row.names(type1error.zl) <-  treesizes
colnames(type1error.zl) <- lambdas

#
for (i in 1:length(treesizes)){
  for (j in 1:length(lambdas)){
    dist.zk <- matrix(NA, nrow= 50, ncol = 50)
    dist.zl <- matrix(NA, nrow= 50, ncol = 50)
    
    for (l in 1:nsim) {
      dist.zk[,l] <- (as.numeric(allresults.zk[,j,i])-as.numeric(allresults.zk[l,j,i]))/sqrt(2) 
      dist.zl[,l] <- (as.numeric(allresults.zl[,j,i])-as.numeric(allresults.zl[l,j,i]))/sqrt(2) 
    }
    
    ps.zk <- pnorm(dist.zk[lower.tri(dist.zk)], lower.tail = F) # calculate p value from Z_12 K stat produced above
    ps.zl <- pnorm(dist.zl[lower.tri(dist.zl)], lower.tail = F) # calculate p value from Z_12 lambda stat produced above
    
    type1error.zk[i,j] <- length(which(ps.zk < 0.05))/(length(ps.zk))
    type1error.zl[i,j] <- length(which(ps.zl < 0.05))/(length(ps.zl))
    
  }
}


write.csv(type1error.zk, "results/pb.error.sametrees.zk.csv")
write.csv(type1error.zl, "results/pb.error.sametrees.zl.csv")

type1error.zk <- read.csv("results/pb.error.sametrees.zk.csv") # this is in a slightly diff format than the difftrees file equivalent
type1error.zl <- read.csv("results/pb.error.sametrees.zl.csv") # this is in a slightly diff format than the difftrees file equivalent


zk.fd.sum <- colMeans(type1error.zk[,-1])
zk.fd.sum[1] # 0.04544218 <- type 1 error for Z_K on the same trees
mean(zk.fd.sum[-c(1)]) # 0.0690340 <- mean of false discovery rate for Z_K on the same trees

zl.fd.sum <- colMeans(type1error.zl[,-1]) # same data for Z_lambda
zl.fd.sum[1] # 0.00081632
mean(zl.fd.sum[-c(1)]) # 0.1764898



##### new figures 

# all Z_K lambda on 1 plot
type1error.zk <- read.csv("pb.error.sametrees.zk.csv") # this is in a slightly diff format than the difftrees file equivalent
mag <- 11

par(mar = c(5,5,1,5))

KZ <- sapply(1:length(sim.sets), function(j) by(sim.sets[[j]]$kappa.z, sim.sets[[j]]$lambda.in, mean))
plot(lambdas, KZ[,1], ylim = c(0, mag), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{K}$"),type = "l", cex.lab=1.25, cex.axis= 1.25)
for(i in 2:ncol(KZ)) points(lambdas, KZ[,i], type = "l", lwd = i)

points(x = lambdas, y = (type1error.zk[1,-1])*mag, type = "l", lwd = 2, col = "purple") # 32s
points(x = lambdas, y = (type1error.zk[2,-1])*mag, type = "l", lwd = 2, col = "blue") # 64s
points(x = lambdas, y = (type1error.zk[3,-1])*mag, type = "l", lwd = 2, col = "darkgreen") # 128s
points(x = lambdas, y = (type1error.zk[4,-1])*mag, type = "l", lwd = 2, col = "green") # 256s
points(x = lambdas, y = (type1error.zk[5,-1])*mag, type = "l", lwd = 2, col = "yellow") # 512s
points(x = lambdas, y = (type1error.zk[6,-1])*mag, type = "l", lwd = 2, col = "orange") # 1024s
abline(h = mag*.05, col = "red") # acceptable type 1 error rates
axis(side = 4, seq(0, mag, mag/10), seq(0,1, .1))#
mtext(side = 4, "False Discovery Rate", line = 3)

legend("topleft", as.character(treesizes), lwd = 1:length(treesizes), col =  c("purple", "blue", "darkgreen", "green", "yellow", "orange"))


# Just error/false discovery data

png("Manuscript/fig.S21.png", width = 500, 
    height = 600, pointsize = 10)
par(mar = c(5,5,1,3))

plot(lambdas, (type1error.zk[1,-1]), ylim = c(0, 1), lwd = 2, col = "purple",
     ylab = latex2exp::TeX("Type 1 Error/False Discovery Rate"), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     type = "l")

points(x = lambdas, y = (type1error.zk[2,-1]), type = "l", lwd = 2, col = "blue") # 64s
points(x = lambdas, y = (type1error.zk[3,-1]), type = "l", lwd = 2, col = "darkgreen") # 128s
points(x = lambdas, y = (type1error.zk[4,-1]), type = "l", lwd = 2, col = "green") # 256s
points(x = lambdas, y = (type1error.zk[5,-1]), type = "l", lwd = 2, col = "yellow") # 512s
points(x = lambdas, y = (type1error.zk[6,-1]), type = "l", lwd = 2, col = "orange") # 1024s
points(x = lambdas, y = zk.fd.sum, type = "l", lwd = 4, col = "black") # mean false discovery rate
abline(h = .05, col = "red") # acceptable type 1 error rates

legend("topleft", as.character(treesizes), fill =  c("purple", "blue", "darkgreen", "green", "yellow", "orange"))

dev.off()
