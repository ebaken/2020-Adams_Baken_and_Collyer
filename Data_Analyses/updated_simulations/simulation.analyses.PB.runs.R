# Simulation analyses

# general phy-signal analyses

library(readxl)
library(latex2exp)
source("simulation.source.R")

#### Pure-birth tree analyses ------------------------------------------

### Do this now: set working directory to same location as this file
### Can do that with drop-drown in R Studio

### Parameters

nsim <- 10 # change to 50 for final runs
treesizes <- 2^(4:7) # change to 5:10 for final runs
lambdas <- seq(0, 1, 0.1) # change to seq(0, 1, 0.05) for final runs

### Prepare results

location <- getwd() 
destination <- paste(getwd(), "/results", sep = "")

### Arguments for pure-birth tree analyses

sim.args <- list(treesize = treesizes[1], tree.type = "pb",
                 nsim = nsim, lambdas = lambdas, poly.per = 0,
                 beta = 1, mod.type = "mean")

setwd(destination)

### RUNS!  change effects.simulations.PP to effects.simulations below to avoid parallel processing

sim.sets <- lapply(1:length(treesizes), function(j){
  sim.args$treesize = treesizes[j]
  R <- do.call(effects.simulations, sim.args)
  #write.csv(R$Result, file = R$file.name, row.names = FALSE)
  cat(j, "tree of", length(treesizes), "trees, completed\n")
  R
})

names(sim.sets) <- treesizes

####  save files as Rdata (optional) / .csv files already saved

nm <- paste(sim.args$tree.type, ".sets.Rdata", sep = "")
save(sim.sets, file = nm)
setwd(location)


### plot lambdas
library(readxl)
library(latex2exp)
source("Data_Analyses/updated_simulations/simulation.source.R")
nsim <- 50 # change to 50 for final runs
treesizes <- 2^(5:10) # change to 5:10 for final runs
lambdas <- seq(0, 1, 0.05) # change to seq(0, 1, 0.05) for final runs
sim.sets_32 <- read.csv("Data_Analyses/updated_simulations/results/pb.32.50.sim.0.poly.mean.type.csv")
sim.sets_64 <- read.csv("Data_Analyses/updated_simulations/results/pb.64.50.sim.0.poly.mean.type.csv")
sim.sets_128 <- read.csv("Data_Analyses/updated_simulations/results/pb.128.50.sim.0.poly.mean.type.csv")
sim.sets_256 <- read.csv("Data_Analyses/updated_simulations/results/pb.256.50.sim.0.poly.mean.type.csv")
sim.sets_512 <- read.csv("Data_Analyses/updated_simulations/results/pb.512.50.sim.0.poly.mean.type.csv")
sim.sets_1024 <- read.csv("Data_Analyses/updated_simulations/results/pb.1024.50.sim.0.poly.mean.type.csv")

sim.sets <- list(sim.sets_32, sim.sets_64, sim.sets_128, sim.sets_256, sim.sets_512, sim.sets_1024)
names(sim.sets) <- treesizes

png(filename = "Manuscript/fig.1.png", width = 800, 
    height = 1000, units = "px", pointsize = 18)
par(mfrow = c(3, 2), 
    mar = c(5,5,1,2))
for(i in 1:length(treesizes)) {
  
  R <- list("Result" = sim.sets[[i]], "file.name" = "filler")
  Rr <- R$Result
  lambda.plot(R, pch = 19,  col = gray.colors(1, start = 0.6, end = 0.6, alpha = 0.5),
              xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
              ylab = latex2exp::TeX("Estimated Phylogenetic Signal Attributes"),
              ylim = c(0, 1.1), cex.lab=1.25, cex.axis= 1.25)
  
  points(lambdas, by(Rr$lambda.opt.lambda, Rr$lambda.in, mean),
         type = "l", lwd = 3)
  abline(0,1)
  
  sds <-  by(Rr$lambda.opt.lambda, Rr$lambda.in, function(x) sd(x))
  
  points(lambdas, sds,
         type = "l", lwd = 3, col = 2)
  
  Ws <- by(Rr$lambda.opt.lambda, Rr$lambda.in, function(x) if(length(unique(x)) == 1) 1 else
    shapiro.test(x)$statistic)
  
  points(lambdas, Ws,
         type = "l", lwd = 3, col = 4)
  
  legend("topleft", paste("N =", treesizes[i]), bty = "n", cex = 1.5)
}
dev.off()


### plot kappas
png(filename = "Manuscript/fig.2.png", width = 800, 
    height = 1000, units = "px", pointsize = 18)
par(mfrow = c(3,2), 
    mar = c(5,5,1,2))
for(i in 1:length(treesizes)) {
  R <- list("Result" = sim.sets[[i]], "file.name" = "filler")
  Rr <- R$Result
  kappa.plot(R, pch = 19,  col = gray.colors(1, start = 0.6, end = 0.6, alpha = 0.5),
             xlab = latex2exp::TeX("Input Phylogenetic Signal($\\lambda_{in}$)"),
             ylab = latex2exp::TeX("Estimated Phylogenetic Signal Attributes"),
             ylim = c(0, 2), cex.lab=1.25, cex.axis= 1.25)
  
  points(lambdas, by(Rr$phy.signal, Rr$lambda.in, mean),
         type = "l", lwd = 3)
  
  sds <-  by(Rr$phy.signal, Rr$lambda.in, function(x) sd(x))
  
  points(lambdas, sds,
         type = "l", lwd = 3, col = 2)
  
  Ws <- by(Rr$phy.signal, Rr$lambda.in, function(x) if(length(unique(x)) == 1) 1 else
    shapiro.test(x)$statistic)
  
  points(lambdas, Ws,
         type = "l", lwd = 3, col = 4)
  
  legend("topleft", paste("N =", treesizes[i]), bty = "n", cex = 1.5)
}
dev.off()

### plot z profiles
### only plot all points for smallest and largest, followed by means trends for all
png(filename = "Manuscript/fig.3.png", width = 800, 
    height = 330, units = "px", pointsize = 10)

par(mfrow = c(1,2), 
    mar = c(5,5,1,2))

LZ <- sapply(1:length(sim.sets), function(j) by(sim.sets[[j]]$lambda.z, sim.sets[[j]]$lambda.in, mean))

# Note! change ylim below, if needed.  Too difficult to automate!

plot(lambdas, LZ[,1], ylim = c(-3, 50), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{\\lambda}$"), type = "l", cex.lab=1.25, cex.axis= 1.25)
for(i in 2:ncol(LZ)) points(lambdas, LZ[,i], type = "l", lwd = i)

legend("topleft", as.character(treesizes), lwd = 1:length(treesizes))

# Note! change ylim below, if needed.  Too difficult to automate!

KZ <- sapply(1:length(sim.sets), function(j) by(sim.sets[[j]]$kappa.z, sim.sets[[j]]$lambda.in, mean))
plot(lambdas, KZ[,1], ylim = c(-3, 50), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{K}$"),type = "l", cex.lab=1.25, cex.axis= 1.25)
for(i in 2:ncol(KZ)) points(lambdas, KZ[,i], type = "l", lwd = i)

legend("topleft", as.character(treesizes), lwd = 1:length(treesizes))

par(mfcol = c(1,1))

dev.off()



png(filename = "Manuscript/fig.S1.png", width = 800, 
    height = 1000, units = "px", pointsize = 18)

par(mfrow = c(3,2), mar = c(5,5,1,2))
zref <- c(1, length(sim.sets))

for(i in 1:2) {
  R <- list("Result" = sim.sets[[zref[i]]], "file.name" = "filler")
  maxZ <- max(na.omit(cbind(R$Result$lambda.z, R$Result$kappa.z)))
  minZ <- min(na.omit(cbind(R$Result$lambda.z, R$Result$kappa.z)))
  
  lambda.z.plot(R, pch = 19,col = gray.colors(1, start = 0.6, end = 0.6, alpha = 0.5),
                xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
                ylab = latex2exp::TeX("$Z_{\\lambda}$"),
                mean.col = 1, mean.lwd = 3, ylim = c(-3, 50), cex.lab=1.25, cex.axis= 1.25)
  abline(h = sqrt(qchisq(0.95,1)))
  
  legend("topleft", paste("N =", treesizes[zref[i]]), bty = "n", cex = 1.5)
  
  kappa.z.plot(R, pch = 19,col = gray.colors(1, start = 0.6, end = 0.6, alpha = 0.5),
               xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
               ylab = latex2exp::TeX("$Z_{K}$"),
               mean.col = 1, mean.lwd = 3, ylim = c(-3, 50), cex.lab=1.25, cex.axis= 1.25)
  
  abline(h = sqrt(qchisq(0.95,1)))
  
  legend("topleft", paste("N =", treesizes[zref[i]]), bty = "n", cex = 1.5)
}
LZ <- sapply(1:length(sim.sets), function(j) by(sim.sets[[j]]$lambda.z, sim.sets[[j]]$lambda.in, mean))

# Note! change ylim below, if needed.  Too difficult to automate!

plot(lambdas, LZ[,1], ylim = c(-3, 50), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{\\lambda}$"), type = "l", cex = 1)
for(i in 2:ncol(LZ)) points(lambdas, LZ[,i], type = "l", lwd = i)

legend("topleft", as.character(treesizes), lwd = 1:length(treesizes))

# Note! change ylim below, if needed.  Too difficult to automate!

KZ <- sapply(1:length(sim.sets), function(j) by(sim.sets[[j]]$kappa.z, sim.sets[[j]]$lambda.in, mean))
plot(lambdas, KZ[,1], ylim = c(-3, 50), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{K}$"),type = "l", cex = 1)
for(i in 2:ncol(KZ)) points(lambdas, KZ[,i], type = "l", lwd = i)

legend("topleft", as.character(treesizes), lwd = 1:length(treesizes))

dev.off()

# Z_K tree comparison
treesizes <- 2^(5:10)
lambdas <- seq(0,1,0.05)
nsim <- 50


treecomp.mat.zk <- matrix(NA, nrow = 21, ncol = 21)
colnames(treecomp.mat.zk) <- lambdas
row.names(treecomp.mat.zk) <- c("32vs32", "32vs64", "32vs128", "32vs256", "32vs512", "32vs1024",
                                "64vs64", "64vs128", "64vs256", "64vs512", "64vs1024",
                                "128vs128", "128vs256", "128vs512", "128vs1024",
                                "256vs256", "256vs512", "256vs1024",
                                "512vs512", "512vs1024",
                                "1024vs1024")
treecomp.mat.zl <- matrix(NA, nrow = 21, ncol = 21)
colnames(treecomp.mat.zl) <- lambdas
row.names(treecomp.mat.zl) <- c("32vs32", "32vs64", "32vs128", "32vs256", "32vs512", "32vs1024",
                                "64vs64", "64vs128", "64vs256", "64vs512", "64vs1024",
                                "128vs128", "128vs256", "128vs512", "128vs1024",
                                "256vs256", "256vs512", "256vs1024",
                                "512vs512", "512vs1024",
                                "1024vs1024")

dist.zk <- matrix(NA, nrow = nsim, ncol = nsim)
dist.zl <- matrix(NA, nrow = nsim, ncol = nsim)

for (i in 1:length(treesizes)){
  #i <- 6
  
  for (j in 1:length(treesizes)) {
    #j <- 6
    if (j>=i){ # only do 1 triangle
      if (i == 1 & j == 1) {row.insert <- 1}
      if (i == 1 & j == 2) {row.insert <- 2}
      if (i == 1 & j == 3) {row.insert <- 3}
      if (i == 1 & j == 4) {row.insert <- 4}
      if (i == 1 & j == 5) {row.insert <- 5}
      if (i == 1 & j == 6) {row.insert <- 6}
      if (i == 2 & j == 2) {row.insert <- 7}
      if (i == 2 & j == 3) {row.insert <- 8}
      if (i == 2 & j == 4) {row.insert <- 9}
      if (i == 2 & j == 5) {row.insert <- 10}
      if (i == 2 & j == 6) {row.insert <- 11}
      if (i == 3 & j == 3) {row.insert <- 12}
      if (i == 3 & j == 4) {row.insert <- 13}
      if (i == 3 & j == 5) {row.insert <- 14}
      if (i == 3 & j == 6) {row.insert <- 15}
      if (i == 3 & j == 4) {row.insert <- 13}
      if (i == 3 & j == 5) {row.insert <- 14}
      if (i == 3 & j == 6) {row.insert <- 15}
      if (i == 4 & j == 4) {row.insert <- 16}
      if (i == 4 & j == 5) {row.insert <- 17}
      if (i == 4 & j == 6) {row.insert <- 18}
      if (i == 5 & j == 5) {row.insert <- 19}
      if (i == 5 & j == 6) {row.insert <- 20}
      if (i == 6 & j == 6) {row.insert <- 21}
      
      for (k in 1:length(lambdas)){
        #k <- 21
        data.i <- sim.sets[[i]][which(as.character(sim.sets[[i]]$lambda.in) == as.character(lambdas[k])),]
        data.j <- sim.sets[[j]][which(as.character(sim.sets[[j]]$lambda.in) == as.character(lambdas[k])),]
        
        for (l in 1:nsim) {
          dist.zk[,l] <- abs(data.j$kappa.z-data.i$kappa.z[l])/sqrt(2) 
          dist.zl[,l] <- abs(data.j$lambda.z-data.i$lambda.z[l])/sqrt(2) 
        }
        ps.zk <- pnorm(dist.zk, lower.tail = F)
        ps.zl <- pnorm(dist.zl, lower.tail = F)
        treecomp.mat.zk[row.insert,k] <- length(which(ps.zk < 0.05) == "TRUE")/(length(ps.zk))
        treecomp.mat.zl[row.insert,k] <- length(which(ps.zl < 0.05) == "TRUE")/(length(ps.zl))
      }
    }
  }
}

#type1error.zk <- matrix(NA, nrow = 21, ncol = 6)
#colnames(type1error.zk) <-  treesizes
#row.names(type1error.zk) <- lambdas
#
#type1error.zl <- matrix(NA, nrow = 21, ncol = 6)
#colnames(type1error.zl) <-  treesizes
#row.names(type1error.zl) <- lambdas
#
#for (i in 1:6) {
#  for (k in 1:length(lambdas)){
#    data.i <- sim.sets[[i]][which(as.character(sim.sets[[i]]$lambda.in) == as.character(lambdas[k])),]
#    dist.zk <- dist(data.i$kappa.z)/sqrt(2)
#    dist.zl <- dist(data.i$lambda.z)/sqrt(2)
#    ps.zk <- pnorm(dist.zk, lower.tail = F)
#    ps.zl <- pnorm(dist.zl, lower.tail = F)
#    type1error.zk[k,i] <- length(which(ps.zk < 0.05) == "TRUE")/(length(ps.zk))
#    type1error.zl[k,i] <- length(which(ps.zl < 0.05) == "TRUE")/(length(ps.zl))
#  }
#}
#
#i <- 1
#data.i <- sim.sets[[i]][which(as.character(sim.sets[[i]]$lambda.in) == as.character(lambdas[k])),]
#hist(data.i$kappa.z, breaks = 20)    
#pnorm(data.i$kappa.z)
#
#dist <- dist(data.i$kappa.z)/sqrt(2)
#length(which(pnorm(dist, lower.tail = F)<.05))/(length(dist))

# new fig 3

par(mfrow = c(1,2), 
    mar = c(5,5,1,2))

LZ <- sapply(1:length(sim.sets), function(j) by(sim.sets[[j]]$lambda.z, sim.sets[[j]]$lambda.in, mean))

# Note! change ylim below, if needed.  Too difficult to automate!

plot(lambdas, LZ[,1], ylim = c(-3, 50), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{\\lambda}$"), type = "l", cex.lab=1.25, cex.axis= 1.25)
for(i in 2:ncol(LZ)) points(lambdas, LZ[,i], type = "l", lwd = i)

points(x = lambdas, y = (treecomp.mat.zl[1,])*50,  type = "l", lwd = 3, col = "purple") # 32s
points(x = lambdas, y = (treecomp.mat.zl[7,])*50,  type = "l", lwd = 3, col = "blue") # 64s
points(x = lambdas, y = (treecomp.mat.zl[12,])*50, type = "l", lwd = 3, col = "darkgreen") # 128s
points(x = lambdas, y = (treecomp.mat.zl[16,])*50, type = "l", lwd = 3, col = "green") # 256s
points(x = lambdas, y = (treecomp.mat.zl[19,])*50, type = "l", lwd = 3, col = "yellow") # 512s
points(x = lambdas, y = (treecomp.mat.zl[21,])*50, type = "l", lwd = 3, col = "orange") # 1024s
abline(h = 2.5, col = "red") # acceptable type 1 error rates
axis(side = 4, seq(0,50, 5), seq(0,1, .1))
#mtext("Type 1 Error Rate", side = 4)

legend("topleft", as.character(treesizes), lwd = 1:length(treesizes), col =  c("purple", "blue", "darkgreen", "green", "yellow", "orange"))

# Note! change ylim below, if needed.  Too difficult to automate!

KZ <- sapply(1:length(sim.sets), function(j) by(sim.sets[[j]]$kappa.z, sim.sets[[j]]$lambda.in, mean))
plot(lambdas, KZ[,1], ylim = c(-3, 50), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{K}$"),type = "l", cex.lab=1.25, cex.axis= 1.25)
for(i in 2:ncol(KZ)) points(lambdas, KZ[,i], type = "l", lwd = i)

points(x = lambdas, y = (treecomp.mat.zk[1,])*50,  type = "l", lwd = 3, col = "purple") # 32s
points(x = lambdas, y = (treecomp.mat.zk[7,])*50,  type = "l", lwd = 3, col = "blue") # 64s
points(x = lambdas, y = (treecomp.mat.zk[12,])*50, type = "l", lwd = 3, col = "darkgreen") # 128s
points(x = lambdas, y = (treecomp.mat.zk[16,])*50, type = "l", lwd = 3, col = "green") # 256s
points(x = lambdas, y = (treecomp.mat.zk[19,])*50, type = "l", lwd = 3, col = "yellow") # 512s
points(x = lambdas, y = (treecomp.mat.zk[21,])*50, type = "l", lwd = 3, col = "orange") # 1024s
abline(h = 2.5, col = "red") # acceptable type 1 error rates
axis(side = 4, seq(0,50, 5), seq(0,1, .1))


legend("topleft", as.character(treesizes), lwd = 1:length(treesizes), col =  c("purple", "blue", "darkgreen", "green", "yellow", "orange"))

par(mfcol = c(1,1))