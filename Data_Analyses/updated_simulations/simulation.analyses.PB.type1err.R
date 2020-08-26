# simulation.pb.type1error

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

type1error.zk <- matrix(NA, nrow = 21, ncol = 6)
colnames(type1error.zk) <-  treesizes
row.names(type1error.zk) <- lambdas

type1error.zl <- matrix(NA, nrow = 21, ncol = 6)
colnames(type1error.zl) <-  treesizes
row.names(type1error.zl) <- lambdas

dist.zk <- matrix(NA, nrow= 50, ncol = 50)
dist.zl <- matrix(NA, nrow= 50, ncol = 50)

for (i in 1:length(treesizes)) { # for each tree size
  # i <- 1
  for (k in 1:length(lambdas)){ # for each lambda input
    # k <- 20
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

# new fig 3
# same code as what generated original fig 3 except lines 52:59 and 69:76, plus edit to legends to add color
png("Data_Analyses/updated_simulations/fig3.type1err.PB.png", width = 800, 
    height = 600, units = "px", pointsize = 10)
par(mfrow = c(1,2), 
    mar = c(5,5,1,2))

LZ <- sapply(1:length(sim.sets), function(j) by(sim.sets[[j]]$lambda.z, sim.sets[[j]]$lambda.in, mean))

plot(lambdas, LZ[,1], ylim = c(-3, 50), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{\\lambda}$"), type = "l", cex.lab=1.25, cex.axis= 1.25)
for(i in 2:ncol(LZ)) points(lambdas, LZ[,i], type = "l", lwd = i)

points(x = lambdas, y = (type1error.zl[,1])*50, type = "l", lwd = 2, col = "purple") # 32s
points(x = lambdas, y = (type1error.zl[,2])*50, type = "l", lwd = 2, col = "blue") # 64s
points(x = lambdas, y = (type1error.zl[,3])*50, type = "l", lwd = 2, col = "darkgreen") # 128s
points(x = lambdas, y = (type1error.zl[,4])*50, type = "l", lwd = 2, col = "green") # 256s
points(x = lambdas, y = (type1error.zl[,5])*50, type = "l", lwd = 2, col = "yellow") # 512s
points(x = lambdas, y = (type1error.zl[,6])*50, type = "l", lwd = 2, col = "orange") # 1024s
abline(h = 2.5, col = "red") # acceptable type 1 error rates
axis(side = 4, seq(0,50, 5), seq(0,1, .1))

legend("topleft", as.character(treesizes), lwd = 1:length(treesizes), col =  c("purple", "blue", "darkgreen", "green", "yellow", "orange"))

KZ <- sapply(1:length(sim.sets), function(j) by(sim.sets[[j]]$kappa.z, sim.sets[[j]]$lambda.in, mean))
plot(lambdas, KZ[,1], ylim = c(-3, 50), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{K}$"),type = "l", cex.lab=1.25, cex.axis= 1.25)
for(i in 2:ncol(KZ)) points(lambdas, KZ[,i], type = "l", lwd = i)

points(x = lambdas, y = (type1error.zk[,1])*50, type = "l", lwd = 2, col = "purple") # 32s
points(x = lambdas, y = (type1error.zk[,2])*50, type = "l", lwd = 2, col = "blue") # 64s
points(x = lambdas, y = (type1error.zk[,3])*50, type = "l", lwd = 2, col = "darkgreen") # 128s
points(x = lambdas, y = (type1error.zk[,4])*50, type = "l", lwd = 2, col = "green") # 256s
points(x = lambdas, y = (type1error.zk[,5])*50, type = "l", lwd = 2, col = "yellow") # 512s
points(x = lambdas, y = (type1error.zk[,6])*50, type = "l", lwd = 2, col = "orange") # 1024s
abline(h = 2.5, col = "red") # acceptable type 1 error rates
axis(side = 4, seq(0,50, 5), seq(0,1, .1))

legend("topleft", as.character(treesizes), lwd = 1:length(treesizes), col =  c("purple", "blue", "darkgreen", "green", "yellow", "orange"))

par(mfcol = c(1,1))

dev.off()
