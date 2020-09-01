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

### Arguments for pure-birth tree analyses. Supporting information data came from changing these arguments

sim.args <- list(treesize = treesizes[1], tree.type = "pb", # can change tree.type to "bal" for balanced trees or "pect" for pectinate trees
                 nsim = nsim, lambdas = lambdas, poly.per = 0, # can change poly.per to 0.2 for 20% polytomies
                 beta = 1, mod.type = "mean") # can adjust beta to a vector of input betas for investigating ANOVA and regression analyses (e.g., c(0, .25, .5, .75, 1))
# can change mod.type to "anova" or "reg" for investigating ANOVA and regression analyses respectively

setwd(destination)

### RUNS!  change effects.simulations to effects.simulations.PP below if parallel processing is available and preferred

sim.sets <- lapply(1:length(treesizes), function(j){
  sim.args$treesize = treesizes[j]
  R <- do.call(effects.simulations, sim.args)
  write.csv(R$Result, file = R$file.name, row.names = FALSE)
  cat(j, "tree of", length(treesizes), "trees, completed\n")
  R
})

names(sim.sets) <- treesizes

####  save files as Rdata (optional) / .csv files already saved

nm <- paste(sim.args$tree.type, ".sets.Rdata", sep = "")
save(sim.sets, file = nm)
setwd(location)


### plot lambdas

sim.sets_32 <- read.csv("results/pb.32.50.sim.0.poly.mean.type.csv")
sim.sets_64 <- read.csv("results/pb.64.50.sim.0.poly.mean.type.csv")
sim.sets_128 <- read.csv("results/pb.128.50.sim.0.poly.mean.type.csv")
sim.sets_256 <- read.csv("results/pb.256.50.sim.0.poly.mean.type.csv")
sim.sets_512 <- read.csv("results/pb.512.50.sim.0.poly.mean.type.csv")
sim.sets_1024 <- read.csv("results/pb.1024.50.sim.0.poly.mean.type.csv")
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