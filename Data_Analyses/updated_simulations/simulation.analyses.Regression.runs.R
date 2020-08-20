# Simulation analyses

# general phy-signal analyses

library(readxl)
library(latex2exp)
source("simulation.source.R")

#### Pure-birth tree analyses  with Regression ---------------------------------------

### Do this now: set working directory to same location as this file
### Can do that with drop-drown in R Studio

### Parameters

nsim <- 10 # change to 50 for final runs
treesizes <- 2^(3:4) # change to 5:10 for final runs
lambdas <- seq(0, 1, 0.1) # change to seq(0, 1, 0.05) for final runs

### Prepare results

location <- getwd() 
destination <- paste(getwd(), "/results", sep = "")

### Arguments for pure-birth tree analyses

sim.args <- list(treesize = treesizes[1], tree.type = "pb",
                nsim = nsim, lambdas = lambdas, poly.per = 0,
                beta = c(0, .25, .5, .75, 1), 
                mod.type = "reg")

setwd(destination)

### RUNS!  change effects.simulations.PP to effects.simulations below to avoid parallel processing

sim.sets <- lapply(1:length(treesizes), function(j){
  sim.args$treesize = treesizes[j]
  R <- do.call(effects.simulations.PP, sim.args)
  write.csv(R$Result, file = R$file.name, row.names = FALSE)
  cat(j, "tree of", length(treesizes), "trees, completed\n")
  R
})

names(sim.sets) <- treesizes

####  save files as Rdata (optional) / .csv files already saved

nm <- paste(sim.args$tree.type, ".sets.Rdata", sep = "")
save(sim.sets, file = nm)
setwd(location)


### plot lambdas for beta = 0.5

library(readxl)
library(latex2exp)
source("Data_Analyses/updated_simulations/simulation.source.R")
nsim <- 50 # change to 50 for final runs
treesizes <- 2^(5:10) # change to 5:10 for final runs
lambdas <- seq(0, 1, 0.05) # change to seq(0, 1, 0.05) for final runs
sim.sets_32 <- read.csv("Data_Analyses/updated_simulations/results/pb.32.50.sim.0.poly.reg.type.csv")
sim.sets_64 <- read.csv("Data_Analyses/updated_simulations/results/pb.64.50.sim.0.poly.reg.type.csv")
sim.sets_128 <- read.csv("Data_Analyses/updated_simulations/results/pb.128.50.sim.0.poly.reg.type.csv")
sim.sets_256 <- read.csv("Data_Analyses/updated_simulations/results/pb.256.50.sim.0.poly.reg.type.csv")
sim.sets_512 <- read.csv("Data_Analyses/updated_simulations/results/pb.512.50.sim.0.poly.reg.type.csv")
sim.sets_1024 <- read.csv("Data_Analyses/updated_simulations/results/pb.1024.50.sim.0.poly.reg.type.csv")

sim.sets <- list(sim.sets_32, sim.sets_64, sim.sets_128, sim.sets_256, sim.sets_512, sim.sets_1024)
names(sim.sets) <- treesizes

png(filename = "Manuscript/SupportingInfo_Figs/fig.1.reg.png",width = 800, height = 1000, units = "px", pointsize = 18)
par(mfrow = c(3, 2), 
    mar = c(5,5,1,2))
for(i in 1:length(treesizes)) {
  R <- list("Result" = sim.sets[[i]], "file.name" = "filler")
  Rr <- R$Result
  lambda.plot(R, pch = 19, beta = 0.5, col = gray.colors(1, start = 0.6, end = 0.6, alpha = 0.5),
              xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
              ylab = latex2exp::TeX("Estimated Phylogenetic Attributes"),
              ylim = c(0, 1.1))
  
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

### plot kappas for beta = 0.5

png(filename = "Manuscript/SupportingInfo_Figs/fig.2.reg.png",width = 800, 
    height = 1000, units = "px", pointsize = 18)
par(mfrow = c(3,2), 
    mar = c(5,5,1,2))
for(i in 1:length(treesizes)) {
  R <- list("Result" = sim.sets[[i]], "file.name" = "filler")
  Rr <- R$Result
  kappa.plot(R, pch = 19, beta = 0.5, col = gray.colors(1, start = 0.6, end = 0.6, alpha = 0.5),
              xlab = latex2exp::TeX("Input Phylogenetic Signal($\\lambda_{in}$)"),
              ylab = latex2exp::TeX("Estimated Phylogenetic Attributes"),
              ylim = c(0, 2))
  
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

png(filename = "Manuscript/SupportingInfo_Figs/fig.3.reg.png",width = 800, 
    height = 1000, units = "px", pointsize = 18)

par(mfrow = c(3,2), 
    mar = c(5,5,1,2))
zref <- c(1, length(sim.sets))
for(i in 1:2) {
  R <- list("Result" = sim.sets[[i]], "file.name" = "filler")
  maxZ <- max(na.omit(cbind(R$Result$lambda.z, R$Result$kappa.z)))
  minZ <- min(na.omit(cbind(R$Result$lambda.z, R$Result$kappa.z)))
  
  lambda.z.plot(R, pch = 19, col = gray.colors(1, start = 0.6, end = 0.6, alpha = 0.5),
                xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
                ylab = latex2exp::TeX("$Z_{\\lambda}$"),
                mean.col = 1, mean.lwd = 3, ylim = c(-3, 25))
  abline(h = sqrt(qchisq(0.95,1)))
  
  legend("topleft", paste("N =", treesizes[zref[i]]), bty = "n", cex = 1.5)
  
  kappa.z.plot(R, pch = 19, col = gray.colors(1, start = 0.6, end = 0.6, alpha = 0.5),
                xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
                ylab = latex2exp::TeX("$Z_{K}$"),
                mean.col = 1, mean.lwd = 3, ylim = c(-3, 25))
  
  abline(h = sqrt(qchisq(0.95,1)))
  
  legend("topleft", paste("N =", treesizes[zref[i]]), bty = "n", cex = 1.5)
  
}

LZ <- sapply(1:length(sim.sets), function(j) by(sim.sets[[j]]$lambda.z, sim.sets[[j]]$lambda.in, mean))

# Note! change ylim below, if needed.  Too difficult to automate!

plot(lambdas, LZ[,1], ylim = c(-3, 25), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{\\lambda}$"), type = "l")
for(i in 2:ncol(LZ)) points(lambdas, LZ[,i], type = "l", lwd = i)

legend("topleft", as.character(treesizes), lwd = 1:length(treesizes))

# Note! change ylim below, if needed.  Too difficult to automate!

KZ <- sapply(1:length(sim.sets), function(j) by(sim.sets[[j]]$kappa.z, sim.sets[[j]]$lambda.in, mean))
plot(lambdas, KZ[,1], ylim = c(-3, 25), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{K}$"),type = "l")
for(i in 2:ncol(KZ)) points(lambdas, KZ[,i], type = "l", lwd = i)

legend("topleft", as.character(treesizes), lwd = 1:length(treesizes))


dev.off()

### regression precision plot
png(filename = "Manuscript/SupportingInfo_Figs/regression.precision.reg.png",width = 800, 
    height = 1000, units = "px", pointsize = 18)

par(mfrow = c(3,2), 
    mar = c(5,5,1,2))

for(i in 1:length(treesizes)) {
  R <- list("Result" = sim.sets[[i]], "file.name" = "filler")
  Rr <- R$Result
  regression.precision.plot(R, pch = 19, col = gray.colors(1, start = 0.6, end = 0.6, alpha = 0.5),
             xlab = latex2exp::TeX("Input Effect ($\\beta_{in}$)"),
             ylab = latex2exp::TeX("Estimated Slope ($\\beta_{est}$)"),
             ylim = c(0, 2))
  
  legend("topleft", paste("N =", treesizes[i]), bty = "n", cex = 1.5)
}
dev.off()

par(mfcol = c(1,1))


### power curves
library(dplyr)

png(filename = "Manuscript/fig.S19.png", width = 800, 
    height = 1000, units = "px", pointsize = 18)
par(mfrow = c(3,2), 
    mar = c(5,5,1,2))
for(i in 1:length(treesizes)) {
  df <- sim.sets[[i]][,c("beta.in", "lambda.in", "Pval.x")] 
  df$beta.in <- as.factor(df$beta.in)
  df$lambda.in <- as.factor(df$lambda.in)
  new.df <- df %>% group_by(lambda.in, beta.in) %>% summarise(power = length(which(Pval.x < 0.05))/nsim)
  new.df$beta.in <- as.character(new.df$beta.in)
  new.df$beta.in <- as.numeric(new.df$beta.in)
  lambda.cols <- rep(gray.colors(nlevels(new.df$lambda.in), start = 0.2, end = 0.9, alpha = 0.5), each = 5)
  
  plot(new.df$beta.in, new.df$power, pch = 19,  col = lambda.cols,
             xlab = latex2exp::TeX("Input Effect ($\\beta_{in}$)"),
             ylab = latex2exp::TeX("Statistical Power"),
             ylim = c(0, 1))

  for(j in 1:21){ points(new.df$beta.in[(j*5-4):(j*5)], new.df$power[(j*5-4):(j*5)], type = "l", 
                         col = lambda.cols[(j*5-4):(j*5)])}

 legend(title =  latex2exp::TeX("$\\lambda_{in}$"), x = .85, y = .47, title.adj = .5,
        legend = levels(new.df$lambda.in)[c(21, 15, 11, 7, 1)], fill = unique(lambda.cols)[c(21, 15, 11, 7, 1)], cex = .8)
 legend("bottomright", paste("N =", treesizes[i]), bty = "n", cex = 1.5)
}
dev.off()
