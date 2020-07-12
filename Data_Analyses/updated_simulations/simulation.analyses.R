# Simulation analyses

# general phy-signal analyses

library(readxl)
library(latex2exp)

#### pure-birth tree analyses
source("simulation.source.R")
nsim <- 50
treesizes <- 2^(5:8)
lambdas <- seq(0, 1, 0.05)
location <- getwd() # should be same as this file location
destination <- paste(getwd(), "/results", sep = "")

pb.args <- list(treesize = treesizes[1], tree.type = "pb",
                nsim = nsim, lambdas = lambdas, poly.per = 0,
                beta = 1, mod.type = "mean")

setwd(destination)

pb.sets <- lapply(1:length(treesizes), function(j){
  pb.args$treesize = treesizes[j]
  R <- do.call(effects.simulations.PP, pb.args)
  write.csv(R$Result, file = R$file.name, row.names = FALSE)
  cat(j, "tree of", length(treesizes), "trees, completed\n")
  R
})

names(pb.sets) <- treesizes

# save files

save(pb.sets, file = "pb.sets.Rdata")
setwd(location)


# plot lambdas
par(mfrow = c(2, length(treesizes)/2))
for(i in 1:length(treesizes)) {
  R <- pb.sets[[i]]
  Rr <- R$Result
  lambda.plot(R, pch = 19, col = gray.colors(1, start = 0.6, end = 0.6, alpha = 0.5),
              xlab = latex2exp::TeX("Input Phylogenetic Signal($\\lambda_{in}$)"),
              ylab = latex2exp::TeX("Estimated Phylogenetic Signal ($\\hat{lambda}$), $s$, or $W$"),
              ylim = c(0, 1.1))
  
  points(lambdas, by(Rr$lambda.opt.lambda, Rr$lambda.in, mean),
         type = "l", lwd = 3)
  abline(0,1)
  
  sds <-  by(Rr$lambda.opt.lambda, Rr$lambda.in, function(x) sd(x))
  
  points(lambdas, sds,
         type = "l", lwd = 3, col = 2)
  
  
  Ws <- by(Rr$lambda.opt.lambda, Rr$lambda.in, function(x) shapiro.test(x)$statistic)
  
  points(lambdas, Ws,
         type = "l", lwd = 3, col = 4)
  
  legend("topleft", paste("N =", treesizes[i]), bty = "n", cex = 1.5)
}

# plot kappas
par(mfrow = c(2, length(treesizes)/2))
for(i in 1:length(treesizes)) {
  R <- pb.sets[[i]]
  Rr <- R$Result
  kappa.plot(R, pch = 19, col = gray.colors(1, start = 0.6, end = 0.6, alpha = 0.5),
              xlab = latex2exp::TeX("Input Phylogenetic Signal($\\lambda_{in}$)"),
              ylab = latex2exp::TeX("Estimated Phylogenetic Signal ($\\kappa$),, $s$, or $W$"),
              ylim = c(0, 2))
  
  points(lambdas, by(Rr$phy.signal, Rr$lambda.in, mean),
         type = "l", lwd = 3)
  
  sds <-  by(Rr$phy.signal, Rr$lambda.in, function(x) sd(x))
  
  points(lambdas, sds,
         type = "l", lwd = 3, col = 2)
  
  Ws <- by(Rr$phy.signal, Rr$lambda.in, function(x) shapiro.test(x)$statistic)
  
  points(lambdas, Ws,
         type = "l", lwd = 3, col = 4)
  
  legend("topleft", paste("N =", treesizes[i]), bty = "n", cex = 1.5)
}



par(mfcol = c(2, length(treesizes) + 1))
zref <- 1:length(pb.sets)
for(i in seq_len(length(zref))) {
  R <- pb.sets[[zref[i]]]
  maxZ <- max(na.omit(cbind(R$Result$lambda.z, R$Result$kappa.z)))
  minZ <- min(na.omit(cbind(R$Result$lambda.z, R$Result$kappa.z)))
  
  lambda.z.plot(R, pch = 19,col = gray.colors(1, start = 0.6, end = 0.6, alpha = 0.5),
                xlab = latex2exp::TeX("Input Phylogenetic Signal($\\lambda_{in}$)"),
                ylab = latex2exp::TeX("$Z_{\\lambda}$"),
                mean.col = 1, mean.lwd = 3, ylim = c(-3, 25))
  abline(h = sqrt(qchisq(0.95,1)))
  


  legend("topleft", paste("N =", treesizes[zref[i]]), bty = "n", cex = 1.5)
  
  kappa.z.plot(R, pch = 19,col = gray.colors(1, start = 0.6, end = 0.6, alpha = 0.5),
                xlab = latex2exp::TeX("Input Phylogenetic Signal($\\lambda_{in}$)"),
                ylab = latex2exp::TeX("$Z_{\\kappa}$"),
                mean.col = 1, mean.lwd = 3, ylim = c(-3, 25))
  
  abline(h = sqrt(qchisq(0.95,1)))
  

  
  legend("topleft", paste("N =", treesizes[zref[i]]), bty = "n", cex = 1.5)
  
}

LZ <- sapply(1:length(pb.sets), function(j) by(pb.sets[[j]]$Result$lambda.z, pb.sets[[j]]$Result$lambda.in, mean))

plot(lambdas, LZ[,1], ylim = c(-3, 25), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{\\lambda}$"), type = "l")
for(i in 2:ncol(LZ)) points(lambdas, LZ[,i], type = "l", lwd = i)

legend("topleft", as.character(treesizes), lwd = 1:length(treesizes))

KZ <- sapply(1:length(pb.sets), function(j) by(pb.sets[[j]]$Result$kappa.z, pb.sets[[j]]$Result$lambda.in, mean))
plot(lambdas, KZ[,1], ylim = c(-3, 25), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("$Z_{\\lambda}$"),type = "l")
for(i in 2:ncol(KZ)) points(lambdas, KZ[,i], type = "l", lwd = i)

legend("topleft", as.character(treesizes), lwd = 1:length(treesizes))




