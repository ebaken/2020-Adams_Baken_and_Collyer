# Simulating pgls inference errors from estimated lambdas

# Set up
# setwd("/Users/ericabaken/Documents/School/Projects/Lambda_Eval/Manuscript/2020-Adams_and_Baken/")

library(phytools)
library(geomorph)
library(geiger)
library(caper)

 # extract kappa values for comparison

treesizes <-  2^(5:10) # 32,64,128,256,512,1024
nsim <- 200 # number of simulations
lambdas <- seq(0,1,.05) # starting with just 3 lambda values
lambdas_expanded <- rep(lambdas, each = nsim)
fullsimlength <- length(lambdas_expanded)

n <- treesizes[2] # number of tips: adjust each iteration

Perc <- 0 # Percentage of the nodes of the given “true” chronogram to be collapsed
Total<- n - 2
To_collapse <- round((Total*Perc)/100)

gg <- function(tree, n){ # From Molina-Venegas and Rodríguez (2017)
  if(n > 0) {
    for(i in 1:n){
      ii <- sample(2:tree$Nnode,1) + length(tree$tip.label)
      ll <- tree$edge.length[which(tree$edge[,2] == ii)]
      tree$edge.length[which(tree$edge[,1] == ii)]<-
        tree$edge.length[which(tree$edge[,1] == ii)] + ll
      tree$edge.length[which(tree$edge[,2] == ii)] <- 0
      tree<-di2multi(tree)
    }
  }
  tree
}



lambda.Z <- function(pgls.mod) {
  FP <- pgls.mod$param
  pars <- pgls.likelihood(optimPar = FP, fixedPar = pgls.mod$param, 
                          y = pgls.mod$y, x = pgls.mod$x, V = pgls.mod$Vt, 
                          optim.output = FALSE)
  FP[2] <- 0
  pars0 <- pgls.likelihood(optimPar = FP, fixedPar = FP, 
                           y = pgls.mod$y, x = pgls.mod$x, V = pgls.mod$Vt, 
                           optim.output = FALSE)
  
  Q0 <- -2*(pars0$ll - pars$ll)
  Z <- sqrt(abs(Q0))

  
  Z
  
}

# New Rep #####
DataTable_lambda <- data.frame(lambda.input = lambdas_expanded, 
                               lambda.est = rep(NA, nsim*length(lambdas)), 
                               kappa = rep(NA, nsim*length(lambdas)),
                               kappa.z = rep(NA, nsim*length(lambdas)),
                               CI.lower = rep(NA, nsim*length(lambdas)),
                               CI.upper = rep(NA, nsim*length(lambdas)),
                               lambda.z0 = rep(NA, nsim*length(lambdas))) 

pb <- txtProgressBar(1, fullsimlength, style=1)

Start.time<-Sys.time()
Start.time
for (j in 1:fullsimlength) {
  while(TRUE){
 
    tree<-pbtree(n=n,scale=1) 
    
    if(lambdas_expanded[j] > 0) {
    tree_scaled <- rescale(tree, model = "lambda", lambda = lambdas_expanded[j])
    }
    
    if(lambdas_expanded[j] == 0) { # without this bit, lambda always estimated as 0.5, and it doesn't make sense
    tree_scaled <- rescale(tree, model = "lambda", lambda = 0.01) 
    }
    
    tree_scaled <- gg(tree_scaled, To_collapse)
    
    X <- sim.char(tree_scaled, par = 1, nsim = 1)[,,1] 
    
    comp_data_x <- comparative.data(tree, data.frame(X.df=X, Species = "y"), Species, vcv = F, vcv.dim = 2)
    
    pgls_est <- try(pgls(X.df ~ 1, comp_data_x, lambda='ML'),silent = T)
    kappa_est <- try(physignal(X, phy = tree, iter = 999, print.progress = F), silent = T)
    
    if(!is(pgls_est, 'try-error') & !is(kappa_est, 'try-error')) break  }

DataTable_lambda[j,2] <- pgls_est$param[2]
DataTable_lambda[j,3] <- kappa_est$phy.signal 
DataTable_lambda[j,4] <- kappa_est$Z

DataTable_lambda[j,5] <- pgls_est$param.CI$lambda$ci.val[1] # lower CI
DataTable_lambda[j,6] <- pgls_est$param.CI$lambda$ci.val[2] # upper CI
DataTable_lambda[j,7] <- lambda.Z(pgls_est)

setTxtProgressBar(pb, j)
} 
Sys.time()-Start.time

# fixing the lower CI NAs (which are 0) and the upper CI NAs (which are 1)

DataTable_lambda[which(is.na(DataTable_lambda$CI.lower)=="TRUE"),5] <- 0
DataTable_lambda[which(is.na(DataTable_lambda$CI.upper)=="TRUE"),6] <- 1

anyNA(DataTable_lambda) # want FALSE

DataTable_lambda[which(is.na(DataTable_lambda$kappa.z)=="TRUE"),]
length(which(is.na(DataTable_lambda$kappa.z)=="TRUE"))

#  data tables  ####
DataTableMat <- DataTable_lambda


# Test Plots ####

DTM <- DataTableMat

par(mfrow = c(1,2))

all.z <- na.omit(c(DTM$lambda.z0, DTM$kappa.z))
Ylim <- c(min(all.z), max(all.z))

plot(DTM$lambda.z0 ~ DTM$lambda.input, pch = 19, ylim = Ylim,
     xlab = "Input Lambda", ylab = "Estimated Lambda Effect Size or Variance of Z ",
     main = "Z0 Lambda")
vars <- by(DTM$lambda.z0, DTM$lambda.input, var)
vars[is.na(sds)] <- 0
points(lambdas, vars, type = "l", col = "red", lwd = 2)



plot(DTM$kappa.z ~ DTM$lambda.input, pch = 19,  ylim = Ylim,
     xlab = "Input Lambda", ylab = "Estimated Kappa Effect Size or Variance of Z",
     main = "Z Kappa")
vars <- by(DTM$kappa.z, DTM$lambda.input, var)
points(lambdas, vars, type = "l", col = "red", lwd = 2)

par(mfrow = c(1, 1))

par(mfrow = c(2, 2))

plot(DTM$lambda.est ~ DTM$lambda.input, pch = 19, 
     xlab = "Input Lambda", ylab = "Lambda Estimated")
plot(DTM$kappa ~ DTM$lambda.input, pch = 19, 
     xlab = "Input Lambda", ylab = "Kappa Estimated")
plot(DTM2$lambda.z ~ DTM2$lambda.input, pch = 19, 
     xlab = "Input Lambda", ylab = "Estimated Lambda Effect Size")
plot(DTM$kappa.z ~ DTM$lambda.input, pch = 19, 
     xlab = "Input Lambda", ylab = "Estimated Kappa Effect Size")

par(mfrow = c(1, 1))

coefvar <- function(x) sd(x)/mean(x) * 100
by(DTM2$lambda.z, DTM2$lambda.input, var)
by(DTM$kappa.z, DTM$lambda.input, var)

# Writing Output Files ####
file_name <- paste("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_for.polytomies_", n, ".csv", sep = "")
write.csv(DataTable_lambda, file_name, row.names = F)
