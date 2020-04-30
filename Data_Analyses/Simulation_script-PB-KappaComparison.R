# Simulating pgls inference errors from estimated lambdas

# Set up
# setwd("/Users/ericabaken/Documents/School/Projects/Lambda_Eval/Manuscript/2020-Adams_and_Baken/")

library(phytools)
library(geomorph)
library(geiger)
library(caper)

 # extract kappa values for comparison

treesizes <-  2^(5:10) # 32,64,128,256,512,1024
nsim <- 50 # number of simulations
lambdas <- seq(0,1,.05) # starting with just 3 lambda values
lambdas_expanded <- rep(lambdas, each = nsim)
fullsimlength <- length(lambdas_expanded)

n <- treesizes[1] # number of tips: adjust each iteration

# New Rep #####
DataTable_lambda <- data.frame(lambda.input = lambdas_expanded, 
                               lambda.est = rep(NA, nsim*length(lambdas)), 
                               kappa = rep(NA, nsim*length(lambdas)),
                               kappa.z = rep(NA, nsim*length(lambdas)),
                               CI.lower = rep(NA, nsim*length(lambdas)),
                               CI.upper = rep(NA, nsim*length(lambdas))) 

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
    
    X <- sim.char(tree_scaled, par = 1, nsim = 1)[,,1] 
    
    comp_data_x <- comparative.data(tree, data.frame(X.df=X, Species = names(X)), Species, vcv = F, vcv.dim = 2)
    
    pgls_est <- try(pgls(X.df ~ 1, comp_data_x, lambda='ML'),silent = T)
    kappa_est <- try(physignal(X, phy = tree, iter = 999, print.progress = F), silent = T)
    
    if(!is(pgls_est, 'try-error') & !is(kappa_est, 'try-error')) break  }

DataTable_lambda[j,2] <- pgls_est$param[2]
DataTable_lambda[j,3] <- kappa_est$phy.signal 
DataTable_lambda[j,4] <- kappa_est$Z

DataTable_lambda[j,5] <- pgls_est$param.CI$lambda$ci.val[1] # lower CI
DataTable_lambda[j,6] <- pgls_est$param.CI$lambda$ci.val[2] # upper CI

setTxtProgressBar(pb, j)
} 
Sys.time()-Start.time

# fixing the lower CI NAs (which are 0) and the upper CI NAs (which are 1)

DataTable_lambda[which(is.na(DataTable_lambda$CI.lower)=="TRUE"),5] <- 0
DataTable_lambda[which(is.na(DataTable_lambda$CI.upper)=="TRUE"),6] <- 1

anyNA(DataTable_lambda) # want FALSE

DataTable_lambda[which(is.na(DataTable_lambda$kappa.z)=="TRUE"),]
length(which(is.na(DataTable_lambda$kappa.z)=="TRUE"))

# Writing Output Files ####
file_name <- paste("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_", n, ".csv", sep = "")
write.csv(DataTable_lambda, file_name, row.names = F)
