# Pectinate

# Set up
# setwd("/Users/ericabaken/Documents/School/Projects/Lambda_Eval/Manuscript/2020-Adams_and_Baken/")

library(phytools)
library(geomorph)
library(geiger)
library(caper)

treesizes <-  2^(5:10) # 32,64,128,256,512,1024
nsim <- 50 # number of simulations
lambdas <- seq(0,1,.05) # starting with just 3 lambda values
lambdas_expanded <- rep(lambdas, each = nsim)
fullsimlength <- length(lambdas_expanded)

# Making all the trees ####
library(TreeSimGM)

asymmetric_tree <- sim.age(age=63, numbsim=1, waitsp=function()2, waitext=function()2.5, symmetric=FALSE)
length(asymmetric_tree[[1]]$tip.label)
tree_ladder_32 <- force.ultrametric(asymmetric_tree[[1]], method = "extend")
for (i in 5:12) {tree_ladder_32[[5]] <- NULL}

asymmetric_tree <- sim.age(age=127, numbsim=1, waitsp=function()2, waitext=function()2.5, symmetric=FALSE)
length(asymmetric_tree[[1]]$tip.label)
tree_ladder_64 <- force.ultrametric(asymmetric_tree[[1]], method = "extend")
for (i in 5:12) {tree_ladder_64[[5]] <- NULL}

asymmetric_tree <- sim.age(age=255, numbsim=1, waitsp=function()2, waitext=function()2.5, symmetric=FALSE)
length(asymmetric_tree[[1]]$tip.label)
tree_ladder_128 <- force.ultrametric(asymmetric_tree[[1]], method = "extend")
for (i in 5:12) {tree_ladder_128[[5]] <- NULL}

asymmetric_tree <- sim.age(age=511, numbsim=1, waitsp=function()2, waitext=function()2.5, symmetric=FALSE)
length(asymmetric_tree[[1]]$tip.label)
tree_ladder_256 <- force.ultrametric(asymmetric_tree[[1]], method = "extend")
for (i in 5:12) {tree_ladder_256[[5]] <- NULL}

asymmetric_tree <- sim.age(age=1023, numbsim=1, waitsp=function()2, waitext=function()2.5, symmetric=FALSE)
length(asymmetric_tree[[1]]$tip.label)
tree_ladder_512 <- force.ultrametric(asymmetric_tree[[1]], method = "extend")
for (i in 5:12) {tree_ladder_512[[5]] <- NULL}

asymmetric_tree <- sim.age(age=2048, numbsim=1, waitsp=function()2, waitext=function()2.5, symmetric=FALSE)
length(asymmetric_tree[[1]]$tip.label)
tree_ladder_1024 <- force.ultrametric(asymmetric_tree[[1]], method = "extend")
for (i in 5:12) {tree_ladder_1024[[5]] <- NULL}

phy.list <- list(tree_ladder_32, tree_ladder_64, tree_ladder_128, tree_ladder_256, tree_ladder_512, tree_ladder_1024)

# Change tree size ####

tree<- phy.list[[6]] # pectinate
n <- treesizes[6]

# New Rep #####
DataTable_lambda <- data.frame(lambda.input = lambdas_expanded, 
                               lambda.est = rep(NA, nsim*length(lambdas)), 
                               kappa = rep(NA, nsim*length(lambdas)),
                               kappa.z = rep(NA, nsim*length(lambdas)),
                               CI.lower = rep(NA, nsim*length(lambdas)),
                               CI.upper = rep(NA, nsim*length(lambdas))) 

pb <- txtProgressBar(153, fullsimlength, style=1)

Start.time<-Sys.time()
Start.time
for (j in 153:fullsimlength) {
  while(TRUE){

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

#DataTable_lambda[which(is.na(DataTable_lambda$kappa.z)=="TRUE"),]
#length(which(is.na(DataTable_lambda$kappa.z)=="TRUE"))

# Test Plot ####

plot(DataTable_lambda$lambda.est~DataTable_lambda$lambda.input, pch = 19)

# Writing Output Files ####
file_name <- paste("Data_Analyses/Sim_Data/Pect_lambda_kappacomparison_", n, "-153andUp.csv", sep = "")
write.csv(DataTable_lambda, file_name, row.names = F)
