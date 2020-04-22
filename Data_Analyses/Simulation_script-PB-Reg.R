# Simulating pgls inference errors from estimated lambdas

# Set up
setwd("/Users/ericabaken/Documents/School/Projects/Lambda_Eval/Manuscript/2020-Adams_and_Baken/")

library(phytools)
library(geomorph)
library(geiger)
library(caper)

treesizes <-  2^(5:10) # 32,64,128,256,512,1024
nsim <- 50 # number of simulations
lambdas <- seq(0,1,.05) # starting with just 3 lambda values
lambdas_expanded <- rep(lambdas, each = nsim)
beta <- c(0, .25, .5, .75, 1) # true slope
fullsimlength <- length(lambdas_expanded)

n <- treesizes[4] # number of tips: adjust each iteration

# New Rep #####
DataTable <- array(NA, dim=c(length(beta), 6, length(lambdas)*nsim), 
                   dimnames = list(beta, c("AIC", "F", "P", "lambda.input", "lambda.est.y", "slope")))

DataTable[,4,] <- rep(lambdas_expanded, each = length(beta)) # putting in input lambda

pb <- txtProgressBar(0, fullsimlength, style=1)

for (j in 1:(length(lambdas_expanded))) {

  while(TRUE){
 
    tree<-pbtree(n=n,scale=1) 
    
    if(lambdas_expanded[j] > 0) {
    tree_scaled <- rescale(tree, model = "lambda", lambda = lambdas_expanded[j])
    }
    
    if(lambdas_expanded[j] == 0) { # without this bit, lambda always estimated as 0.5, and it doesn't make sense
    tree_scaled <- rescale(tree, model = "lambda", lambda = 0.01) 
    }
    
    X <- sim.char(tree_scaled, par = 1, nsim = 1)[,,1] 
    Y_correlated <- lapply(1:length(beta), function(i) {X*beta[i] + sim.char(tree_scaled, par = 1, nsim = 1)[,,1]})
        
    comp_data_y <- lapply(1:length(beta), function(i) {
      comparative.data(tree, data.frame(Y=Y_correlated[[i]], X=X, Species = names(Y_correlated[[i]])), Species, vcv = F, vcv.dim = 2) } )

    pgls_results <- try(lapply(1:length(beta), function(i)  pgls(Y ~ X, comp_data_y[[i]], lambda='ML')), silent = T) # estimating lambda using pgls

    if(!is(pgls_results, 'try-error')) break  }

DataTable[,5,j] <- unlist(lapply(1:length(beta), function(i) (pgls_results[[i]]$param[2]))) # estimated lambdas for y
DataTable[,6,j] <- unlist(lapply(1:length(beta), function(i) pgls_results[[i]]$model$coef[2,1])) # slopes
DataTable[,3,j] <- unlist(lapply(1:length(beta), function(i) summary(pgls_results[[i]])$coefficients[2,4])) # p values on slope
DataTable[,2,j] <- unlist(lapply(1:length(beta), function(i) summary(pgls_results[[i]])$fstatistic[[1]])) # f stat
DataTable[,1,j] <- unlist(lapply(1:length(beta), function(i) pgls_results[[i]]$aic)) # AIC

setTxtProgressBar(pb, j)
} 


# Munging data tables  ####
DataTableMat <- rbind(DataTable[,,1], DataTable[,,2])
for (i in 3:length(lambdas_expanded)) {DataTableMat <- rbind(DataTableMat, DataTable[,,i])}
DataTableMat <- as.data.frame(DataTableMat)
DataTableMat$beta <- rep(beta,(length(lambdas)*nsim))

# Test Plots ####


#jpeg("Figures/FigS4.jpeg", res = 60, quality = 100, width = 440, height = 660)
plot(DataTableMat$lambda.est.y ~ DataTableMat$lambda.input, pch = 19, xlab = "Input Lambda", ylab = "Lambda Estimated on Y (reg)")
#dev.off()


# Writing Output Files ####
file_name <- paste("Data_Analyses/Sim_Data/PB_reg-", n, ".csv", sep = "")
write.csv(DataTableMat, file_name, row.names = F)
