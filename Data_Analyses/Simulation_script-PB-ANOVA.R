# Simulating pgls inference errors from estimated lambdas

# look at Revell ~2009 where he investigated Y~X simulating error in X and in Y or one or the other under BM
setwd("/Users/ericabaken/Documents/School/Projects/Lambda_Eval/Manuscript/2020-Adams_and_Baken/")

library(phytools)
library(geomorph)
library(geiger)
library(svMisc)
library(caper)

treesizes <-  2^(5:10) # 32,64,128,256,512,1024
n <- treesizes[5] # number of tips (not simulations)
nsim <- 50 # number of simulations
lambdas <- seq(0,1,.05) # starting with just 3 lambda values
lambdas_expanded <- rep(lambdas, each = nsim)
beta <- c(0, .25, .5, .75, 1) # true slope

DataTable_lambda_0 <- array(NA, dim=c(length(beta), 6, length(lambdas)*nsim), 
                   dimnames = list(beta, c("Rsq", "F", "P", "lambda.input", "lambda.est", "slope")))

DataTable_lambda_ml <- array(NA, dim=c(length(beta), 6, length(lambdas)*nsim), 
                   dimnames = list(beta, c("Rsq", "F", "P", "lambda.input", "lambda.est", "slope")))

DataTable_lambda_0[,4,] <- rep(lambdas_expanded, each = length(beta)) # putting in input lambda
DataTable_lambda_ml[,4,] <- rep(lambdas_expanded, each = length(beta)) # putting in input lambda

Start.time<-Sys.time()
Start.time
for (j in 1:(length(lambdas_expanded))) {

  while(TRUE){
 
    tree<-pbtree(n=n,scale=1) 
    
    if(lambdas_expanded[j] > 0) {
    tree_scaled <- rescale(tree, model = "lambda", lambda = lambdas_expanded[j])
    }
    
    if(lambdas_expanded[j] == 0) { # without this bit, lambda always estimated as 0.5, and it doesn't make sense
    tree_scaled <- rescale(tree, model = "lambda", lambda = 0.01) 
    }
    
    X <- rTraitDisc(tree_scaled, model = "ER", states = c(0,1))
    #X <- sim.char(trees_scaled, par = 1, nsim = 1)[,,1] # simmap for binary? email from dean, luke harmon
           # We need to simulate a discrete variable:https://lukejharmon.github.io/pcm/chapter7_introdiscrete/
            # Threshold model: http://blog.phytools.org/2013/09/threshold-character-evolution.html
    #comp_data_x <- comparative.data(tree, data.frame(X=X, Species = names(X)), Species, vcv = F, vcv.dim = 2)
    
    #pgls_x_lambda_0 <- try(pgls(X ~ 1, comp_data_x, lambda=0.01), silent = T)
    #pgls_x_lambda_ml <- try(pgls(X ~ 1, comp_data_x, lambda='ML'), silent = T)
    
    Y_correlated <- lapply(1:length(beta), function(i) {as.numeric(X)*beta[i] + rnorm(n = n)})
        
    comp_data_y <- lapply(1:length(beta), function(i) {
      comparative.data(tree, data.frame(Y=Y_correlated[[i]], X=X, Species = names(X)), Species, vcv = F, vcv.dim = 2) 
      } )
     
    pgls_xy_lambda_0 <- try(lapply(1:length(beta), function(i)  pgls(Y ~ X, comp_data_y[[i]], lambda=0.01)), silent = T)
    pgls_xy_lambda_ml <- try(lapply(1:length(beta), function(i)  pgls(Y ~ X, comp_data_y[[i]], lambda='ML')), silent = T) # estimating lambda using pgls

    if(#!is(pgls_x_lambda_0, 'try-error') & !is(pgls_x_lambda_ml, 'try-error') &
       !is(pgls_xy_lambda_0, 'try-error') & !is(pgls_xy_lambda_ml, 'try-error')) break  }
  
DataTable_lambda_0[,6,j] <- unlist(lapply(1:length(beta), function(i) pgls_xy_lambda_0[[i]]$model$coef[2,1])) # slopes
DataTable_lambda_ml[,6,j] <- unlist(lapply(1:length(beta), function(i) pgls_xy_lambda_ml[[i]]$model$coef[2,1])) # slopes

DataTable_lambda_0[,5,j] <- unlist(lapply(1:length(beta), function(i) pgls_xy_lambda_0[[i]]$param[2])) # estimated lambdas
DataTable_lambda_ml[,5,j] <- unlist(lapply(1:length(beta), function(i) pgls_xy_lambda_ml[[i]]$param[2])) # estimated lambdas

DataTable_lambda_0[,3,j] <- unlist(lapply(1:length(beta), function(i) summary(pgls_xy_lambda_0[[i]])$coefficients[2,4])) # p values on slope
DataTable_lambda_ml[,3,j] <- unlist(lapply(1:length(beta), function(i) summary(pgls_xy_lambda_ml[[i]])$coefficients[2,4])) # p values on slope

DataTable_lambda_0[,2,j] <- unlist(lapply(1:length(beta), function(i) summary(pgls_xy_lambda_0[[i]])$fstatistic[[1]])) # f stat
DataTable_lambda_ml[,2,j] <- unlist(lapply(1:length(beta), function(i) summary(pgls_xy_lambda_ml[[i]])$fstatistic[[1]])) # f stat

DataTable_lambda_0[,1,j] <- unlist(lapply(1:length(beta), function(i) summary(pgls_xy_lambda_0[[i]])$r.squared)) # rsq
DataTable_lambda_ml[,1,j] <- unlist(lapply(1:length(beta), function(i) summary(pgls_xy_lambda_ml[[i]])$r.squared)) # rsq

progress(j, length(lambdas_expanded), progress.bar = T, init = F, gui = T)
} 

Sys.time()-Start.time # 14 min for tips = 16

# munging data tables 
DataTableMat_lambda_0 <- rbind(DataTable_lambda_0[,,1], DataTable_lambda_0[,,2])
for (i in 3:length(lambdas_expanded)) {
  DataTableMat_lambda_0 <- rbind(DataTableMat_lambda_0, DataTable_lambda_0[,,i])
}
DataTableMat_lambda_0 <- as.data.frame(DataTableMat_lambda_0)
DataTableMat_lambda_0$beta <- rep(beta,(length(lambdas)*nsim))


DataTableMat_lambda_ml <- rbind(DataTable_lambda_ml[,,1], DataTable_lambda_ml[,,2])
for (i in 3:length(lambdas_expanded)) {
  DataTableMat_lambda_ml <- rbind(DataTableMat_lambda_ml, DataTable_lambda_ml[,,i])
}
DataTableMat_lambda_ml <- as.data.frame(DataTableMat_lambda_ml)
DataTableMat_lambda_ml$beta <- rep(beta,(length(lambdas)*nsim))








plot(DataTableMat_lambda_ml$lambda.est ~ DataTableMat_lambda_ml$lambda.input, pch = 19)






# writing files
file_name <- paste("Data_Analyses/Sim_Data/PB_lambda_ml-", n, "-ANOVA.csv", sep = "")
write.csv(DataTableMat_lambda_ml, file_name, row.names = F)

file_name <- paste("Data_Analyses/Sim_Data/PB_lambda_0-", n, "-ANOVA.csv", sep = "")
write.csv(DataTableMat_lambda_0, file_name, row.names = F)


# 32 
# 64 
# 128  
# 256 
# 512 
# 1024 

data <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-512-ANOVA.csv") # run another example where you make lmanba input onto y instead of x
# go reread liams paper on lambda input on both x and y and both
data <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-32.csv") # this is producing same correlation, fix

data_beta <- subset(data, data$beta == 0)

plot(data_beta$lambda.est ~ data_beta$lambda.input, pch = 19)

cor(data_beta$lambda.est, data_beta$lambda.input)


plot(data$slope ~ data$beta, pch = 19)
abline(a=c(0,1)) # parameter estimate in the phylo regression outputs an expected value, with broad range of accuracy. 
# that distribution decreases with increasing phylo size
# look at how the range mean std of slope estimates changes with ml vs 0 lambda inputs
# incorporating lambda does not affect inference in the 
# error rates (x axis = f in pgls lambda 0, y axis = f in pgls labmda ml)

data0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-32.csv")
dataml <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-32.csv")

plot(dataml$'F'~data0$'F')
plot(dataml$'P'~data0$'P')
plot(dataml$slope~data0$slope)
# conclusion holds: interpreting the regression is not appreciably different (in terms of slopes and f values)


plot(DataTableMat_lambda_0$slope ~ DataTableMat_lambda_0$beta)
