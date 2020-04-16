# Simulating pgls inference errors from estimated lambdas

# Set up
setwd("/Users/ericabaken/Documents/School/Projects/Lambda_Eval/Manuscript/2020-Adams_and_Baken/")

library(phytools)
library(geomorph)
library(geiger)
#library(svMisc)
library(caper)

 # extract kappa values for comparison

treesizes <-  2^(5:10) # 32,64,128,256,512,1024
nsim <- 50 # number of simulations
lambdas <- seq(0,1,.05) # starting with just 3 lambda values
lambdas_expanded <- rep(lambdas, each = nsim)
fullsimlength <- length(lambdas_expanded)

n <- treesizes[3] # number of tips: adjust each iteration

# New Rep #####
DataTable_lambda <- data.frame(lambda.input = lambdas_expanded, 
                               lambda.est = rep(NA, nsim*length(lambdas)), 
                               kappa = rep(NA, nsim*length(lambdas)),
                               kappa.z = rep(NA, nsim*length(lambdas))) 

pb <- txtProgressBar(0, fullsimlength, style=1)

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
    kappa_est <- try(physignal(X, phy = tree, iter = 999, print.progress = F), silent = F)
    
    if(!is(pgls_est, 'try-error') & !is(kappa_est, 'try-error')) break  }

DataTable_lambda[j,2] <- pgls_est$param[2]
DataTable_lambda[j,3] <- kappa_est$phy.signal 

kappa_sum <- capture.output(summary(kappa_est)) # this capture.output was necessary to silence the print function that cluttered my console. not analytically necessary.
text_z <- kappa_sum[11]
DataTable_lambda[j,4] <- as.numeric(gsub("[^0-9.-]", "", text_z)) # gettin rid of text

setTxtProgressBar(pb, j)
} 
Sys.time()-Start.time

# plot
plot(x = DataTable_lambda$lambda.input, y = DataTable_lambda$lambda.est, pch = 19)
plot(x = DataTable_lambda$lambda.input, y = DataTable_lambda$kappa, pch = 19, ylim = c(0,2.5))

jpeg("Figures/FigS5.jpeg", res = 80, quality = 100, width = 440, height = 990)
par(mfrow = c(3,1))
plot(x = DataTable_lambda$lambda.input, y = DataTable_lambda$lambda.est, pch = 19, xlab = "Input Lambda", ylab = "Estimated Lambda X~1")
plot(x = DataTable_lambda$lambda.input, y = DataTable_lambda$kappa, pch = 19, ylim = c(0,2.5), xlab = "Input Lambda", ylab = "Estimated Kappa")
plot(x = DataTable_lambda$lambda.input, y = DataTable_lambda$kappa.z, pch = 19, xlab = "Input Lambda", ylab = "Kappa Z Scores") # look at this one
dev.off()


# Writing Output Files ####
file_name <- paste("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_", n, ".csv", sep = "")
write.csv(DataTable_lambda, file_name, row.names = F)
