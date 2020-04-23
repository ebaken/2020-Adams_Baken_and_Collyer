# Compare variation across runs

treesizes <-  2^(5:10)

lambda <- read.csv("Data_Analyses/Munged_Data/Lambda_est_vars.csv", header=TRUE)
Z.k <- read.csv("Data_Analyses/Munged_Data/Kappa_Z_norm_vars.csv", header=TRUE)
 colnames(lambda) <- colnames(Z.k) <- c("lambda.in","32","64","128","256","512","1024")
kappa_128 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_128.csv")

#Plot for SM
library(latex2exp)
jpeg("Fig2.new.jpeg", res = 120, quality = 100, width = 1000, height = 1000)
par(mfrow=c(2,2))
plot(kappa_128$lambda.est ~ kappa_128$lambda.input, 
     xlab = TeX('$\\lambda_{in}$'), ylab = TeX('$\\lambda_{est}$'), pch = 19)
title("A", adj=0)
text(min(kappa_128$lambda.input),max(kappa_128$lambda.est),
     "n=128",adj= c(0.1,1))
plot(kappa_128$kappa.z ~ kappa_128$lambda.input, 
     xlab = TeX('$\\lambda_{in}$'), ylab = TeX('$Z_K$'),pch = 19)
title("B", adj=0)
text(min(kappa_128$lambda.input),max(kappa_128$kappa.z),
     "n=128",adj= c(0.1,1))
boxplot(lambda[,-1],ylim=c(0,.1),ylab = TeX('$\\sigma^2$'), xlab = "n")
title("C", adj=0)
boxplot(Z.k[,-1],ylim=c(0,.1),ylab = TeX('$\\sigma^2$'), xlab = "n")
title("D", adj=0)
par(mfrow=c(1,1))

dev.off()
