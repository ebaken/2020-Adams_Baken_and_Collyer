#Playing with CI for lambda effect size

library(geiger)
library(phytools)
library(caper)
library(geomorph)

ntaxa <-  2^(5:10) # 32,64,128,256,512,1024
nsim <- 25 # number of simulations
lambdas <- seq(0,1,.05); lambdas[1] <- 0.005
n <- ntaxa[3]

tree<-pbtree(n=n,scale=1) 
tree.scaled <- lapply(1:length(lambdas), function(j) rescale(tree, 
            model = "lambda", lambda = lambdas[j]))
lambda.est <- Z.lambda <- k.val <- Z.k <- matrix(NA, nrow=nsim, ncol=length(lambdas))

X <- lapply(1:length(lambdas), function(j) sim.char(tree.scaled[[j]], par = 1, nsim = nsim) )
for (i in 21:length(lambdas)){
  data_x <- lapply(1:nsim, function(j) comparative.data(tree, data.frame(X.df=X[[i]][,,j], 
                                                                       Species = names(X[[i]][,,j])), Species, vcv = F, vcv.dim = 2))
  pgls_est <- lapply(1:nsim, function(j) pgls(X.df ~ 1, data_x[[j]], lambda='ML'))
  K <- lapply(1:nsim, function(j) physignal(X[[i]][,,j],tree, print.progress = FALSE))
  
  lambda.est[,i] <- unlist(lapply(1:nsim, function(j) pgls_est[[j]]$param[2]))
  CI <- lapply(1:nsim, function(j) pgls_est[[j]]$param.CI$lambda$ci.val)
  
  lambda.std <- unlist(lapply(1:nsim, function(j) {
    if(is.na(CI[[j]])[1]==TRUE){lambda.std <- (CI[[j]][2]-lambda.est[j,i])*sqrt(n)/qnorm(0.975)} 
    else {lambda.std <- (lambda.est[j,i]-CI[[j]][1])*sqrt(n)/qnorm(0.975)}
  }))

  #lambda.std <- (lambda.est[,i]-CI.lo)*sqrt(n)/qnorm(0.975)  #simple version unlooped
  Z.lambda[,i] <- lambda.est[,i]/lambda.std
  
  k.val[,i] <- unlist(lapply(1:nsim, function(j) K[[j]]$phy.signal))
  Z.k[,i] <- unlist(lapply(1:nsim, function(j) K[[j]]$Z))
 
}

#make looped plots
par(mfrow=c(2,2))
matplot(lambdas,t(lambda.est),pch=21,bg="black", col="black",
        xlab="Lambda.in", ylab="lambda.est")
matplot(lambdas,t(k.val),pch=21,bg="black", col="black",
        xlab="Lambda.in", ylab="K.est")
matplot(lambdas,t(Z.lambda),pch=21,bg="black", col="black",
        xlab="Lambda.in", ylab="Z.lambda")
matplot(lambdas,t(Z.k),pch=21,bg="black", col="black" , 
        xlab="Lambda.in", ylab="Z.k")
par(mfrow=c(1,1))

Z.l.mn <- apply(Z.lambda,2,mean)
Z.k.mn <- apply(Z.k.norm,2,mean)
par(mfrow=c(1,2))
plot(lambdas,Z.l.mn)
plot(lambdas,Z.k.mn)
par(mfrow=c(1,1))

apply(Z.lambda,2,var)
apply(Z.k,2,var)

#Coefficient of Variation
CV.l <- sd(apply(Z.lambda,2,var))/   mean(apply(Z.lambda,2,var))*100
CV.k <- sd(apply(Z.k,2,var))/   mean(apply(Z.k,2,var))*100
#Do this for each level of species richness!!! LAMBDA THEN kAPPA.
barplot(c(CV.l,CV.k))

