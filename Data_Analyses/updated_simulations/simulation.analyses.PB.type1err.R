
treesizes <- 2^(5:7)
lambdas <- seq(0, 1, 0.05)
ntrees <- 10
ndatasets <- 50
physig_iter <- 999

##### functions from develop branch ######
center <- function(x){
  if(is.vector(x)) x - mean(x) else {
    x <- as.matrix(x)
    dims <- dim(x)
    fast.center(x, dims[1], dims[2])
  }
}
box.cox <- function(y, eps = 0.02) {
  if(any(y <= 0)) y = y - min(y) + 0.0001
  n <- length(y)
  yy <- y/exp(mean(log(y)))
  logy <- log(yy)
  lambda <- seq(-2.4, 2.4, 0.4)
  m <- length(lambda)
  
  loglik <- sapply(1:m, function(j){ # same as MASS::boxcox loglik f
    la <- lambda[j]
    yt <- if(abs(la) > eps) yt <- (yy^la - 1)/la else
      logy * (1 + (la * logy)/2 * (1 + (la * logy)/3 * (1 + (la * logy)/4)))
    
    -n/2 * log(sum(center(yt)^2))
  })
  
  lambda.opt <- lambda[which.max(loglik)][[1]]
  
  if(abs(lambda.opt) == 2.4) {
    lambda <- seq(-3.2, 3.2, 0.2)
    m <- length(lambda)
    
    loglik <- sapply(1:m, function(j){ # same as MASS::boxcox loglik f
      la <- lambda[j]
      yt <- if(abs(la) > eps) yt <- (yy^la - 1)/la else
        logy * (1 + (la * logy)/2 * (1 + (la * logy)/3 * (1 + (la * logy)/4)))
      
      -n/2 * log(sum(center(yt)^2))
    })
    lambda.opt <- lambda[which.max(loglik)][[1]]
  }
  
  sp <- spline(lambda, loglik, n = 100)
  lambda.opt <- sp$x[which.max(sp$y)]
  if(abs(lambda.opt) < eps) lambda.opt <- 0
  res <- if(lambda.opt == 0) log(y) else (y^lambda.opt - 1)/lambda.opt
  list(opt.lambda = lambda.opt, transformed = res, lambda = sp$x, loglik = sp$y)
}

#### Runs ####

results <- lapply(1:length(treesizes), function(i){ 
  starttime <- Sys.time()
  file.name <- paste("Data_Analyses/updated_simulations/results/pb", treesizes[i], "error.csv", sep = ".")
  pb <- txtProgressBar(min = 0, max = length(lambdas), initial = 0, style=3)
  
  tree <- pbtree(n=treesizes[i], nsim = ntrees)
  lambdas_output <- lapply(1:length(lambdas), function(j) {
    tree_scaled <- rescale(tree[[i]], model = "lambda",  max(lambdas[j], 0.01))
    resid_all <- sim.char(tree_scaled, par = 1, nsim = ndatasets)
    datasets_output <-  lapply(1:ndatasets,  function(k) {
      kappa_est <- physignal(as.matrix(resid_all[,,k]), phy = tree_scaled, iter = physig_iter, print.progress = F)
      random.K.trans <- center(box.cox(kappa_est$random.K)$transformed)
      mu <- mean(random.K.trans)
      sd <- sd(random.K.trans)*sqrt((physig_iter - 1)/physig_iter)
      
      list("K" = kappa_est$phy.signal, "mu" = mu, "sd" = sd)
    })
    setTxtProgressBar(pb, j)
    for (jj in 1:ndatasets) names(datasets_output)[jj] <- paste("dataset", jj, sep = ".")
    datasets_output
    
  })
  
  endtime <- Sys.time()
  cat(paste(c("\nSimulations for treesize n =", treesizes[i],"completed in: "), sep =" "))
  diff <- endtime - starttime
  cat(diff, attr(diff, "units"))
  cat("\n")
  
  lambda_output_mat <- matrix(unlist(lambdas_output), ncol = 3, byrow = TRUE)
  lambda_output_mat <- cbind(rep(lambdas, each = ndatasets), lambda_output_mat)
  head(lambda_output_mat)
  colnames(lambda_output_mat) <- c("lambda_input", "K", "mu", "sd")
  write.csv(lambda_output_mat, file.name)
  
  for(jj in 1:length(lambdas)) {names(lambdas_output)[jj] <- paste("lambda", lambdas[jj], sep = ".")}
  lambdas_output
  
})
names(results) <- treesizes



pair.comps <- combn(ndatasets, 2)

z12 <- lapply(1:length(treesizes), function(i){
  lambdas.out <- lapply(1:length(lambdas), function(j) {
    pair.out <- lapply(1:ncol(pair.comps), function(k){
      abs((results[[i]][[j]][[pair.comps[1,k]]]$K - 
             results[[i]][[j]][[pair.comps[1,k]]]$mu) - 
            (results[[i]][[j]][[pair.comps[2,k]]]$K - 
               results[[i]][[j]][[pair.comps[2,k]]]$mu)) / 
        sqrt(results[[i]][[j]][[pair.comps[1,k]]]$sd^2 + results[[i]][[j]][[pair.comps[2,k]]]$sd^2)
    })
    for (jj in 1:ncol(pair.comps)) { names(pair.out)[jj] <- paste(pair.comps[1,jj], pair.comps[2,jj], sep = "-") }
    pair.out
  })
  names(lambdas.out) <- lambdas
  lambdas.out
})
names(z12) <- treesizes

type1error.zk <- lapply(1:length(treesizes), function(i){
  error <- lapply(1:length(lambdas), function(j) {
    ps.zk <- pnorm(unlist(z12[[i]][[j]]), lower.tail = F) # calculate p value from Z_12 K stat produced above
    length(which(ps.zk < 0.05))/(length(ps.zk))
  })
  names(error) <- lambdas
  error
})
names(type1error.zk) <- treesizes
type1error.zk



plot(x = lambdas, y = unlist(type1error.zk$`32`), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("Type 1 Error/False Discovery Rates"), 
     type = "l", ylim = c(0,1), col = "purple")
points(x = lambdas, y = unlist(type1error.zk$`64`), type = "l", col =  "blue")
points(x = lambdas, y = unlist(type1error.zk$`128`), type = "l", col = "darkgreen")
abline(h = 0.05, col = "red")

legend("topleft", as.character(treesizes), fill =  c("purple", "blue", "darkgreen", "green", "yellow", "orange")[1:length(treesizes)])
