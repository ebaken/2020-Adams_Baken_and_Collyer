#### functions from develop branch ######
center <- function(x){
  if(is.vector(x)) x - mean(x) else {
    x <- as.matrix(x)
    dims <- dim(x)
    fast.center(x, dims[1], dims[2])
  }
}
#Adjust box.cox to use distribution minus observed, then transform observed
box.cox <- function(y, eps = 0.02) {
  y.obs <- y[1];   y <- y[-1] #pull out observed
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
  y <- c(y.obs,y) #put back in here
  res <- if(lambda.opt == 0) log(y) else (y^lambda.opt - 1)/lambda.opt
  list(opt.lambda = lambda.opt, transformed = res, lambda = sp$x, loglik = sp$y)
}