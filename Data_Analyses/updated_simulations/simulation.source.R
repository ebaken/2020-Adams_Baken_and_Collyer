
# Source file for all needed simulation functions


effects.simulations <- function(treesize , tree.type = c("bal", "pect", "pb"),
                                nsim = 50, lambdas, poly.per = 0,
                                beta = 1, mod.type = c("mean", "reg", "anova")) { 
  # poly.per is franction (percent) of nodes that have polytomies
  # beta should be 1, unless varied
  
  require(phytools)
  require(geomorph)
  require(geiger)
  require(caper)
  
  starttime <- Sys.time()
  pb <- txtProgressBar(min = 0, max = nsim, initial = 0, style=3)
  
  L <- as.list(environment())
  file.name <- paste(paste(L[c("tree.type", "treesize", "nsim")], collapse = "."), "sim",
        L$poly.per, "poly", L$mod.type, "type.csv", sep = ".")
  
  tree.type <- match.arg(tree.type)
  mod.type <- match.arg(mod.type)
  
  n <- treesize

  x <- if(mod.type == "mean" || mod.type == "reg") matrix(1, n) else 
    sample(matrix(rep(c(0, 1), each  = n/2)))
  
  stat.args <- list(tree = 0, x = x, beta = beta[1], lambda = lambdas[1],
                    type = mod.type)
  tree.args <- list(n = n, type = tree.type, poly.per = poly.per, lambda = lambdas[1])
 
  Result <- lapply(1:nsim, function(j) {
    
    result <- lapply(1:length(lambdas), function(jj){
      lambda <- lambdas[jj]
      stat.args$lambda <- tree.args$lambda <- lambda
      stat.args$tree  <- do.call(make.tree, tree.args)
    
      stats <- sapply(1:length(beta), function(jjj){
        TE <- TRUE
        stat.args$beta <- beta[jjj]
        while(TE) {
          res <- do.call(get.stats, stat.args)
          TE <- inherits(res, "try-error")
        }
        res <- c(lambda.in = stat.args$lambda, beta.in = stat.args$beta, res)
        res
      })
      t(stats)
    })
    
    result <- do.call(rbind, result)
    setTxtProgressBar(pb, j)
    result
    
  })
  
  Result <- do.call(rbind, Result)
  
  close(pb)
  endtime <- Sys.time()
  cat("\nSimulations completed in: ")
  diff <- endtime - starttime
  cat(diff, attr(diff, "units"))
  cat("\n")

  list(Result = as.data.frame(Result), file.name = file.name)
  
}

effects.simulations.PP <- function(treesize , tree.type = c("bal", "pect", "pb"),
                                nsim = 50, lambdas, poly.per = 0,
                                beta = 1, mod.type = c("mean", "reg", "anova")) { 
  require(parallel)
  require(phytools)
  require(geomorph)
  require(geiger)
  require(caper)
  
  starttime <- Sys.time()
  
  L <- as.list(environment())
  file.name <- paste(paste(L[c("tree.type", "treesize", "nsim")], collapse = "."), "sim",
                     L$poly.per, "poly", L$mod.type, "type.csv", sep = ".")
  
  tree.type <- match.arg(tree.type)
  mod.type <- match.arg(mod.type)
  
  n <- treesize
  
  x <- if(mod.type == "mean" || mod.type == "reg") matrix(1, n) else matrix(rep(c(0, 1), each  = n/2)) 
  
  stat.args <- list(tree = 0, x = x, beta = beta[1], lambda = lambdas[1],
                    type = mod.type)
  tree.args <- list(n = n, type = tree.type, poly.per = poly.per, lambda = lambdas[1])
  cores <- detectCores() - 1
  
  Result <- mclapply(1:nsim, function(j) {
    
    result <- lapply(1:length(lambdas), function(jj){
      lambda <- lambdas[jj]
      stat.args$lambda <- tree.args$lambda <- lambda
      stat.args$tree  <- do.call(make.tree, tree.args)
      
      stats <- sapply(1:length(beta), function(jjj){
        TE <- TRUE
        stat.args$beta <- beta[jjj]
        while(TE) {
          res <- do.call(get.stats, stat.args)
          TE <- inherits(res, "try-error")
        }
        res <- c(lambda.in = stat.args$lambda, beta.in = stat.args$beta, res)
        res
      })
      t(stats)
    })
    
    result <- do.call(rbind, result)
    setTxtProgressBar(pb, j)
    result
    
  }, mc.cores = cores)
  
  Result <- do.call(rbind, Result)
  
  endtime <- Sys.time()
  cat("\nSimulations completed in: ")
  diff <- endtime - starttime
  cat(diff, attr(diff, "units"))
  cat("\n")
  
  list(Result = as.data.frame(Result), file.name = file.name)
  
}


make.tree <- function(n, type, poly.per, lambda) {
  tree <- if(type == "bal") compute.brlen(stree(n = n, type = "balanced")) else
      if(type == "pect") compute.brlen(stree(n = n, type = "left")) else
        pbtree(n = n, scale = 1)
  
  gg <- function(tree, poly.per){ # From Molina-Venegas and Rodríguez (2017)
    k <- round((Ntip(tree) -2) * poly.per) 
    if(k > 0) {
      for(i in 1:k){
        ii <- sample(2:tree$Nnode,1) + length(tree$tip.label)
        ll <- tree$edge.length[which(tree$edge[,2] == ii)]
        tree$edge.length[which(tree$edge[,1] == ii)] <-
          tree$edge.length[which(tree$edge[,1] == ii)] + ll
        tree$edge.length[which(tree$edge[,2] == ii)] <- 0
        tree<-di2multi(tree)
      }
    }
    tree
  }
  
  tree <- gg(tree, poly.per)
  tree_scaled <- rescale(tree, model = "lambda",  max(lambda, 0.01))
  
  c(tree=tree, tree_scaled = tree_scaled)
  
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

get.stats <- function(tree, x, beta, lambda, type) {
  tree_scaled <- tree$tree_scaled
  tree <- tree$tree
  Sig <- matrix(0, 2, 2)
  diag(Sig) <- 1
  resid <- if(type == "reg") sim.char(tree_scaled, par = Sig)[,,1] else 
    sim.char(tree_scaled, par = 1)[,,1]
  resid <- as.matrix(resid)
  if(type == "reg") x <- resid[,2]
  mu <- x * beta
  y <- if(type != "mean") mu + resid[,1] else resid[,1]
  y <- as.matrix(y)
  x <- as.matrix(x)
  resid <- as.matrix(resid)
  rownames(y) <- rownames(x) <- rownames(resid)
  
  df <- data.frame(x = x, y = y, Species = rownames(y))
  df <- comparative.data(tree, df, "Species", vcv = F, vcv.dim = 2, force.root = TRUE)
  
  pgls_est  <- if(type == "mean") try(pgls(y ~ 1, data = df, lambda='ML'), silent = T) else
    try(pgls(y ~ x, data = df, lambda='ML'), silent = T)
  kappa_est <- physignal(cbind(y,0), phy = tree, iter = 999, print.progress = F)
  if(!inherits(pgls_est, "try-error")) {
    summ <- summary(pgls_est)
    Pval <- unlist(summ$coefficients[,4])
    names(Pval) <- paste("Pval", names(Pval), sep = ".")
    res <- c(unlist(pgls_est$model), Pval, unlist(pgls_est$param.CI)[c(1,6,7)],
             lambda.z = lambda.Z(pgls_est), 
             unlist(kappa_est[c("phy.signal", "Z")]))
    names(res)[length(res)] <- "kappa.z"
  } else res <- pgls_est
  
  res
  
}


lambda.plot <- function(R, beta = 1, ...) {
  df <- R$Result
  df <- df[df$beta.in == beta,]
  plot(df$lambda.in, df$lambda.opt.lambda, ...)
  
}

kappa.plot <- function(R, beta = 1, ...) {
  df <- R$Result
  df <- df[df$beta.in == beta,]
  plot(df$lambda.in, df$phy.signal, ...)
  
}

lambda.z.plot <- function(R, mean.col = 2, mean.lwd = 2, beta = 1, ...) {
  df <- R$Result
  df <- df[df$beta.in == beta,]
  plot(df$lambda.in, df$lambda.z, ...)
  means <- by(df$lambda.z, df$lambda.in, mean)
  points(unique(df$lambda.in), means, type = "l", col = mean.col, lwd = mean.lwd)
  
}

kappa.z.plot <- function(R, mean.col = 2, mean.lwd = 2, beta = 1, ...) {
  df <- R$Result
  df <- df[df$beta.in == beta,]
  plot(df$lambda.in, df$kappa.z, ...)
  means <- by(df$kappa.z, df$lambda.in, mean)
  points(unique(df$lambda.in), means, type = "l", col = mean.col, lwd = mean.lwd)
  
}

regression.precision.plot <- function(R, mean.col = 2, mean.lwd = 2,  ...) {
  df <- R$Result
  plot(df$beta.in, df$coef.x, ...)
  means <- by(df$coef.x, df$beta.in, mean)
  points(unique(df$beta.in), means, type = "l", col = mean.col, lwd = mean.lwd)
  abline(0, 1)
}
  
