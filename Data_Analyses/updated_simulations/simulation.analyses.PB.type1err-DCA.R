library(geiger) 
library(phytools)
library(geomorph)

##DCA: adjusted BC to not use observed in calculating the transformation,

treesizes <- 2^(8) # 2^(5:7) done, 2^8 being run by EKB, still need to run 2^(9:10)
lambdas <- seq(0, 1, 0.05)
ntrees <- 2
ndatasets <- 50
physig_iter <- 999

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

#### Run the simulations ####

results <- lapply(1:length(treesizes), function(i){ 
  starttime <- Sys.time()
  filenames <- list.files(path="results-DCA/")
  other_reps <- length(which(startsWith(filenames, paste("pb.", treesizes[i], ".", ntrees, "trees.", ndatasets, "datasets", sep= ""))=="TRUE"))
  repname <- letters[other_reps+1] # adding repname avoids overwriting so we can do it piecemeal 
  file.name <- paste("results-DCA/pb.", treesizes[i], ".", ntrees, "trees.", ndatasets, "datasets.", repname, ".error.csv", sep = "")
  pb <- txtProgressBar(min = 0, max = length(lambdas), initial = 0, style=3)
  
  tree <- pbtree(n=treesizes[i], nsim = ntrees)
  lambdas_output <- lapply(1:length(lambdas), function(j) {
    ind_tree_output <- lapply(1:ntrees, function(jj) {
      tree_scaled <- rescale(tree[[jj]], model = "lambda",  max(lambdas[j], 0.01))
      resid_all <- sim.char(tree_scaled, par = 1, nsim = ndatasets)
      datasets_output <-  lapply(1:ndatasets,  function(k) {
        kappa_est <- physignal(as.matrix(resid_all[,,k]), phy = tree_scaled, iter = physig_iter, print.progress = F)
        random.K.trans <- (box.cox(kappa_est$random.K)$transformed)
#        random.K.trans <- log(kappa_est$random.K)  #LOG-transformed
#        random.K.trans <- (kappa_est$random.K)  #ORIGINAL
        mu <- mean(random.K.trans)
        sd <- sd(random.K.trans)*sqrt((physig_iter - 1)/physig_iter)
        
        list("K" = kappa_est$phy.signal, "mu" = mu, "sd" = sd)
      })
      for (jjj in 1:ndatasets) names(datasets_output)[jjj] <- paste("dataset", jjj, sep = ".")
      datasets_output
    })
    setTxtProgressBar(pb, j)
    for (jjj in 1:ntrees) names(ind_tree_output)[jjj] <- paste("rand.tree", (other_reps*ntrees) + jjj, sep = ".")
    ind_tree_output
  })
  
  endtime <- Sys.time()
  cat(paste(c("\nSimulations for treesize n =", treesizes[i],"completed in: "), sep =" "))
  diff <- endtime - starttime
  cat(diff, attr(diff, "units"))
  cat("\n")
  
  lambda_output_mat <- matrix(unlist(lambdas_output), ncol = 3, byrow = TRUE)
  rand.tree.labels <- (other_reps*ntrees) + 1:ntrees # this only labels trees correctly if using the same ntrees every time. 
  lambda_output_mat <- cbind(rep(treesizes[[i]], each = ndatasets*ntrees*length(lambdas)), rep(lambdas, each = ndatasets*ntrees), rep(rand.tree.labels, each = ndatasets),  lambda_output_mat)
  colnames(lambda_output_mat) <- c("treesize", "lambda_input", "rand.tree", "K", "mu", "sd")
  write.csv(lambda_output_mat, file.name, row.names = F)
  
  
  for(jj in 1:length(lambdas)) {names(lambdas_output)[jj] <- paste("lambda", lambdas[jj], sep = ".")}
  lambdas_output
  
})
names(results) <- treesizes

results_mat <- matrix(unlist(results), ncol = 3, byrow = TRUE)
results_mat <- cbind(rep(treesizes, each = ndatasets*ntrees*length(lambdas)), rep(lambdas, each = ndatasets*ntrees), rep(1:ntrees, each = ndatasets),  results_mat)
colnames(results_mat) <- c("treesize", "lambda_input", "rand.tree", "K", "mu", "sd")
results_df <- as.data.frame(results_mat)

#### Calculate Z12s ####

# if simulations already done, read in simulation data and munge together:
treesizes <- 2^(5:10)
lambdas <- seq(0, 1, 0.05)
ntrees <- 20
ndatasets <- 50
files <- list.files(path = "results-DCA", pattern = "*.error.csv", full.names = T)
all_files <- lapply(files[-c(2,6,10)], read.csv) # excluding the 1000 tree files 
results_in <- do.call(rbind,all_files)
results_df <- as.data.frame(results_in) # edit rand.tree names manually if diff reps had diff ntrees***

# reordering results to match completely
library(dplyr)
results_df_ordered <- arrange(results_df, rand.tree)
results_df_ordered <- arrange(results_df_ordered, lambda_input)
results_df_ordered <- arrange(results_df_ordered, treesize)

identical(results_df_ordered[which(results_df_ordered$treesize == 32),"lambda_input"], results_df_ordered[which(results_df_ordered$treesize == 512),"lambda_input"])

pair.comps <- combn(ndatasets, 2)
z12 <- lapply(1:length(treesizes), function(i){
  results_by_treesize <- results_df_ordered[which(results_df_ordered$treesize == treesizes[i]),]
  lambdas.out <- lapply(1:length(lambdas), function(j) {
    results_by_lambda <- results_by_treesize[which(results_by_treesize$lambda_input == as.character(lambdas[j])),]
    rand.tree.out <- lapply(1:ntrees, function(k) {
      results_by_rand_tree <- results_by_lambda[which(results_by_lambda$rand.tree == k),]
      
      dataset.out <- lapply(1:ncol(pair.comps), function(l){
        abs((results_by_rand_tree[pair.comps[1,l],]$K - 
               results_by_rand_tree[pair.comps[1,l],]$mu) - 
              (results_by_rand_tree[pair.comps[2,l],]$K - 
                 results_by_rand_tree[pair.comps[2,l],]$mu)) / 
          sqrt(results_by_rand_tree[pair.comps[1,l],]$sd^2 + results_by_rand_tree[pair.comps[2,l],]$sd^2)
      })
      for (jj in 1:ncol(pair.comps)) { names(dataset.out)[jj] <- paste(pair.comps[1,jj], pair.comps[2,jj], sep = "-") }
      dataset.out
    })
    names(rand.tree.out) <- 1:ntrees
    rand.tree.out
  })
  names(lambdas.out) <- lambdas
  lambdas.out
}) # 23 min for full dataset
names(z12) <- treesizes
anyNA(z12)

z12_mat <- matrix(unlist(z12), ncol = 1, byrow = F)
z12_mat <- cbind(rep(treesizes, each = length(lambdas)*ntrees*ncol(pair.comps)), 
                 rep(rep(lambdas, each = ntrees*ncol(pair.comps)), length(treesizes)), 
                 rep(rep(1:ntrees, each = ncol(pair.comps)), length(treesizes)*length(lambdas)), 
                 z12_mat)
colnames(z12_mat) <- c("treesize", "lambda_input", "rand.tree", "z12")
write.csv(z12_mat, "results-DCA/pb.alltreesizes.20trees.50datasets.z12s.csv", row.names = F)


#### Calculating Type 1 Error & False Discovery Rates ####
z12_mat <- read.csv("results-DCA/pb.alltreesizes.20trees.50datasets.z12s.csv")


z12_mat <- as.data.frame(z12_mat)

type1error.zk <- lapply(1:length(treesizes), function(i){
  error <- lapply(1:length(lambdas), function(j) {
    z12_subset <- z12_mat[which(z12_mat[,1] == treesizes[i] & z12_mat[,2] == as.character(lambdas[j])),4]
    ps.zk <- pnorm(z12_subset, lower.tail = F) # calculate p value from Z_12 K stat produced above
    length(which(ps.zk < 0.05))/(length(ps.zk))
  })
  names(error) <- lambdas
  error
})
names(type1error.zk) <- treesizes

#### Plot ####

t1err_mat <- matrix(unlist(type1error.zk), ncol = 21, byrow = T)
colnames(t1err_mat) <- lambdas
rownames(t1err_mat) <- treesizes

errormeans <- colMeans(t1err_mat)
cis <- lapply(1:length(lambdas), function(i) qnorm(0.975)*(sd(t1err_mat[,i])/sqrt(6)))

library(gplots)
png("results-DCA/fig.s21.png")
plotCI(x = lambdas, y = errormeans, gap = 0.5,
       xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
       ylab = latex2exp::TeX("Type 1 Error/False Discovery Rates"), 
       ylim = c(0,1), uiw = unlist(cis), pch = 19)
points(x = lambdas, y = errormeans, type = "l", col =  "black")
abline(h = 0.05, col = "red")
dev.off()




# 
png("results-DCA/fig.S21-NotUsed.png")
plot(x = lambdas, y = unlist(type1error.zk$`32`), 
     xlab = latex2exp::TeX("Input Phylogenetic Signal ($\\lambda_{in}$)"),
     ylab = latex2exp::TeX("Type 1 Error/False Discovery Rates"), 
     type = "l", ylim = c(0,1), col = "purple")
points(x = lambdas, y = unlist(type1error.zk$`64`), type = "l", col =  "blue")
points(x = lambdas, y = unlist(type1error.zk$`128`), type = "l", col = "darkgreen")
#points(x = lambdas, y = unlist(type1error.zk$`256`), type = "l", col =  "green")
#points(x = lambdas, y = unlist(type1error.zk$`512`), type = "l", col = "yellow")
#points(x = lambdas, y = unlist(type1error.zk$`1024`), type = "l", col = "orange")
abline(h = 0.05, col = "red")

legend("topleft", as.character(treesizes), fill =  c("purple", "blue", "darkgreen", "green", "yellow", "orange")[1:length(treesizes)])
dev.off()