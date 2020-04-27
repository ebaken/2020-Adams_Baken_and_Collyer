#' Comparisons of Empirically-Derived Effect Sizes from Multiple Analyses
#'
#' Function performs an analysis to compare the empirically-generated effect sizes of two or more analyses
#'
#' The function statistically compares the effect sizes derived from empirical sampling distributions for 
#' two or more statistical analyses. The approach is generalized from effect sizes derived from residual 
#' randomization permutation procedures (RRPP: Collyer et al. 2015: see also Adams and Collyer 2016; 2019). 
#' The analysis calculates effect sizes as standard deviates, z, from empirically-derived sampling distributions.
#' These provide a measure of the strength of signal for a particular analysis (e.g., the strength of associaton 
#' from a partial least squares analysis, the strength of the effect of a factor in a linear model, etc.).
#' Effect sizes are then statistically compared using a two-sample z-test, based on the pooled 
#' standard error from the sampling distributions of the analyses. The method follows suggestions in 
#' Adams and Collyer (2019), which generalized the use of empirically-derived sampling distributions from 
#' permutation as a means to compare the strength of patterns in multivariate datasets. 
#' 
#' To use this function, simply perform a series of analyses in geomorph or RRPP. Any number of objects can be used as 
#' input. The function will then determine the class of each object: "CR", "evolrate", "lm.rrpp", 
#' "physignal" or "pls", and will extract the appropriate components. It will then estimate effect sizes for each object
#' and perform pairwise comparisons of the strength of signal of those objects based on a two-sample z-score (NOTE: 
#' If one linear model is used as input the function will compare effect sizes across factors within that linear
#' model). Finally, one may perform the comparison as either a one-tailed or a two-tailed (default) test.
#' 
#' @param ... saved analyses of any RRPP or geomorph class
#' @param two.tailed A logical value to indicate whether a two-tailed test (typical and default) should be performed.
#' @keywords analysis
#' @export
#' @author Dean Adams and Michael Collyer
#' @return An object of class compare.Z, returns a list of the following
#' \item{sample.z}{A vector of effect sizes for each sample.}
#' \item{sample.r.sd}{A vector of standard deviations for each sampling distribution.}
#' \item{pairwise.z}{A matrix of pairwise, two-sample z scores between all pairs of effect sizes.}
#' \item{pairwise.p}{A matrix of corresponding P-values.}
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357-365.
#' @references Adams, D.C. and M.L. Collyer. 2016.  On the comparison of the strength of morphological integration across morphometric 
#' datasets. Evolution. 70:2623-2631.
#' @references Adams, D.C. and M.L. Collyer. 2019.  Comparing the strength of modular signal, and evaluating alternative modular hypotheses,
#' using covariance ratio effect sizes for morphometric data. *Evolution*. (Accepted).

compare.Z <- function(...,two.tailed = TRUE){
  dots <- list(...)
  tails <- if(two.tailed) 2 else 1
  if(length(dots) == 1) n <- length(dots[[1]]) else n <- length(dots)
  if(n == 1) stop("At least two objects are needed")
  if(length(dots) == 1) {
    list.names <- names(dots[[1]]) 
    dots <- lapply(1:n, function(j) dots[[1]][[j]])
    names(dots) <- list.names
  } else list.names <- names(dots)
  if(length(dots) < 2) stop("At least two objects are needed")
  if(is.null(list.names)) list.names <- paste("Analysis", 1:n, sep = ".")
  names(dots) <- list.names
  random.val <- lapply(1:length(dots), function(j) { #Check class; log random values as needed
    if (is(dots[[j]],"lm.rrpp")  && is.null(dots[[j]]$ANOVA$SS))  { rep(0,dim(dots[[j]]$Y)[1])}  #For intercept model
    else if (is(dots[[j]],"lm.rrpp") && !is.null(dots[[j]]$ANOVA$SS))  {log(dots[[j]]$ANOVA$SS)}  #NOTE: SS  only for now
    else if (is(dots[[j]],"CR")) {t(as.matrix(dots[[j]]$random.CR))} 
    else if (is(dots[[j]],"pls")) {t(as.matrix(dots[[j]]$random.r))} 
    else if (is(dots[[j]],"physignal")) {t(as.matrix(log(dots[[j]]$random.K)))} 
    else if (is(dots[[j]],"evolrate")) {t(as.matrix(log(dots[[j]]$random.sigma)))}
  })
  sdn <- function(x) sqrt(sum((x-mean(x))^2)/length(x))
  getdrs <- function(fit){  
    if(is.null(fit)) 0 else
      fit[1] - mean(fit[-1])
  }
  getsds <- function(fit){
    if(is.null(fit)) 0 else
      sdn(fit[-1])
  }
  getZ <- function(fit) {
    if(is.null(fit)) 0 else
      RRPP:::effect.size(as.vector(fit))
  }
  #added unlist b/c could have multiple factors in lm.rrpp
  list.drs <- unlist(sapply(1:length(random.val), function(j) { apply(random.val[[j]],1,getdrs) }))
  list.sds <- unlist(sapply(1:length(random.val), function(j) { apply(random.val[[j]],1,getsds) }))
  list.zs <- unlist(sapply(1:length(random.val), function(j) { apply(random.val[[j]],1,getZ) }))

  k <- length(list.drs)
  k.combn <- combn(k,2)  
  z12 <- sapply(1:ncol(k.combn), function(j) {
    a <- k.combn[1, j]
    b <- k.combn[2, j]
    r1 <- list.drs[a]
    r2 <- list.drs[b]
    sd1 <- list.sds[a]
    sd2 <- list.sds[b]
    (r1 - r2)/sqrt(sd1^2 + sd2^2)
  })

  z12.p <- sapply(1:length(z12), function(j) pnorm(abs(z12[[j]]), 
                                                   lower.tail = FALSE) * tails)

  #NOTE: names are not correct when there are factors in lm.rrpp  
  d <- rep(0,k); names(d) <- list.names
  D <-dist(d)
  z12.pw <- p12.pw <- D
  for(i in 1:length(z12)) z12.pw[i] <-z12[i]
  for(i in 1:length(z12)) p12.pw[i] <-z12.p[i]
  names(list.zs) <- names(list.sds) <-list.names
  pairwise.z <- as.matrix(z12.pw)
  pairwise.P <- as.matrix(p12.pw)
  diag(pairwise.P) <- 1
  
  out <- list(sample.z = list.zs,
              sample.r.sd.pk.stand = list.sds,
              pairwise.z = abs(pairwise.z),
              pairwise.P = pairwise.P)
  class(out) <- "compare.Z"
  out
}
