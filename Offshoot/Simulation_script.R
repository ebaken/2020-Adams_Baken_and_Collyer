# Testing consistency of null dist for estimated z scores


# Sim 3 traits on 5 on each tree, look at histograms colored by tree

    # on those trees, simulate bm, trend, ou, and do it again.

library(phytools)
library(geiger)
library(geomorph)


# Plot_A_32: 5 trees, 3 traits on each, all BM #####

ntips <- 2^(5:10)
ntrees <- 5
ntraits <- 3

kappa.obs.null <- array(NA, dim = c(100, ntraits, ntrees))
phy.sim <- lapply(1:ntrees, function(i) pbtree(n = ntips[1]))
cols <-brewer.pal(ntrees, "Dark2")

plot(1, type="n", axes=T, xlab="", ylab="", xlim = c(0,1), ylim = c(0,35))

for (i in 1:ntrees) {
traits.sim <- sim.char(phy = phy.sim[[i]], par = 1, nsim = ntraits)
kappa.obs.null <- lapply(1:ntraits, function(i) physignal(traits.sim[,,i], phy = phy.sim[[i]], print.progress = F, iter = 99)$random.K)
    for(j in 1:ntraits) {
    hist(kappa.obs.null[[j]][-1], col=  alpha(cols[i], .5), breaks = seq(0,1,.01), add = T)
      abline(v = kappa.obs.null[[j]][1], col = cols[i])
    }

}
legend("topright", c("Tree 1" ,"Tree 2", "Tree 3", "Tree 4", "Tree 5"),fill=cols, inset = .03)

# Plot_B_32: 1 tree, 5 traits (1 BM and 4 OU of alphas .25, .5, .75, and 1) ####
ntips <- 2^(5:10)
ntrees <- 1

kappa.obs.null <- array(NA, dim = c(100, ntraits, ntrees))
phy.sim <- pbtree(n = ntips[1])
models <- c("BM", "OU", "OU", "OU", "OU")
alphas <- seq(0,1,.25)

cols <-brewer.pal(5, "Accent")

plot(1, type="n", axes=T, xlab="", ylab="", xlim = c(0,1), ylim = c(0,25))

for (i in 1:5) {

traits.sim <- rTraitCont(phy = phy.sim, model = models[i], alpha = alphas[i])
kappa.obs.null <- physignal(traits.sim, phy = phy.sim, print.progress = F, iter = 99)$random.K

hist(kappa.obs.null[-1], col=  alpha(cols[i], .5), breaks = seq(0,1,.01), add = T)
abline(v = kappa.obs.null[1], col = cols[i], lwd =3)
}

legend("topright", c("BM" ,"OU (alpha = .25)", "OU (alpha = .5)", "OU (alpha = .75)", "OU (alpha = 1)"), fill=cols, inset = .03)

