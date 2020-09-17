# Compare Phylogenetic Signal: SA:V vs. BW/SVL 

# Libraries #### 
library(stringr)
library(ape)
library(geomorph) # Morphological analyses
library(geiger) # for pruning trees and datasets to each other (function datatree)

Tree <- read.tree("Bonett_and_Blair_Phylo.tre") # Full Tree
LinMeas <- read.csv("LinMeas.SpeciesSummaries.csv")

row.names(LinMeas) <- LinMeas$Species
TreeDataPruned <- treedata(Tree, LinMeas) 
TreePruned <- TreeDataPruned$phy
LinMeasOrdered <- LinMeas[match(TreePruned$tip.label, LinMeas$Species),]
SVL <- LinMeasOrdered$SVL; names(SVL) <- rownames(LinMeasOrdered)
BWL <- LinMeasOrdered$BWL /SVL; names(BWL) <- rownames(LinMeasOrdered)

DF <- read.csv("SAV.Disparity.csv", header = TRUE)
SAtoV <- DF$SAtoV; names(SAtoV) <-DF$X
DATA <- treedata(TreePruned, SAtoV, warnings= T)

###### Analyses
res.BWL <-physignal(A = BWL, phy = TreePruned)
summary(res.BWL)
plot(res.BWL)
length(BWL)            

res.SAToV<- physignal(A = DATA$data, phy = DATA$phy)
summary(res.SAToV)
plot(res.SAToV)

### Compare Z12
source('newBC.r')
source('CompareZ.r')

compare.Z(res.SAToV,res.BWL)

