# Phylogenetic signal of SA:V
library(stringr)
library(geiger)
library(geomorph)
library(ape)

DF <- read.csv("SAV.Disparity.csv")
row.names(DF) <- DF$X
Tree <- read.tree("BB.PrunedToLMData.tre")
  Tree$tip.label <- str_replace_all(Tree$tip.label, "_", " ")
DATA <- treedata(Tree, DF, warnings= T)
  Tree.SAtoV <- DATA$phy
  Data.SAtoV <- DATA$data[match(Tree.SAtoV$tip.label, row.names(DATA$data)),-1]
SAtoV <- as.numeric(Data.SAtoV[,1])
names(SAtoV) <- row.names(Data.SAtoV)

res.SAToV<- physignal(A = SAtoV, phy = Tree.SAtoV)
summary(res.SAToV)
plot(res.SAToV)

length(Tree.SAtoV$tip.label)
plot(Tree.SAtoV)

