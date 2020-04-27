# Phylogenetic Signal of Foot Shape and SVL

#TRY A RATIO

# Libraries #### 
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
SE <- LinMeasOrdered$SE  /SVL; names(SE) <- rownames(LinMeasOrdered)
HL <- LinMeasOrdered$HL  /SVL; names(HL) <- rownames(LinMeasOrdered)
FLL <- LinMeasOrdered$FLL  /SVL; names(FLL) <- rownames(LinMeasOrdered)
HLL <- LinMeasOrdered$HLL  /SVL; names(HLL) <- rownames(LinMeasOrdered)
TL <- LinMeasOrdered$TL  /SVL; names(TL) <- rownames(LinMeasOrdered)


physignal(A = SVL, phy = TreePruned)
physignal(A = BWL, phy = TreePruned)
physignal(A = SE, phy = TreePruned)
physignal(A = HL, phy = TreePruned)
physignal(A = FLL, phy = TreePruned)
physignal(A = HLL, phy = TreePruned)
physignal(A = TL, phy = TreePruned)

res.BWL <-physignal(A = BWL, phy = TreePruned)
summary(res.BWL)
plot(res.BWL)
length(BWL)            

# Foot Shape Analyses ####
Tree <- read.tree("Bonett_and_Blair_Phylo.tre") # Full Tree
Microhabitats <- read.csv("Microhabitats.csv")  # Full Tree Microhabitat Data
FootLMs <- readland.tps("FootShape.SpeciesSummaries.tps", specID = "ID")
dimnames(FootLMs)[[3]] <- stringr::str_replace_all(dimnames(FootLMs)[[3]], " ", "_")
MicrohabitatsPruned <- Microhabitats[match(dimnames(FootLMs)[[3]],Microhabitats$Species),] 
row.names(MicrohabitatsPruned) <- as.character(MicrohabitatsPruned$Species)

PrunedPhyData <- treedata(Tree, MicrohabitatsPruned)
TreePruned.Foot <- PrunedPhyData$phy

res.Foot <- physignal(A = FootLMs, phy = TreePruned.Foot)
summary(res.Foot)
plot(res.Foot)
