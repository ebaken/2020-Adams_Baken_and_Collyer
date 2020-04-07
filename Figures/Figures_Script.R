# Figures
treesizes <-  2^(5:10)

ANOVA.data <- read.csv("Data_Analyses/Munged_Data/ANOVA.fulldataset.csv")
ANOVA.data <- as.data.frame(ANOVA.data)

Regression.data <- read.csv("Data_Analyses/Munged_Data/Regression.fulldataset.csv")
Regression.data <- as.data.frame(Regression.data)

## Figure 1 ######

library(cowplot)
library(ggplot2)

PlotList <- lapply(1:6, function(j){
   reduced.data <- Regression.data[which(Regression.data$tree.size == treesizes[[j]]),]
   ggplot(reduced.data, aes(x = lambda.input, y = lambda.est)) +
      geom_point() + theme(legend.position= "none", panel.background = element_rect("transparent")) + 
      xlab("Input Lambda") + ylab("Estimated Lambda") + ggtitle(treesizes[[j]])
})

jpeg("Figures/Fig1.jpeg", res = 120, quality = 100, width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()

jpeg("Manuscript/Fig1.jpeg", res = 120, quality = 100, width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()

## Figure 2 ######

ANOVA.data$beta <- as.factor(ANOVA.data$beta)

ANOVA.data.ml <- ANOVA.data[which(ANOVA.data$method == "ml"),]

PlotList <- lapply(1:6, function(j){
   reduced.data <- ANOVA.data.ml[which(ANOVA.data.ml$tree.size == treesizes[[j]]),]
   ggplot(reduced.data, aes(x = beta, y = slope)) + #geom_abline(intercept = 0, slope = .2) +
      geom_boxplot() + theme(legend.position= "none", panel.background = element_rect("transparent")) + 
      xlab("Input Slope") + ylab("Estimated Slope") + ggtitle(treesizes[[j]]) + expand_limits(y=c(-3.5,4))
})

jpeg("Figures/Fig2.jpeg", res = 120, quality = 100, width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()

jpeg("Manuscript/Fig2.jpeg", res = 120, quality = 100, width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()

# other plot options for Figure 2 

ANOVA.data.0 <- ANOVA.data[which(ANOVA.data$method == "0"),]
 
ANOVA.data.modified <- cbind(ANOVA.data.ml$tree.size, ANOVA.data.ml$lambda.input, ANOVA.data.ml$'F', ANOVA.data.0$'F')
ANOVA.data.modified <- as.data.frame(ANOVA.data.modified)
colnames(ANOVA.data.modified) <- c("tree.size", "lambda.input", "F.ml", "F.0")


PlotList <- lapply(1:6, function(j){
   reduced.data <- ANOVA.data.modified[which(ANOVA.data.modified$tree.size == treesizes[[j]]),]
   ggplot(reduced.data, aes(x = F.ml, y = F.0)) +
      geom_point() + theme(legend.position= "none", panel.background = element_rect("transparent")) + 
      xlab("F from PGLS (est lambda)") + ylab("F from OLS") + ggtitle(treesizes[[j]]) 
})

plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])


## Figure 3 #####
library(readxl)
review <- read_excel("Literature_Review/1Summary.xlsx", sheet = "Since 2019 Lambda Vals", col_types = "numeric")

published.lambdas <- review$`Estimated Lambdas Simple`
published.lambdas.pruned <- published.lambdas[-which(published.lambdas > 1 | published.lambdas < 0)]
published.lambdas.df <- as.data.frame(published.lambdas.pruned)
colnames(published.lambdas.df) <- "Published_Lambdas"

jpeg("Figures/Fig3.jpeg", res = 120, quality = 100, width = 1000, height = 660)
    ggplot(published.lambdas.df, aes(Published_Lambdas)) + 
      geom_histogram(binwidth = 0.02) + theme_bw() +
      xlab("Published Lambda Values") + ylab("Frequency")
dev.off()

jpeg("Manuscript/Fig3.jpeg", res = 120, quality = 100, width = 1000, height = 660)
    ggplot(published.lambdas.df, aes(Published_Lambdas)) + 
      geom_histogram(binwidth = 0.02) + theme_bw() +
      xlab("Published Lambda Values") + ylab("Frequency")
dev.off()


## Figure S1 #####

PlotList <- lapply(1:6, function(j){
   reduced.data <- ANOVA.data[which(ANOVA.data$tree.size == treesizes[[j]]),]
   ggplot(reduced.data, aes(x = lambda.input, y = lambda.est)) +
      geom_point() + theme(legend.position= "none", panel.background = element_rect("transparent")) + 
      xlab("Input Lambda") + ylab("Estimated Lambda") + ggtitle(treesizes[[j]])
})

jpeg("Figures/FigS1.jpeg", res = 120, quality = 100, width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()


# Rmarkdown PDF rendering ####

rmarkdown::render("Manuscript/Manuscript.Rmd", output_dir = "Manuscript/")
