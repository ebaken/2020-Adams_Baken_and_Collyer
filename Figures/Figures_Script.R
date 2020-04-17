# Figures
treesizes <-  2^(5:10)

ANOVA.data <- read.csv("Data_Analyses/Munged_Data/ANOVA.fulldataset.csv")
ANOVA.data <- as.data.frame(ANOVA.data)

Regression.data <- read.csv("Data_Analyses/Munged_Data/Regression.fulldataset.csv")
Regression.data <- as.data.frame(Regression.data)

## Figure 1 ######

library(cowplot)
library(ggplot2)
Regression.data.ml <- Regression.data[(Regression.data$method == "ml"),]

Regression.data.ml

data_list <- list(Regression.data.ml[which(Regression.data.ml$tree.size == treesizes[1]),],
                  Regression.data.ml[which(Regression.data.ml$tree.size == treesizes[2]),],
                  Regression.data.ml[which(Regression.data.ml$tree.size == treesizes[3]),],
                  Regression.data.ml[which(Regression.data.ml$tree.size == treesizes[4]),],
                  Regression.data.ml[which(Regression.data.ml$tree.size == treesizes[5]),],
                  Regression.data.ml[which(Regression.data.ml$tree.size == treesizes[6]),])

fit_list <- lapply(1:6, function(i) lm(data_list[[i]]$lambda.est.x~data_list[[i]]$lambda.input))
slopes <- unlist(lapply(1:6, function(i) coef(fit_list[[i]])[2]))
intercepts <- unlist(lapply(1:6, function(i) coef(fit_list[[i]])[1]))

PlotList <- lapply(1:6, function(j){
   reduced.data <- Regression.data.ml[which(Regression.data.ml$tree.size == treesizes[[j]]),]
   ggplot(reduced.data, aes(x = lambda.input, y = lambda.est.x)) +
      geom_point() + theme(legend.position= "none", panel.background = element_rect("transparent")) + 
      xlab("Input Lambda") + ylab("Estimated Lambda") + ggtitle(treesizes[[j]]) + 
     geom_abline(slope = slopes[[j]], intercept = intercepts[[j]], color = "red")
})

jpeg("Figures/Fig1.jpeg", res = 120, quality = 100, width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()

jpeg("Manuscript/Fig1.jpeg", res = 120, quality = 100, width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()


## Figure 2 ####
kappa_128 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_128.csv")

kappa.norm <- (kappa_128$kappa-min(kappa_128$kappa))/(range(kappa_128$kappa)[2] - range(kappa_128$kappa)[1])

kappa.z.norm <- (kappa_128$kappa.z-min(kappa_128$kappa.z))/(range(kappa_128$kappa.z)[2] - range(kappa_128$kappa.z)[1])


jpeg("Figures/Fig2.jpeg", res = 120, quality = 100, width = 1000, height = 1000)
layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE))
par(mar = c(4,4,2,2))
plot(kappa_128$lambda.est ~ kappa_128$lambda.input, pch = 19)
plot(kappa_128$kappa ~ kappa_128$lambda.input, pch = 19)
plot(kappa.norm ~ kappa_128$lambda.input, pch = 19)
plot(kappa_128$kappa.z ~ kappa_128$lambda.input, pch = 19)
plot(kappa.z.norm ~ kappa_128$lambda.input, pch = 19)
dev.off()



## Figure 3 #####
library(readxl)
review <- read_excel("Literature_Review/1Summary.xlsx", sheet = "Since 2019 Lambda Vals", col_types = "numeric")
colnames(review)
published.lambdas <- review[,c(3,4)]
published.lambdas.pruned <- published.lambdas[-which(published.lambdas > 1 | published.lambdas < 0),]
published.lambdas.df <- as.data.frame(published.lambdas.pruned)
colnames(published.lambdas.df) <- c("Published_Lambdas", "Taxa")

jpeg("Figures/Fig3_bw.jpeg", res = 120, quality = 100, width = 1000, height = 660)
    ggplot(published.lambdas.df, aes(Published_Lambdas)) + 
      geom_histogram(binwidth = 0.02) + theme_bw() +
      xlab("Published Lambda Values") + ylab("Frequency")
dev.off()


# Combining published data to simulated data

published.lambdas.df$Taxa_Bins <- as.character(published.lambdas.df$Taxa)
published.lambdas.df$Taxa_Bins[which(published.lambdas.df$Taxa<=tree.sizes[[1]])] <- "Up To 32"
published.lambdas.df$Taxa_Bins[which(published.lambdas.df$Taxa>tree.sizes[[1]] & published.lambdas.df$Taxa<=tree.sizes[[2]])] <- "33 to 64"
published.lambdas.df$Taxa_Bins[which(published.lambdas.df$Taxa>tree.sizes[[2]] & published.lambdas.df$Taxa<=tree.sizes[[3]])] <- "65 to 128"
published.lambdas.df$Taxa_Bins[which(published.lambdas.df$Taxa>tree.sizes[[3]] & published.lambdas.df$Taxa<=tree.sizes[[4]])] <- "129 to 256"
published.lambdas.df$Taxa_Bins[which(published.lambdas.df$Taxa>tree.sizes[[4]] & published.lambdas.df$Taxa<=tree.sizes[[5]])] <- "257 to 512"
published.lambdas.df$Taxa_Bins[which(published.lambdas.df$Taxa>tree.sizes[[5]] & published.lambdas.df$Taxa<=tree.sizes[[6]])] <- "513 to 1024"
published.lambdas.df$Taxa_Bins[which(published.lambdas.df$Taxa>tree.sizes[[6]])] <- "Over 1024"

published.lambdas.df$Taxa_Bins <- as.factor(published.lambdas.df$Taxa_Bins)
published.lambdas.df$Taxa_Bins <- relevel(published.lambdas.df$Taxa_Bins, "513 to 1024")
published.lambdas.df$Taxa_Bins <- relevel(published.lambdas.df$Taxa_Bins, "257 to 512")
published.lambdas.df$Taxa_Bins <- relevel(published.lambdas.df$Taxa_Bins, "129 to 256")
published.lambdas.df$Taxa_Bins <- relevel(published.lambdas.df$Taxa_Bins, "65 to 128")
published.lambdas.df$Taxa_Bins <- relevel(published.lambdas.df$Taxa_Bins, "33 to 64")
published.lambdas.df$Taxa_Bins <- relevel(published.lambdas.df$Taxa_Bins, "Up To 32")


length(which(published.lambdas.df$Published_Lambdas > 0.95))/nrow(published.lambdas.pruned)*100 # 18.23%
length(which(published.lambdas.df$Published_Lambdas < 0.05))/nrow(published.lambdas.pruned)*100 # 25.32%

library(RColorBrewer)

jpeg("Figures/Fig3_col.jpeg", res = 120, quality = 100, width = 1000, height = 660)
ggplot(published.lambdas.df, aes(x = Published_Lambdas, fill = Taxa_Bins)) + 
      geom_histogram(binwidth = 0.02) + labs(fill='Tree Size') + 
      theme(legend.position = c(.85,.77), legend.background = element_rect(fill = "white", color = "grey50"),
            panel.background = element_rect(fill = "white"), panel.grid.major = element_line(colour = "grey80")) +
      xlab("Published Lambda Values") + ylab("Frequency") + 
      scale_fill_manual(values = brewer.pal(7, "YlOrRd")[7:1]) 
dev.off()

jpeg("Manuscript/Fig3.jpeg", res = 120, quality = 100, width = 1000, height = 660)
ggplot(published.lambdas.df, aes(x = Published_Lambdas, fill = Taxa_Bins)) + 
      geom_histogram(binwidth = 0.02) + labs(fill='Tree Size') + 
      theme(legend.position = c(.85,.77), legend.background = element_rect(fill = "white", color = "grey50"),
            panel.background = element_rect(fill = "white"), panel.grid.major = element_line(colour = "grey80")) +
      xlab("Published Lambda Values") + ylab("Frequency") + 
      scale_fill_manual(values = brewer.pal(7, "YlOrRd")[7:1]) 
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


## Figure S2 #####

ANOVA.data.ml$lambda.input.bins <- ANOVA.data.ml$lambda.input
ANOVA.data.ml$lambda.input.bins[ANOVA.data.ml$lambda.input < .31] <- "Low"
ANOVA.data.ml$lambda.input.bins[ANOVA.data.ml$lambda.input > .31] <- "Medium"
ANOVA.data.ml$lambda.input.bins[ANOVA.data.ml$lambda.input > .69] <- "High"
ANOVA.data.ml$lambda.input.bins <- as.factor(ANOVA.data.ml$lambda.input.bins)
summary(ANOVA.data.ml$lambda.input.bins)

PlotList <- lapply(1:6, function(j){
   reduced.data <- ANOVA.data.ml[which(ANOVA.data.ml$tree.size == treesizes[[j]]),]
   ggplot(reduced.data, aes(x = beta, y = slope, fill = lambda.input.bins)) + #geom_abline(intercept = 0, slope = .2) +
      geom_boxplot() + theme(legend.position= "none", panel.background = element_rect("transparent")) + 
      xlab("Input Slope") + ylab("Estimated Slope") + ggtitle(treesizes[[j]]) + expand_limits(y=c(-3.5,4))
})

jpeg("Figures/FigS2.jpeg", res = 120, quality = 100, width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()


## Figure S3 ####
library(cowplot)
library(ggplot2)
treesizes <-  2^(5:10)

kappadata_1 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_32.csv")
kappadata_2 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_64.csv")
kappadata_3 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_128.csv")
kappadata_4 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_256.csv")
kappadata_5 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_512.csv")
kappadata_6 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_512.csv")

kappa_all <- rbind(kappadata_1, kappadata_2, kappadata_3, kappadata_4, kappadata_5, kappadata_6)
kappa_all$treesize <- rep(treesizes, each = 1050)


PlotList <- lapply(1:6, function(j){
   reduced.data <- kappa_all[which(kappa_all$treesize == treesizes[[j]]),]
   ggplot(reduced.data, aes(x = lambda.input, y = kappa)) +
      geom_point() + theme(legend.position= "none", panel.background = element_rect("transparent")) + 
      xlab("Input Lambda") + ylab("Estimated Kappa") + ggtitle(treesizes[[j]])
})

jpeg("Figures/FigS3.jpeg", res = 120, quality = 100, width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()



## Sup Figure: ANOVA beta ~ slope  ######

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

# other plot options for Figure 2 - alt1


PlotList <- lapply(1:6, function(j){
   reduced.data <- ANOVA.data[which(ANOVA.data$tree.size == treesizes[[j]]),]
   ggplot(reduced.data, aes(x = beta, y = slope, fill = method)) + #geom_abline(intercept = 0, slope = .2) +
      geom_boxplot() + theme(legend.position= "none", panel.background = element_rect("transparent")) + 
      xlab("Input Slope") + ylab("Estimated Slope") + ggtitle(treesizes[[j]]) + expand_limits(y=c(-3.5,4))
})

jpeg("Figures/Fig2_alt.jpeg", res = 120, quality = 100, width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()


# another plot options for Figure 2 -alt2

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

jpeg("Figures/Fig2_alt2.jpeg", res = 120, quality = 100, width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()





# Rmarkdown PDF rendering ####

rmarkdown::render("Manuscript/Manuscript.Rmd", output_dir = "Manuscript/")
