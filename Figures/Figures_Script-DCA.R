# Figures


# Fig 1: lit hist ####

library(readxl)

tree.sizes <-  treesizes <- 2^(5:10)

review <- read_excel("Literature_Review/1Summary.xlsx", sheet = "Since 2019 Lambda Vals", col_types = "numeric")
colnames(review)
published.lambdas <- review[,c(3,4)]
published.lambdas.pruned <- published.lambdas[-which(published.lambdas > 1 | published.lambdas < 0),]
published.lambdas.df <- as.data.frame(published.lambdas.pruned)
colnames(published.lambdas.df) <- c("Published_Lambdas", "Taxa")

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

length(which(published.lambdas.df$Published_Lambdas > 0.95))/nrow(published.lambdas.pruned)*100 # 18.509%
length(which(published.lambdas.df$Published_Lambdas < 0.05))/nrow(published.lambdas.pruned)*100 # 25.321%

library(RColorBrewer)
library(ggplot2)

tiff("Figures/Fig1.tiff", width = 1000, height = 660, pointsize = 16)
ggplot(published.lambdas.df, aes(x = Published_Lambdas, fill = Taxa_Bins)) + 
      geom_histogram(binwidth = 0.02) + labs(fill='Tree Size') + 
      theme(legend.position = c(.85,.77), legend.background = element_rect(fill = "white", color = "grey50"),
            panel.background = element_rect(fill = "white"), panel.grid.major = element_line(colour = "grey80"),
            text = element_text(size = 16)) +
      xlab("Published Lambda Values") + ylab("Frequency") + 
      scale_fill_manual(values = brewer.pal(7, "YlOrRd")[7:1]) 
dev.off()

png("Figures/Fig1.png", width = 1000, height = 660, pointsize = 16)
ggplot(published.lambdas.df, aes(x = Published_Lambdas, fill = Taxa_Bins)) + 
      geom_histogram(binwidth = 0.02) + labs(fill='Tree Size') + 
      theme(legend.position = c(.85,.77), legend.background = element_rect(fill = "white", color = "grey50"),
            panel.background = element_rect(fill = "white"), panel.grid.major = element_line(colour = "grey80"),
            text = element_text(size = 16)) +
      xlab("Published Lambda Values") + ylab("Frequency") + 
      scale_fill_manual(values = brewer.pal(7, "YlOrRd")[7:1]) 
dev.off()


# Fig 2: est~input at each sample size ####

library(cowplot)
library(ggplot2)

Data32 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_32.csv")
Data64 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_64.csv")
Data128 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_128.csv")
Data256 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_256.csv")
Data512 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_512.csv")
Data1024_1 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_1024_1.csv")
Data1024_2 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_1024_2.csv")
Data1024 <- rbind(Data1024_1,Data1024_2)

data_list <- list(Data32, Data64, Data128,
                  Data256, Data512, Data1024)

PlotList <- lapply(1:6, function(j){
   ggplot(data_list[[j]], aes(x = lambda.input, y = lambda.est)) +
      geom_point() + theme(legend.position= "none", panel.background = element_rect("transparent"),
                           panel.border = element_rect(color = "black", fill = NA), text = element_text(size = 12)) + 
      xlab("Input Lambda") + ylab("Estimated Lambda") + ggtitle(treesizes[[j]]) 
})

tiff("Figures/Fig2.tiff", width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()

png("Figures/Fig2.png", width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()

# Fig 3: same thing with regression beta = 0.5 ####
Data32 <- read.csv("Data_Analyses/Sim_Data/PB_reg-32.csv")
Data64 <- read.csv("Data_Analyses/Sim_Data/PB_reg-64.csv")
Data128 <- read.csv("Data_Analyses/Sim_Data/PB_reg-128.csv")
Data256 <- read.csv("Data_Analyses/Sim_Data/PB_reg-256.csv")
Data512_1 <- read.csv("Data_Analyses/Sim_Data/PB_reg-512-1.csv")
Data512_2 <- read.csv("Data_Analyses/Sim_Data/PB_reg-512-2.csv")
Data512 <- rbind(Data512_1, Data512_2)
Data1024_1 <- read.csv("Data_Analyses/Sim_Data/PB_reg-1024-1.csv")
Data1024_2 <- read.csv("Data_Analyses/Sim_Data/PB_reg-1024-2.csv")
Data1024_3 <- read.csv("Data_Analyses/Sim_Data/PB_reg-1024-3.csv")
Data1024_4 <- read.csv("Data_Analyses/Sim_Data/PB_reg-1024-4.csv")
Data1024_5 <- read.csv("Data_Analyses/Sim_Data/PB_reg-1024-5.csv")
Data1024 <- rbind(Data1024_1,Data1024_2,Data1024_3,Data1024_4,Data1024_5)

Data32_pruned <- Data32[which(Data32$beta == 0.5),]
Data64_pruned <- Data64[which(Data64$beta == 0.5),]
Data128_pruned <- Data128[which(Data128$beta == 0.5),]
Data256_pruned <- Data256[which(Data256$beta == 0.5),]
Data512_pruned <- Data512[which(Data512$beta == 0.5),]
Data1024_pruned <- Data1024[which(Data1024$beta == 0.5),]

data_list <- list(Data32_pruned, Data64_pruned, Data128_pruned,
                  Data256_pruned, Data512_pruned, Data1024_pruned)

PlotList <- lapply(1:6, function(j){
   ggplot(data_list[[j]], aes(x = lambda.input, y = lambda.est.y)) +
      geom_point() + theme(legend.position= "none", panel.background = element_rect("transparent"),
                           panel.border = element_rect(color = "black", fill = NA), text = element_text(size = 12)) + 
      xlab("Input Lambda") + ylab("Estimated Lambda") + ggtitle(treesizes[[j]]) 
})

tiff("Figures/Fig3.tiff", width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()

png("Figures/Fig3.png", width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()

# Fig 4: panel A: z score of lambda ~ lambda input at n = 128, ####
  #REWORKING FOR FIG 4

# READ IN ALL DATA
Data32 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_32.csv")
Data64 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_64.csv")
Data128 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_128.csv")
Data256 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_256.csv")
Data512 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_512.csv")
Data1024 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_1024.csv")


###LOOP THIS OVER EACH LEVEL OF N. 
        # calculating z scores for lambdas 

# add loop here over all tree sizes
        Z.lambda <- matrix(NA, nrow=50, ncol = 21)
        n <- 128
        nsim <- 50
        lambdas <- levels(as.factor(Data128$lambda.input))
        
        for (i in 1:21) {
            Data128_pruned <- Data128[which(Data128$lambda.input == lambdas[i]),]
            
            lambda.std <- unlist(lapply(1:nsim, function(j) {
                if(Data128_pruned$CI.lower[j]==0){lambda.std <- (Data128_pruned$CI.upper[j]-Data128_pruned$lambda.est[j])*sqrt(n)/qnorm(0.975)} 
                else {lambda.std <- (Data128_pruned$lambda.est[j]-Data128_pruned$CI.lower[j])*sqrt(n)/qnorm(0.975)}
              }))
            
            Z.lambda[,i] <- Data128_pruned$lambda.est/lambda.std
        }
        
        Data128$lambda.est.z <- as.vector(Z.lambda)
        
        # end loop

#I THINK for Panel A & B, use N=32

PanelA <- ggplot(Data128, aes(x = lambda.input, y = lambda.est.z)) + 
                geom_point() + theme(legend.position= "none", panel.background = element_rect("transparent"),
                                         panel.border = element_rect(color = "black", fill = NA), text = element_text(size = 12)) + 
                    xlab("Input Lambda") + ylab("Estimated Lambda Effect Size") + ggtitle("A") 

# panel B: z score for kappa ~ lambda at n = 128

PanelB <- ggplot(Data128, aes(x = lambda.input, y = kappa.z)) + 
                geom_point() + theme(legend.position= "none", panel.background = element_rect("transparent"),
                                         panel.border = element_rect(color = "black", fill = NA), text = element_text(size = 12)) + 
                    xlab("Input Lambda") + ylab("Estimated Kappa Effect Size") + ggtitle("B")


# panel C: coefficient of variation across precision estimates for each sample size ~ n, 
# go to all 21 lambda inputs and get variance of lambda z scores, then calculate coefficient of variation of those variances
      
# LOOP THIS OVER EACH LEVEL OF N.  SO obtain 21 lambda.z.vars at each
      lambda.z.vars <- apply(Z.lambda,2,var) # loop this over all tree sizes

#Then for each N do:
      CV.k <- sd(apply(Z.k,2,var))/   mean(apply(Z.k,2,var))*100 # loop that over all tree sizes

# THEN PLOT all 6 CV for K and all 6 CV for lambda.  (Hmm, perhaps this is a single plot?)      

      barplot(lambda.z.vars)
      library(plotrix)
      from <- 1
      to <- 571


gap.barplot(lambda.z.vars, gap=c(from,to), col = rep("darkgray",21), ytics = c(0,1,572, 573), 
                  xaxlab = lambdas, 
                  xlab = "Lambda Input", ylab = "Lambda Z Score Variances", yaxs="i")
axis.break(2, from, breakcol="snow", style="gap")
axis.break(2, from*(1+0.04), breakcol="black", style="slash")
axis.break(4, from*(1+0.04), breakcol="black", style="slash")
axis(2, at=from) 
title(main = list("C", cex = 1, font = 1), adj = 0)

      CV.l <- sd(apply(Z.lambda,2,var))/ mean(apply(Z.lambda,2,var))*100

# panel D: same as panel C except kappa z scores [he is sending me code for this], but barplot all of cv of z of lambda, then all cv of z of kappa
      
      Z.k <- matrix(Data128$kappa.z, nrow = 50, ncol = 21)
      colnames(Z.k) <- lambdas
      
      lambda.k.vars <- apply(Z.k,2,var) # loop over all tree sizes

barplot(lambda.k.vars, col = "darkgray", space = -0.01, axis.lty = 1, ylim = c(0,max(lambda.k.vars)),
        xlab = "Lambda Input", ylab = "Lambda Kappa Score Variances")
box()
title(main = list("D", cex = 1, font = 1), adj = 0)
      
      CV.k <- sd(apply(Z.k,2,var))/   mean(apply(Z.k,2,var))*100 # loop over all tree sizes
      
      barplot(c(CV.l,CV.k))
      
      # bar plot with tree sizes on x axis and CV for lambda z and kappa z side by side
      
# All together
library(gridExtra)
png("Figures/Fig4_AB.png", width = 600, height = 300)
grid.arrange(PanelA, PanelB, nrow=1)
dev.off()

png("Figures/Fig4_CD.png", width = 600, height = 300)
par(mfrow=c(1,2), mar = c(5,4,2,2))
gap.barplot(lambda.z.vars, gap=c(from,to), col = rep("darkgray",21), ytics = c(0,1,572, 573), 
                  xaxlab = lambdas, 
                  xlab = "Lambda Input", ylab = "Lambda Z Score Variances", yaxs="i")
axis.break(2, from, breakcol="snow", style="gap")
axis.break(2, from*(1+0.04), breakcol="black", style="slash")
axis.break(4, from*(1+0.04), breakcol="black", style="slash")
axis(2, at=from) 
title(main = list("C", cex = 1, font = 1), adj = 0)

barplot(lambda.k.vars, col = "darkgray", space = -0.01, axis.lty = 1, ylim = c(0,max(lambda.k.vars)),
        xlab = "Lambda Input", ylab = "Lambda Kappa Score Variances")
box()
title(main = list("D", cex = 1, font = 1), adj = 0)

dev.off()

# Fig 5: Salamander fig ####

# Panel A made in ppt, panels B and C made in Analysis_script

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
