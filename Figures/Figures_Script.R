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

# Fig 4: Edited April 28 2020 ####

  Data32 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_32.csv")
  Data64 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_64.csv")
  Data128 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_128.csv")
  Data256 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_256.csv")
  Data512 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_512.csv")
  Data1024 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_1024.csv")

        # calculating z scores for lambdas 
        tree.sizes <-  2^(5:10)
        lambdas <- levels(as.factor(Data128$lambda.input))
        nsim <- 50
        data_list <- list(Data32, Data64, Data128, Data256, Data512, Data1024)
        Z.lambda <- array(NA, dim = c(50,21,6))
        
        for (k in 1:6) {
          
            n <- tree.sizes[k]
            data <- data_list[[k]]
        
            for (i in 1:21) {
                data_pruned <- data[which(data$lambda.input == lambdas[i]),]
                  lambda.std <- unlist(lapply(1:nsim, function(j) {
                      if(data_pruned$CI.lower[j]==0 & data_pruned$CI.upper[j]==1) {lambda.std <- NA} # skipping all values for which CIs were (0,1)
                    else{
                      if(data_pruned$CI.lower[j]==0){lambda.std <- (data_pruned$CI.upper[j]-data_pruned$lambda.est[j])*sqrt(n)/qnorm(0.975)} 
                      else {lambda.std <- (data_pruned$lambda.est[j]-data_pruned$CI.lower[j])*sqrt(n)/qnorm(0.975)}
                    }
                    }))
                 Z.lambda[,i,k] <- data_pruned$lambda.est/lambda.std
            }
        }
        Z.lambda[,,1]

# Panel A: lambda.z ~ lambda.input at n = 32
        # put data from above loop in long format
        library(tidyr)
        Data <- as.data.frame(Z.lambda[,,1])
        colnames(Data) <- lambdas
        Data_long <- gather(Data, key = "lambda.input", value = "lambda.est.z")
        Data_long$lambda.input <- as.numeric(Data_long$lambda.input)

PanelA <- ggplot(Data_long, aes(x = lambda.input, y = lambda.est.z)) + 
                geom_point() + theme(legend.position= "none", panel.background = element_rect("transparent"),
                                         panel.border = element_rect(color = "black", fill = NA), text = element_text(size = 12)) + 
                    xlab("Input Lambda") + ylab("Estimated Lambda Effect Size") 
        # warning about removing missing values is ok, those are the CIs that were (0,1) and thus Lambda Z could not be calculated

# Panel B: kappa.z ~ lambda.input at n = 32

PanelB <- ggplot(Data32, aes(x = lambda.input, y = kappa.z)) + 
                geom_point() + theme(legend.position= "none", panel.background = element_rect("transparent"),
                                         panel.border = element_rect(color = "black", fill = NA), text = element_text(size = 12)) + 
                    xlab("Input Lambda") + ylab("Estimated Kappa Effect Size")


# panel C: coefficient of variation across precision estimates for each sample size ~ n, 
      lambda.z.vars <- lapply(1:6, function(i) { apply(na.omit(Z.lambda[,,i]), 2, var) }) # variances of z scores across each input lambda for each tree
      CV.lambda <- lapply(1:6, function(i) sd(lambda.z.vars[[i]])/mean(lambda.z.vars[[i]])*100) # CV for each tree size

      data_list <- list(Data32, Data64, Data128, Data256, Data512, Data1024)
      lambda.k.vars <- lapply(1:6, function(i) { 
                    data_wide <- matrix(data_list[[i]]$kappa.z, nrow = 50, byrow = F)
                    apply(data_wide, 2, var) 
                    })
      CV.kappa <- lapply(1:6, function(i) sd(lambda.k.vars[[i]])/mean(lambda.k.vars[[i]])*100)
      # CV.kappa <- sd(apply(Z.lambda,2,var))/ mean(apply(Z.lambda,2,var))*100 # Dean's code
      
      DF <- data.frame(n = rep(tree.sizes, 2), CV = c(unlist(CV.lambda),unlist(CV.kappa)), Statistic = rep(c("Lambda Z", "Kappa Z"), each= 6))

PanelC <- ggplot(DF, aes(x = as.factor(n), y = CV, fill = Statistic)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal() +
  xlab("Tree Size") + ylab("Coefficient of Variance") + 
  scale_fill_manual(values=c("white","darkgray"))

PanelC      
      

      
# All together
library(gridExtra)
lay <- rbind(c(1,1,2,2),
             c(1,1,2,2),
             c(3,3,3,3))
gs <- list(PanelA,PanelB,PanelC)

png("Figures/Fig4.png", width = 600, height = 600)
grid.arrange(grobs = gs, layout_matrix = lay)
dev.off()


# Fig 5: Salamander fig ####

# Panel A made in ppt, panels B and C made in Analysis_script

# Rmarkdown PDF rendering ####

rmarkdown::render("Manuscript/Manuscript.Rmd", output_dir = "Manuscript/")
