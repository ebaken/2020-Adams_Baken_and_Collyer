# Figures

data_reg_32 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-32.csv")
data_reg_64 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-64.csv")
data_reg_128 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-128.csv")
data_reg_256 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-256.csv")
data_reg_512 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-512.csv")
data_reg_1028 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-1028.csv")

data_anova_32 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-32-ANOVA.csv")
data_anova_64 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-64-ANOVA.csv")
data_anova_128 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-128-ANOVA.csv")
data_anova_256 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-256-ANOVA.csv")
data_anova_512 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-512-ANOVA.csv")
data_anova_1028 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-1028-ANOVA.csv")

data_reg_32_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-32.csv")
data_reg_64_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-64.csv")
data_reg_128_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-128.csv")
data_reg_256_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-256.csv")
data_reg_512_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-512.csv")
data_reg_1028_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-1028.csv")

data_anova_32_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-32-ANOVA.csv")
data_anova_64_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-64-ANOVA.csv")
data_anova_128_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-128-ANOVA.csv")
data_anova_256_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-256-ANOVA.csv")
data_anova_512_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-512-ANOVA.csv")
data_anova_1028_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-1028-ANOVA.csv")

treesizes <-  2^(5:10)

## Figure 1 ######

library(cowplot)

DataArray <- array(data_reg_32, data_reg_64, data_reg_128,
                   data_reg_256, data_reg_512, data_reg_1028)

PlotList <- lapply(1:6, function(j){
   ggplot(DataArray[[j]], aes(x = lambda.input, y = lambda.est)) +
      geom_point() + theme(legend.position= "none", panel.background = element_rect("transparent")) + 
      xlab("Input Lambda") + ylab("Estimated Lambda") + ggtitle(TreeSizeLabel)
})

jpeg("Figures/Fig1.jpeg", res = 120, quality = 100, width = 1000, height = 660)
plot_grid(PlotList[[1]], PlotList[[2]], PlotList[[3]],
          PlotList[[4]], PlotList[[5]], PlotList[[6]])
dev.off()


## Figure 2 #######







## Figure 3 ######




plot(dataml$'F'~data0$'F')
plot(dataml$'P'~data0$'P')
plot(dataml$slope~data0$slope)


plot(DataTableMat_lambda_0$slope ~ DataTableMat_lambda_0$beta)
