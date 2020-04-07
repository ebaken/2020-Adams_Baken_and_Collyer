# Munging Data
tree.sizes <- 2^(5:10)

# Regression data

data_reg_32 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-32.csv")
data_reg_64 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-64.csv")
data_reg_128 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-128.csv")
data_reg_256 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-256.csv")
data_reg_512 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-512.csv")
data_reg_1024_1 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-1024-1.csv")
data_reg_1024_2 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-1024-2.csv")
data_reg_1024_3 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-1024-3.csv")
data_reg_1024_4 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-1024-4.csv")
data_reg_1024_5 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-1024-5.csv")

data_reg_32_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-32.csv")
data_reg_64_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-64.csv")
data_reg_128_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-128.csv")
data_reg_256_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-256.csv")
data_reg_512_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-512.csv")
data_reg_1024_0_1 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-1024-1.csv")
data_reg_1024_0_2 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-1024-2.csv")
data_reg_1024_0_3 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-1024-3.csv")
data_reg_1024_0_4 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-1024-4.csv")
data_reg_1024_0_5 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-1024-5.csv")


reg.fulldata <- rbind(data_reg_32, data_reg_64, data_reg_128, data_reg_256, data_reg_512, 
                      data_reg_1024_1, data_reg_1024_2, data_reg_1024_3, data_reg_1024_4, data_reg_1024_5,
                      data_reg_32_0, data_reg_64_0, data_reg_128_0, data_reg_256_0, data_reg_512_0, 
                      data_reg_1024_0_1, data_reg_1024_0_2, data_reg_1024_0_3, data_reg_1024_0_4, data_reg_1024_0_5)
# this following line is just a place holder
                reg.fulldata <- rbind(data_reg_32, data_reg_64, data_reg_128, data_reg_256, data_reg_512, data_reg_512,
                                      data_reg_32_0, data_reg_64_0, data_reg_128_0, data_reg_256_0, data_reg_512_0, data_reg_512_0)
#
        
reg.fulldata <- as.data.frame(reg.fulldata)

reg.fulldata$tree.size <- rep(rep(as.character(tree.sizes), each = 5250), 2)
reg.fulldata$method <- rep(c("ml", "0"), each = length(reg.fulldata$tree.size)/2)

write.csv(reg.fulldata, "Data_Analyses/Munged_Data/Regression.fulldataset.csv", row.names = F)

# ANOVA data

data_anova_32 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-32-ANOVA.csv")
data_anova_64 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-64-ANOVA.csv")
data_anova_128 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-128-ANOVA.csv")
data_anova_256 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-256-ANOVA.csv")
data_anova_512 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-512-ANOVA.csv")
data_anova_1024_1 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-1024-ANOVA-1.csv")
data_anova_1024_2 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-1024-ANOVA-2.csv")
data_anova_1024_3 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-1024-ANOVA-3.csv")
data_anova_1024_4 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-1024-ANOVA-4.csv")
data_anova_1024_5 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_ml-1024-ANOVA-5.csv")

data_anova_32_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-32-ANOVA.csv")
data_anova_64_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-64-ANOVA.csv")
data_anova_128_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-128-ANOVA.csv")
data_anova_256_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-256-ANOVA.csv")
data_anova_512_0 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-512-ANOVA.csv")
data_anova_1024_0_1 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-1024-ANOVA-1.csv")
data_anova_1024_0_2 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-1024-ANOVA-2.csv")
data_anova_1024_0_3 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-1024-ANOVA-3.csv")
data_anova_1024_0_4 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-1024-ANOVA-4.csv")
data_anova_1024_0_5 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_0-1024-ANOVA-5.csv")

anova.fulldata <- rbind(data_anova_32, data_anova_64, data_anova_128, data_anova_256, data_anova_512, 
                      data_anova_1024_1, data_anova_1024_2, data_anova_1024_3, data_anova_1024_4, data_anova_1024_5,
                      data_anova_32_0, data_anova_64_0, data_anova_128_0, data_anova_256_0, data_anova_512_0, 
                      data_anova_1024_0_1, data_anova_1024_0_2, data_anova_1024_0_3, data_anova_1024_0_4, data_anova_1024_0_5)
anova.fulldata <- as.data.frame(anova.fulldata)

anova.fulldata$tree.size <- rep(rep(as.character(tree.sizes), each = 5250), 2)
reg.fulldata$method <- rep(c("ml", "0"), each = length(reg.fulldata$tree.size)/2)

write.csv(anova.fulldata, "Data_Analyses/Munged_Data/ANOVA.fulldataset.csv", row.names = F)