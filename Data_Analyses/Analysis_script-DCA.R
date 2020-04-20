# Analyses

treesizes <-  2^(5:10)

ANOVA.data <- read.csv("Munged_Data/ANOVA.fulldataset.csv")
ANOVA.data <- as.data.frame(ANOVA.data)

Regression.data <- read.csv("Munged_Data/Regression.fulldataset.csv")
Regression.data <- as.data.frame(Regression.data)


# Slopes, Ranges, and Variances for Estimated Lambdas, Kappas, and Kappa Zs #####

kappa_data_32 <- read.csv("Sim_Data/PB_lambda_kappacomparison_32.csv")
kappa_data_64 <- read.csv("Sim_Data/PB_lambda_kappacomparison_64.csv")
kappa_data_128 <- read.csv("Sim_Data/PB_lambda_kappacomparison_128.csv")
kappa_data_256 <- read.csv("Sim_Data/PB_lambda_kappacomparison_256.csv")
kappa_data_512 <- read.csv("Sim_Data/PB_lambda_kappacomparison_512.csv")
kappa_data_1024 <- read.csv("Sim_Data/PB_lambda_kappacomparison_512.csv")

# Getting normalized K and Z vals

kappa_data_32$kappa.norm <- (kappa_data_32$kappa-min(kappa_data_32$kappa))/(range(kappa_data_32$kappa)[2] - range(kappa_data_32$kappa)[1])
kappa_data_32$kappa.z.norm <- (kappa_data_32$kappa.z-min(kappa_data_32$kappa.z))/(range(kappa_data_32$kappa.z)[2] - range(kappa_data_32$kappa.z)[1])

kappa_data_64$kappa.norm <- (kappa_data_64$kappa-min(kappa_data_64$kappa))/(range(kappa_data_64$kappa)[2] - range(kappa_data_64$kappa)[1])
kappa_data_64$kappa.z.norm <- (kappa_data_64$kappa.z-min(na.omit(kappa_data_64$kappa.z)))/(range(na.omit(kappa_data_64$kappa.z))[2] - range(na.omit(kappa_data_64$kappa.z))[1])

kappa_data_128$kappa.norm <- (kappa_data_128$kappa-min(kappa_data_128$kappa))/(range(kappa_data_128$kappa)[2] - range(kappa_data_128$kappa)[1])
kappa_data_128$kappa.z.norm <- (kappa_data_128$kappa.z-min(kappa_data_128$kappa.z))/(range(kappa_data_128$kappa.z)[2] - range(kappa_data_128$kappa.z)[1])

kappa_data_256$kappa.norm <- (kappa_data_256$kappa-min(kappa_data_256$kappa))/(range(kappa_data_256$kappa)[2] - range(kappa_data_256$kappa)[1])
kappa_data_256$kappa.z.norm <- (kappa_data_256$kappa.z-min(kappa_data_256$kappa.z))/(range(kappa_data_256$kappa.z)[2] - range(kappa_data_256$kappa.z)[1])

kappa_data_512$kappa.norm <- (kappa_data_512$kappa-min(kappa_data_512$kappa))/(range(kappa_data_512$kappa)[2] - range(kappa_data_512$kappa)[1])
kappa_data_512$kappa.z.norm <- (kappa_data_512$kappa.z-min(kappa_data_512$kappa.z))/(range(kappa_data_512$kappa.z)[2] - range(kappa_data_512$kappa.z)[1])

kappa_data_1024$kappa.norm <- (kappa_data_1024$kappa-min(kappa_data_1024$kappa))/(range(kappa_data_1024$kappa)[2] - range(kappa_data_1024$kappa)[1])
kappa_data_1024$kappa.z.norm <- (kappa_data_1024$kappa.z-min(kappa_data_1024$kappa.z))/(range(kappa_data_1024$kappa.z)[2] - range(kappa_data_1024$kappa.z)[1])

data_list <- list(kappa_data_32, kappa_data_64, kappa_data_128, kappa_data_256, kappa_data_512, kappa_data_1024)


### Some plotting
pairs(kappa_data_128)

##NOTE: I would like these colored by input lambda somehow (or just pull out a few). 
par(mfrow=c(3,2))
hist(kappa_data_32$lambda.est)
hist(kappa_data_512$lambda.est)

hist(kappa_data_32$kappa.norm)
hist(kappa_data_512$kappa.norm)

hist(kappa_data_32$kappa.z.norm)
hist(kappa_data_512$kappa.z.norm)
par(mfrow=c(1,1))

plot(kappa_data_512$lambda.input,kappa_data_512$lambda.est)

# Lambda Estimates 

fit_list <- lapply(1:6, function(i) lm(data_list[[i]]$lambda.est~data_list[[i]]$lambda.input))
slopes <- unlist(lapply(1:6, function(i) coef(fit_list[[i]])[2]))
names(slopes) <- c("32", "64", "128", "256", "512", "1024")

write.csv(slopes, "Data_Analyses/Munged_Data/Lambda_est_slopes.csv", row.names = T)

lambda_input_options <- unique(kappa_data_32$lambda.input)

range_df <- data.frame(lambda_input = lambda_input_options, n32 = rep(NA, 21), n64 = rep(NA, 21), 
                       n128 = rep(NA, 21), n256 = rep(NA, 21), n512 = rep(NA, 21), n1024 = rep(NA, 21))
range_df$n32 <- unlist(lapply(1:21, function(i) max(kappa_data_32[which(kappa_data_32$lambda.input == lambda_input_options[[i]]),2]) - 
                       min(kappa_data_32[which(kappa_data_32$lambda.input == lambda_input_options[[i]]),2])))
range_df$n64 <- unlist(lapply(1:21, function(i) max(kappa_data_64[which(kappa_data_64$lambda.input == lambda_input_options[[i]]),2]) - 
                       min(kappa_data_64[which(kappa_data_64$lambda.input == lambda_input_options[[i]]),2])))
range_df$n128 <- unlist(lapply(1:21, function(i) max(kappa_data_128[which(kappa_data_128$lambda.input == lambda_input_options[[i]]),2]) - 
                       min(kappa_data_128[which(kappa_data_128$lambda.input == lambda_input_options[[i]]),2])))
range_df$n256 <- unlist(lapply(1:21, function(i) max(kappa_data_256[which(kappa_data_256$lambda.input == lambda_input_options[[i]]),2]) - 
                       min(kappa_data_256[which(kappa_data_256$lambda.input == lambda_input_options[[i]]),2])))
range_df$n512 <- unlist(lapply(1:21, function(i) max(kappa_data_512[which(kappa_data_512$lambda.input == lambda_input_options[[i]]),2]) - 
                       min(kappa_data_512[which(kappa_data_512$lambda.input == lambda_input_options[[i]]),2])))
range_df$n1024 <- unlist(lapply(1:21, function(i) max(kappa_data_1024[which(kappa_data_1024$lambda.input == lambda_input_options[[i]]),2]) - 
                       min(kappa_data_1024[which(kappa_data_1024$lambda.input == lambda_input_options[[i]]),2])))


write.csv(range_df, "Munged_Data/Lambda_est_ranges.csv", row.names = F)


var_df <- data.frame(lambda_input = lambda_input_options, n32 = rep(NA, 21), n64 = rep(NA, 21), 
                       n128 = rep(NA, 21), n256 = rep(NA, 21), n512 = rep(NA, 21), n1024 = rep(NA, 21))
var_df$n32 <- unlist(lapply(1:21, function(i) var(kappa_data_32[which(kappa_data_32$lambda.input == lambda_input_options[[i]]),2])))
var_df$n64 <- unlist(lapply(1:21, function(i) var(kappa_data_64[which(kappa_data_64$lambda.input == lambda_input_options[[i]]),2])))
var_df$n128 <- unlist(lapply(1:21, function(i) var(kappa_data_128[which(kappa_data_128$lambda.input == lambda_input_options[[i]]),2])))
var_df$n256 <- unlist(lapply(1:21, function(i) var(kappa_data_256[which(kappa_data_256$lambda.input == lambda_input_options[[i]]),2])))
var_df$n512 <- unlist(lapply(1:21, function(i) var(kappa_data_512[which(kappa_data_512$lambda.input == lambda_input_options[[i]]),2])))
var_df$n1024 <- unlist(lapply(1:21, function(i) var(kappa_data_1024[which(kappa_data_1024$lambda.input == lambda_input_options[[i]]),2])))

var_df

write.csv(var_df, "Munged_Data/Lambda_est_vars.csv", row.names = F)

# Kappa Estimates

lambda_input_options <- unique(kappa_data_32$lambda.input)

range_df <- data.frame(lambda_input = lambda_input_options, n32 = rep(NA, 21), n64 = rep(NA, 21), 
                       n128 = rep(NA, 21), n256 = rep(NA, 21), n512 = rep(NA, 21), n1024 = rep(NA, 21))

range_df$n32 <- unlist(lapply(1:21, function(i) max(kappa_data_32[which(kappa_data_32$lambda.input == lambda_input_options[[i]]),5]) - 
                       min(kappa_data_32[which(kappa_data_32$lambda.input == lambda_input_options[[i]]),5])))
range_df$n64 <- unlist(lapply(1:21, function(i) max(kappa_data_64[which(kappa_data_64$lambda.input == lambda_input_options[[i]]),5]) - 
                       min(kappa_data_64[which(kappa_data_64$lambda.input == lambda_input_options[[i]]),5])))
range_df$n128 <- unlist(lapply(1:21, function(i) max(kappa_data_128[which(kappa_data_128$lambda.input == lambda_input_options[[i]]),5]) - 
                       min(kappa_data_128[which(kappa_data_128$lambda.input == lambda_input_options[[i]]),5])))
range_df$n256 <- unlist(lapply(1:21, function(i) max(kappa_data_256[which(kappa_data_256$lambda.input == lambda_input_options[[i]]),5]) - 
                       min(kappa_data_256[which(kappa_data_256$lambda.input == lambda_input_options[[i]]),5])))
range_df$n512 <- unlist(lapply(1:21, function(i) max(kappa_data_512[which(kappa_data_512$lambda.input == lambda_input_options[[i]]),5]) - 
                       min(kappa_data_512[which(kappa_data_512$lambda.input == lambda_input_options[[i]]),5])))
range_df$n1024 <- unlist(lapply(1:21, function(i) max(kappa_data_1024[which(kappa_data_1024$lambda.input == lambda_input_options[[i]]),5]) - 
                       min(kappa_data_1024[which(kappa_data_1024$lambda.input == lambda_input_options[[i]]),5])))


write.csv(range_df, "Munged_Data/Kappa_norm_ranges.csv", row.names = F)


var_df <- data.frame(lambda_input = lambda_input_options, n32 = rep(NA, 21), n64 = rep(NA, 21), 
                       n128 = rep(NA, 21), n256 = rep(NA, 21), n512 = rep(NA, 21), n1024 = rep(NA, 21))
var_df$n32 <- unlist(lapply(1:21, function(i) var(kappa_data_32[which(kappa_data_32$lambda.input == lambda_input_options[[i]]),5])))
var_df$n64 <- unlist(lapply(1:21, function(i) var(kappa_data_64[which(kappa_data_64$lambda.input == lambda_input_options[[i]]),5])))
var_df$n128 <- unlist(lapply(1:21, function(i) var(kappa_data_128[which(kappa_data_128$lambda.input == lambda_input_options[[i]]),5])))
var_df$n256 <- unlist(lapply(1:21, function(i) var(kappa_data_256[which(kappa_data_256$lambda.input == lambda_input_options[[i]]),5])))
var_df$n512 <- unlist(lapply(1:21, function(i) var(kappa_data_512[which(kappa_data_512$lambda.input == lambda_input_options[[i]]),5])))
var_df$n1024 <- unlist(lapply(1:21, function(i) var(kappa_data_1024[which(kappa_data_1024$lambda.input == lambda_input_options[[i]]),5])))

var_df

write.csv(var_df, "Munged_Data/Kappa_norm_vars.csv", row.names = F)

# Kappa Z Norms

range_df <- data.frame(lambda_input = lambda_input_options, n32 = rep(NA, 21), n64 = rep(NA, 21), 
                       n128 = rep(NA, 21), n256 = rep(NA, 21), n512 = rep(NA, 21), n1024 = rep(NA, 21))

range_df$n32 <- unlist(lapply(1:21, function(i) max(kappa_data_32[which(kappa_data_32$lambda.input == lambda_input_options[[i]]),6]) - 
                       min(kappa_data_32[which(kappa_data_32$lambda.input == lambda_input_options[[i]]),6])))
range_df$n64 <- unlist(lapply(1:21, function(i) max(kappa_data_64[which(kappa_data_64$lambda.input == lambda_input_options[[i]]),6]) - 
                       min(kappa_data_64[which(kappa_data_64$lambda.input == lambda_input_options[[i]]),6])))
range_df$n128 <- unlist(lapply(1:21, function(i) max(kappa_data_128[which(kappa_data_128$lambda.input == lambda_input_options[[i]]),6]) - 
                       min(kappa_data_128[which(kappa_data_128$lambda.input == lambda_input_options[[i]]),6])))
range_df$n256 <- unlist(lapply(1:21, function(i) max(kappa_data_256[which(kappa_data_256$lambda.input == lambda_input_options[[i]]),6]) - 
                       min(kappa_data_256[which(kappa_data_256$lambda.input == lambda_input_options[[i]]),6])))
range_df$n512 <- unlist(lapply(1:21, function(i) max(kappa_data_512[which(kappa_data_512$lambda.input == lambda_input_options[[i]]),6]) - 
                       min(kappa_data_512[which(kappa_data_512$lambda.input == lambda_input_options[[i]]),6])))
range_df$n1024 <- unlist(lapply(1:21, function(i) max(kappa_data_1024[which(kappa_data_1024$lambda.input == lambda_input_options[[i]]),6]) - 
                       min(kappa_data_1024[which(kappa_data_1024$lambda.input == lambda_input_options[[i]]),6])))

range_df

write.csv(range_df, "Munged_Data/Kappa_Z_norm_ranges.csv", row.names = F)


var_df <- data.frame(lambda_input = lambda_input_options, n32 = rep(NA, 21), n64 = rep(NA, 21), 
                       n128 = rep(NA, 21), n256 = rep(NA, 21), n512 = rep(NA, 21), n1024 = rep(NA, 21))
var_df$n32 <- unlist(lapply(1:21, function(i) var(kappa_data_32[which(kappa_data_32$lambda.input == lambda_input_options[[i]]),6])))
var_df$n64 <- unlist(lapply(1:21, function(i) var(kappa_data_64[which(kappa_data_64$lambda.input == lambda_input_options[[i]]),6])))
var_df$n128 <- unlist(lapply(1:21, function(i) var(kappa_data_128[which(kappa_data_128$lambda.input == lambda_input_options[[i]]),6])))
var_df$n256 <- unlist(lapply(1:21, function(i) var(kappa_data_256[which(kappa_data_256$lambda.input == lambda_input_options[[i]]),6])))
var_df$n512 <- unlist(lapply(1:21, function(i) var(kappa_data_512[which(kappa_data_512$lambda.input == lambda_input_options[[i]]),6])))
var_df$n1024 <- unlist(lapply(1:21, function(i) var(kappa_data_1024[which(kappa_data_1024$lambda.input == lambda_input_options[[i]]),6])))

var_df

write.csv(var_df, "Munged_Data/Kappa_Z_norm_vars.csv", row.names = F)




# Literature numbers ####
library(readxl)
review <- read_excel("Literature_Review/1Summary.xlsx", sheet = "Since 2019 Lambda Vals")
colnames(review)
length(unique(review$Reference)) #182

# Published lambda values
published.lambdas <- review[,c(1,3,4)]
published.lambdas.pruned <- published.lambdas[-which(published.lambdas$`Estimated Lambdas Simple` > 1 | published.lambdas$`Estimated Lambdas Simple` < 0),]
published.lambdas.df <- as.data.frame(published.lambdas.pruned)
nrow(published.lambdas.df) # 1552
1552/182 # average number of reported lambda values per manuscript
colnames(published.lambdas.df) <- c("Refs", "Published_Lambdas", "Taxa")
published.lambdas.df$Refs <- as.factor(published.lambdas.df$Refs)

min(table(published.lambdas.df$Refs)) # 1: lowest number of reported lambdas per paper
max(table(published.lambdas.df$Refs)) # 71: highest number of reported lambdas per paper

length(which(published.lambdas.df$Published_Lambdas < 0.05))/nrow(published.lambdas.df)*100 # 25.43%
length(which(published.lambdas.df$Published_Lambdas > 0.90))/nrow(published.lambdas.df)*100 # 24.74%

length(which(published.lambdas.df$Taxa < 30))/nrow(published.lambdas.df)*100 # 75.32%
length(which(published.lambdas.df$Taxa < 30)) # 348




# Numerical summaries of all papers
review.fullsheet <- read_excel("Literature_Review/1Summary.xlsx", sheet = "Since 2019")
review.fullsheet <- as.data.frame(review.fullsheet)

all.papers <- review.fullsheet[,c("Reference", "Estimated lambdas")]
all.papers$YN <- all.papers$`Estimated lambdas`

all.papers$YN[startsWith(all.papers$'Estimated lambdas', "Yes")] <- "Yes"
all.papers$YN[startsWith(all.papers$'Estimated lambdas', "No")] <- "No"

# cleaning up messier answers
all.papers$YN[which(all.papers$'Estimated lambdas' == "-")] <- "No"
all.papers$YN[which(all.papers$'Estimated lambdas' == "no")] <- "No"
all.papers$YN[which(all.papers$'Estimated lambdas' == "Used Blomberg's K and Pagels lambda")] <- "Yes"

all.papers$YN <- as.factor(all.papers$YN)
levels(all.papers$YN)
summary(all.papers$YN) # 204 yeses out of 341 citations


# taking out any rows with NA in the interpretation column
interpreted.pruned <- interpreted[-which(is.na(interpreted$Interp)=="TRUE"),]

# making a column to simplify these answers
interpreted.pruned$YNKinda <- interpreted.pruned$Interp
interpreted.pruned$YNKinda[startsWith(interpreted.pruned$Interp, "Yes")] <- "Yes"
interpreted.pruned$YNKinda[startsWith(interpreted.pruned$Interp, "No")] <- "No"
interpreted.pruned$YNKinda[startsWith(interpreted.pruned$Interp, "Kinda")] <- "Kinda"

# cleaning up messier categories
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "-")] <- NA
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "no")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Only interpreted kappa")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Only interpreted significance")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Only kappa interpreted")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Lambda not listed in results")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Lambda and kappa have conflicting results regarding significance")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Only interpreted significant")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Only published in SI")] <- "Yes"

interpreted.pruned$YNKinda <- as.factor(interpreted.pruned$YNKinda)
levels(interpreted.pruned$YNKinda)

length(which(interpreted.pruned$YNKinda == "Yes"))/nrow(interpreted.pruned)*100 # 20.49%


# comparing lambdas??
interpreted<- review.fullsheet[,c("Reference", "Compared lambdas without sig tests?")]
colnames(interpreted) <- c("Ref", "Compared")
interpreted$Compared <- as.factor(interpreted$Compared)
summary(interpreted$Compared) # 7 yeses

# Misspecification rates - delete

Regression.lambda0 <- Regression.data[which(Regression.data$lambda.input == 0),]
Regression.lambda0 <- Regression.lambda0[Regression.lambda0$method == "ml",]

Reg.l0.listbytree <- list(Regression.lambda0[Regression.lambda0$tree.size == treesizes[[1]],],
                          Regression.lambda0[Regression.lambda0$tree.size == treesizes[[2]],],
                          Regression.lambda0[Regression.lambda0$tree.size == treesizes[[3]],],
                          Regression.lambda0[Regression.lambda0$tree.size == treesizes[[4]],],
                          Regression.lambda0[Regression.lambda0$tree.size == treesizes[[5]],],
                          Regression.lambda0[Regression.lambda0$tree.size == treesizes[[6]],])

summary(Reg.l0.listbytree[[1]]$lambda.input)

lapply(1:6, function(i) { length(which(Reg.l0.listbytree[[i]]$lambda.est.x > 0.1))/nrow(Reg.l0.listbytree[[i]])*100 })
lapply(1:6, function(i) { length(which(Reg.l0.listbytree[[i]]$P < 0.05))/nrow(Reg.l0.listbytree[[i]])*100 })
