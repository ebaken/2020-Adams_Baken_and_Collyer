# Analyses

treesizes <-  2^(5:10)

ANOVA.data <- read.csv("Data_Analyses/Munged_Data/ANOVA.fulldataset.csv")
ANOVA.data <- as.data.frame(ANOVA.data)

Regression.data <- read.csv("Data_Analyses/Munged_Data/Regression.fulldataset.csv")
Regression.data <- as.data.frame(Regression.data)


# Slopes, Ranges, and Variances for Estimated Lambdas, Kappas, and Kappa Zs #####

kappa_data_32 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_32.csv")
kappa_data_64 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_64.csv")
kappa_data_128 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_128.csv")
kappa_data_256 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_256.csv")
kappa_data_512 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_512.csv")
kappa_data_1024 <- read.csv("Data_Analyses/Sim_Data/PB_lambda_kappacomparison_1024.csv")

# Getting normalized K and Z vals

kappa_data_32$kappa.norm <- (kappa_data_32$kappa-min(kappa_data_32$kappa))/(range(kappa_data_32$kappa)[2] - range(kappa_data_32$kappa)[1])
kappa_data_32$kappa.z.norm <- (kappa_data_32$kappa.z-min(kappa_data_32$kappa.z))/(range(kappa_data_32$kappa.z)[2] - range(kappa_data_32$kappa.z)[1])

kappa_data_64$kappa.norm <- (kappa_data_64$kappa-min(kappa_data_64$kappa))/(range(kappa_data_64$kappa)[2] - range(kappa_data_64$kappa)[1])
kappa_data_64$kappa.z.norm <- (kappa_data_64$kappa.z-min(kappa_data_64$kappa.z))/(range(kappa_data_64$kappa.z)[2] - range(kappa_data_64$kappa.z)[1])

kappa_data_128$kappa.norm <- (kappa_data_128$kappa-min(kappa_data_128$kappa))/(range(kappa_data_128$kappa)[2] - range(kappa_data_128$kappa)[1])
kappa_data_128$kappa.z.norm <- (kappa_data_128$kappa.z-min(kappa_data_128$kappa.z))/(range(kappa_data_128$kappa.z)[2] - range(kappa_data_128$kappa.z)[1])

kappa_data_256$kappa.norm <- (kappa_data_256$kappa-min(kappa_data_256$kappa))/(range(kappa_data_256$kappa)[2] - range(kappa_data_256$kappa)[1])
kappa_data_256$kappa.z.norm <- (kappa_data_256$kappa.z-min(kappa_data_256$kappa.z))/(range(kappa_data_256$kappa.z)[2] - range(kappa_data_256$kappa.z)[1])

kappa_data_512$kappa.norm <- (kappa_data_512$kappa-min(kappa_data_512$kappa))/(range(kappa_data_512$kappa)[2] - range(kappa_data_512$kappa)[1])
kappa_data_512$kappa.z.norm <- (kappa_data_512$kappa.z-min(kappa_data_512$kappa.z))/(range(kappa_data_512$kappa.z)[2] - range(kappa_data_512$kappa.z)[1])

kappa_data_1024$kappa.norm <- (kappa_data_1024$kappa-min(kappa_data_1024$kappa))/(range(kappa_data_1024$kappa)[2] - range(kappa_data_1024$kappa)[1])
kappa_data_1024$kappa.z.norm <- (kappa_data_1024$kappa.z-min(kappa_data_1024$kappa.z))/(range(kappa_data_1024$kappa.z)[2] - range(kappa_data_1024$kappa.z)[1])


data_list <- list(kappa_data_32, kappa_data_64, kappa_data_128, kappa_data_256, kappa_data_512, kappa_data_1024)

# Lambda Estimates 

fit_list <- lapply(1:6, function(i) lm(data_list[[i]]$lambda.est~data_list[[i]]$lambda.input))
conf_int <- lapply(1:6, function(i) confint(fit_list[[i]], 'data_list[[i]]$lambda.input', level=0.95))
slopes <- unlist(lapply(1:6, function(i) coef(fit_list[[i]])[2]))
names(slopes) <- c("32", "64", "128", "256", "512", "1024")

write.csv(slopes, "Data_Analyses/Munged_Data/Lambda_est__input_slope.csv", row.names = T)

# y ranges
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

write.csv(range_df, "Data_Analyses/Munged_Data/Lambda_est_yranges.csv", row.names = F)

# y var

var_df <- data.frame(lambda_input = lambda_input_options, n32 = rep(NA, 21), n64 = rep(NA, 21), 
                       n128 = rep(NA, 21), n256 = rep(NA, 21), n512 = rep(NA, 21), n1024 = rep(NA, 21))
var_df$n32 <- unlist(lapply(1:21, function(i) var(kappa_data_32[which(kappa_data_32$lambda.input == lambda_input_options[[i]]),2])))
var_df$n64 <- unlist(lapply(1:21, function(i) var(kappa_data_64[which(kappa_data_64$lambda.input == lambda_input_options[[i]]),2])))
var_df$n128 <- unlist(lapply(1:21, function(i) var(kappa_data_128[which(kappa_data_128$lambda.input == lambda_input_options[[i]]),2])))
var_df$n256 <- unlist(lapply(1:21, function(i) var(kappa_data_256[which(kappa_data_256$lambda.input == lambda_input_options[[i]]),2])))
var_df$n512 <- unlist(lapply(1:21, function(i) var(kappa_data_512[which(kappa_data_512$lambda.input == lambda_input_options[[i]]),2])))
var_df$n1024 <- unlist(lapply(1:21, function(i) var(kappa_data_1024[which(kappa_data_1024$lambda.input == lambda_input_options[[i]]),2])))

var_df

write.csv(var_df, "Data_Analyses/Munged_Data/Lambda_est_vars.csv", row.names = F)

# x ranges and var

range_df <- data.frame(n32_min = rep(NA, 11), n32_max = rep(NA, 11), n32_var = rep(NA, 11), 
                       n64_min = rep(NA, 11), n64_max = rep(NA, 11), n64_var = rep(NA, 11),
                       n128_min = rep(NA, 11), n128_max = rep(NA, 11), n128_var = rep(NA, 11),
                       n256_min = rep(NA, 11), n256_max = rep(NA, 11), n256_var = rep(NA, 11),
                       n512_min = rep(NA, 11), n512_max = rep(NA, 11), n512_var = rep(NA, 11), 
                       n1024_min = rep(NA, 11), n1024_max = rep(NA, 11), n1024_var = rep(NA, 11))

lambda_est_intervals <- c(0,.05,.15,.25,.35,.45,.55,.65,.75,.85,.95,1)
rownames(range_df) <- c("0-.05","0.05-.15",".15-.25",".25-.35",".35-.45",".45-.55",".55-.65",".65-.75",".75-.85",".85-.95",".95-1")

for(i in 1:11) {
  kappa_data_32_interval <- kappa_data_32[which(kappa_data_32$lambda.est >= lambda_est_intervals[i] & kappa_data_32$lambda.est < lambda_est_intervals[i+1]),]
  range_df$n32_min[i] <- min(kappa_data_32_interval$lambda.input)
  range_df$n32_max[i] <- max(kappa_data_32_interval$lambda.input)
  range_df$n32_var[i] <- var(kappa_data_32_interval$lambda.input)
}

for(i in 1:11) {
  kappa_data_64_interval <- kappa_data_64[which(kappa_data_64$lambda.est >= lambda_est_intervals[i] & kappa_data_64$lambda.est < lambda_est_intervals[i+1]),]
  range_df$n64_min[i] <- min(kappa_data_64_interval$lambda.input)
  range_df$n64_max[i] <- max(kappa_data_64_interval$lambda.input)
  range_df$n64_var[i] <- var(kappa_data_64_interval$lambda.input)
}

for(i in 1:11) {
  kappa_data_128_interval <- kappa_data_128[which(kappa_data_128$lambda.est >= lambda_est_intervals[i] & kappa_data_128$lambda.est < lambda_est_intervals[i+1]),]
  range_df$n128_min[i] <- min(kappa_data_128_interval$lambda.input)
  range_df$n128_max[i] <- max(kappa_data_128_interval$lambda.input)
  range_df$n128_var[i] <- var(kappa_data_128_interval$lambda.input)
}

for(i in 1:11) {
  kappa_data_256_interval <- kappa_data_256[which(kappa_data_256$lambda.est >= lambda_est_intervals[i] & kappa_data_256$lambda.est < lambda_est_intervals[i+1]),]
  range_df$n256_min[i] <- min(kappa_data_256_interval$lambda.input)
  range_df$n256_max[i] <- max(kappa_data_256_interval$lambda.input)
  range_df$n256_var[i] <- var(kappa_data_256_interval$lambda.input)
}

for(i in 1:11) {
  kappa_data_512_interval <- kappa_data_512[which(kappa_data_512$lambda.est >= lambda_est_intervals[i] & kappa_data_512$lambda.est < lambda_est_intervals[i+1]),]
  range_df$n512_min[i] <- min(kappa_data_512_interval$lambda.input)
  range_df$n512_max[i] <- max(kappa_data_512_interval$lambda.input)
  range_df$n512_var[i] <- var(kappa_data_512_interval$lambda.input)
}

for(i in 1:11) {
  kappa_data_1024_interval <- kappa_data_1024[which(kappa_data_1024$lambda.est >= lambda_est_intervals[i] & kappa_data_1024$lambda.est < lambda_est_intervals[i+1]),]
  range_df$n1024_min[i] <- min(kappa_data_1024_interval$lambda.input)
  range_df$n1024_max[i] <- max(kappa_data_1024_interval$lambda.input)
  range_df$n1024_var[i] <- var(kappa_data_1024_interval$lambda.input)
}

range_df[6,]

write.csv(range_df, "Data_Analyses/Munged_Data/Lambda_est_xranges_var.csv", row.names = F)

# Kappa Estimates Norm

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


write.csv(range_df, "Data_Analyses/Munged_Data/Kappa_norm_ranges.csv", row.names = F)


var_df <- data.frame(lambda_input = lambda_input_options, n32 = rep(NA, 21), n64 = rep(NA, 21), 
                       n128 = rep(NA, 21), n256 = rep(NA, 21), n512 = rep(NA, 21), n1024 = rep(NA, 21))
var_df$n32 <- unlist(lapply(1:21, function(i) var(kappa_data_32[which(kappa_data_32$lambda.input == lambda_input_options[[i]]),5])))
var_df$n64 <- unlist(lapply(1:21, function(i) var(kappa_data_64[which(kappa_data_64$lambda.input == lambda_input_options[[i]]),5])))
var_df$n128 <- unlist(lapply(1:21, function(i) var(kappa_data_128[which(kappa_data_128$lambda.input == lambda_input_options[[i]]),5])))
var_df$n256 <- unlist(lapply(1:21, function(i) var(kappa_data_256[which(kappa_data_256$lambda.input == lambda_input_options[[i]]),5])))
var_df$n512 <- unlist(lapply(1:21, function(i) var(kappa_data_512[which(kappa_data_512$lambda.input == lambda_input_options[[i]]),5])))
var_df$n1024 <- unlist(lapply(1:21, function(i) var(kappa_data_1024[which(kappa_data_1024$lambda.input == lambda_input_options[[i]]),5])))

var_df

write.csv(var_df, "Data_Analyses/Munged_Data/Kappa_norm_vars.csv", row.names = F)

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

write.csv(range_df, "Data_Analyses/Munged_Data/Kappa_Z_norm_ranges.csv", row.names = F)


var_df <- data.frame(lambda_input = lambda_input_options, n32 = rep(NA, 21), n64 = rep(NA, 21), 
                       n128 = rep(NA, 21), n256 = rep(NA, 21), n512 = rep(NA, 21), n1024 = rep(NA, 21))
var_df$n32 <- unlist(lapply(1:21, function(i) var(kappa_data_32[which(kappa_data_32$lambda.input == lambda_input_options[[i]]),6])))
var_df$n64 <- unlist(lapply(1:21, function(i) var(kappa_data_64[which(kappa_data_64$lambda.input == lambda_input_options[[i]]),6])))
var_df$n128 <- unlist(lapply(1:21, function(i) var(kappa_data_128[which(kappa_data_128$lambda.input == lambda_input_options[[i]]),6])))
var_df$n256 <- unlist(lapply(1:21, function(i) var(kappa_data_256[which(kappa_data_256$lambda.input == lambda_input_options[[i]]),6])))
var_df$n512 <- unlist(lapply(1:21, function(i) var(kappa_data_512[which(kappa_data_512$lambda.input == lambda_input_options[[i]]),6])))
var_df$n1024 <- unlist(lapply(1:21, function(i) var(kappa_data_1024[which(kappa_data_1024$lambda.input == lambda_input_options[[i]]),6])))

var_df

write.csv(var_df, "Data_Analyses/Munged_Data/Kappa_Z_norm_vars.csv", row.names = F)




# Literature numbers ####

library(readxl)
review <- read_excel("Literature_Review/1Summary.xlsx", sheet = "Since 2019 Lambda Vals")
colnames(review)
length(unique(review$Reference)) # 182

# Published lambda values
published.lambdas <- review[,c(1,3,4)]
published.lambdas.pruned <- published.lambdas[-which(published.lambdas$`Estimated Lambdas Simple` > 1 | published.lambdas$`Estimated Lambdas Simple` < 0),]
published.lambdas.df <- as.data.frame(published.lambdas.pruned)
nrow(published.lambdas.df) # 1556
1556/182 # average number of reported lambda values per manuscript = 8.549451
colnames(published.lambdas.df) <- c("Refs", "Published_Lambdas", "Taxa")
published.lambdas.df$Refs <- as.factor(published.lambdas.df$Refs)

min(table(published.lambdas.df$Refs)) # 1: lowest number of reported lambdas per paper
max(table(published.lambdas.df$Refs)) # 71: highest number of reported lambdas per paper

length(which(published.lambdas.df$Published_Lambdas < 0.05))/nrow(published.lambdas.df)*100 # 25.32%
length(which(published.lambdas.df$Published_Lambdas > 0.90))/nrow(published.lambdas.df)*100 # 24.94%

length(which(published.lambdas.df$Taxa < 30))/nrow(published.lambdas.df)*100 # 22.36504%
length(which(published.lambdas.df$Taxa < 30)) # 348



review_25to75 <- review[which(review$`Estimated Lambdas Simple` >=.25 & review$`Estimated Lambdas Simple` < .7501),]
colnames(review_25to75) <- c("Reference", "Estimated Lambda", "Estimated Lambdas Simple", "Taxa", "__1", "Interpreted", "Quote", "Extra")

review_25to75[which(review_25to75$Interpreted == "Kinda"),6] <- "Yes"
review_25to75$Interpreted <- as.factor(review_25to75$Interpreted)

length(unique(review_25to75$Reference)) # 108 manuscripts published lambdas between .25 and .75
length(which(review_25to75$Interpreted=="Yes"))/length(review_25to75$Interpreted)*100 # 29.21615% of these estimated lambdas are interpreted in terms of magnitude (without statistical test for comparison)

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

# how many papers interpreted lambda with magnitude
keepthese <- all.papers$Reference[which(all.papers$YN=="Yes")]
interpreted<- review.fullsheet[match(keepthese, review.fullsheet$Reference), c("Reference", "Interpreted lambda value as strong or weak beyond significance?")]
colnames(interpreted) <- c("Ref", "Interp")

# taking out any rows with NA in the interpretation column
interpreted.pruned <- interpreted[-which(is.na(interpreted$Interp)=="TRUE"),]

# making a column to simplify these answers
interpreted.pruned$YNKinda <- interpreted.pruned$Interp
interpreted.pruned$YNKinda[startsWith(interpreted.pruned$Interp, "Yes")] <- "Yes"
interpreted.pruned$YNKinda[startsWith(interpreted.pruned$Interp, "No")] <- "No"
interpreted.pruned$YNKinda[startsWith(interpreted.pruned$Interp, "Kinda")] <- "Yes"

# cleaning up messier categories
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Kinda")] <- "Yes"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "-")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "no")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "NA")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Only interpreted kappa")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Only interpreted significance")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Only kappa interpreted")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Lambda not listed in results")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Lambda and kappa have conflicting results regarding significance")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Only interpreted significant")] <- "No"
interpreted.pruned$YNKinda[which(interpreted.pruned$YNKinda == "Only published in SI")] <- "Yes"

interpreted.pruned$YNKinda <- as.factor(interpreted.pruned$YNKinda)
levels(interpreted.pruned$YNKinda)

length(which(interpreted.pruned$YNKinda == "Yes"))/nrow(interpreted.pruned)*100 # 40.59%%



# comparing between lambdas
interpreted$Compared <- as.factor(interpreted$Compared)
summary(interpreted$Compared) # 7 yeses



# Empirical Example: SVL and SAV ####
library(caper)
library(ape)
library(stringr)

sav_dataset <- read.csv("~/Documents/School/Thesis/Data/Pruned/Chapter2/SAtoV.AllSpecimens.csv")
species_sums_sav <- ddply(sav_dataset, ~Species, summarise, SAtoV = mean(SAtoV), SVL = mean(SVL))
row.names(species_sums_sav) <- species_sums_sav[,1]

bwl_dataset <- read.csv("~/Documents/School/Thesis/Data/Pruned/Chapter1/Morpho/LinearMeasurements.IndividualSpecimens.EstimatedMissing.JuvEx.csv")
bwl_dataset <- data.frame(Species = bwl_dataset$Species, BWLrel = bwl_dataset$BWL/bwl_dataset$SVL)
species_sums_bwl <- ddply(bwl_dataset, ~Species, summarise, BWLrel = mean(BWLrel))
species_sums_bwl <- species_sums_bwl[match(row.names(species_sums_sav), species_sums_bwl$Species),] # pruning to the same species as SAV
row.names(species_sums_bwl) <- species_sums_bwl[,1]

phylo <- read.tree("~/Documents/School/Thesis/Data/Pruned/BB.PrunedToLMData.tre")
phylo$tip.label <- str_replace_all(phylo$tip.label, "_", " ")
missing_taxa <- phylo$tip.label[which(phylo$tip.label %in% row.names(species_sums_bwl) == "FALSE")]
phylo_pruned <- drop.tip(phylo, tip = missing_taxa)

species_sums_bwl <- species_sums_bwl[match(phylo_pruned$tip.label, row.names(species_sums_bwl)),] # reordering data to match phylo
species_sums_sav <- species_sums_sav[match(phylo_pruned$tip.label, row.names(species_sums_sav)),] # reordering data to match phylo

# sav vs bwl
sav_named <- species_sums_sav$SAtoV
names(sav_named) <- species_sums_sav$Species
kappa_sav <- physignal(sav_named, phy = phylo_pruned, iter = 999, print.progress = F)
kappa_sav # kappa = 0.7608; z = 8.01; p = 0.001

bwl_named <- species_sums_bwl$BWLrel
names(bwl_named) <- species_sums_bwl$Species
kappa_bwl <- physignal(bwl_named, phy = phylo_pruned, iter = 999, print.progress = F)
kappa_bwl # kappa = 0.2515; z = 2.3597; p = 0.001

faketrait <- sim.char(phy = phylo_pruned, par = 1, nsim = 1)[,,1]
kappa_fake <- physignal(faketrait, phy = phylo_pruned)
plot(kappa_fake)

source("Data_Analyses/CompareZ.r")
sav_vs_bwl <- compare.Z(kappa_sav, kappa_bwl)
sav_vs_bwl$pairwise.P

# Figure 5 B and C
Nulls <- cbind(kappa_sav$random.K[-1], kappa_bwl$random.K[-1])
colnames(Nulls) <- c("SAV.k.null", "BWL.k.null")
Interval<-abs(qnorm(sav_vs_bwl$sample.r.sd.pk.stand)) # CI are represented without 

png("Figures/Fig5_BC.png", width = 600, height = 300)
par(mfrow=c(1,2), mar = c(5,4,2,2))
 # null distribution of kappas with observed kappas
hist(Nulls[,1], col = alpha("red",0.5), xlim = c(0,.9), border = "red", breaks = 15, yaxs = "i", xlab = "Null and Observed Kappas", main = "")
hist(Nulls[,2], col = alpha("black", 0.5), border = "black", xlim = c(0,.9), breaks = 15, yaxs = "i", add = T)
abline(v = kappa_sav$random.K[1], col = "red", lwd = 3)
abline(v = kappa_bwl$random.K[1], col = "black", lwd = 3)
legend("topright",c("SA:V","BW Rel"),fill=c("red","black"), inset = .03)
box()

plotCI(x = c(1,2), y = sav_vs_bwl$sample.z, xlim = c(.8, 2.2), 
       xlab = "Morphological Traits", yaxs="i", axes = F,
       ylim = c(0,10), pch = 19, cex = 2, ylab = "Effect Sizes",
       uiw = Interval, col = c("red", "black"))
axis(2, at=0:10, labels=0:10)
axis(1, at=c(1,2), labels = c("SA:V", "Relative Body Width"))
box()

dev.off()

# Sim 3 traits on 5 on each tree, look at histograms colored by tree

    # on those trees, simulate bm, trend, ou, and do it again.




# tropical sav vs temperate sav; z scores comparison here too

trop_species <- read.csv("~/Documents/School/Thesis/Data/Pruned/TropicalSpecies.csv")
trop_sav <- sav_named[match(trop_species$x, names(sav_named))]
trop_sav <- trop_sav[complete.cases(trop_sav)]
temp_sav <- sav_named[-na.omit(match(trop_species$x, names(sav_named)))]
temp_sav <- temp_sav[complete.cases(temp_sav)]

missing_taxa <- phylo$tip.label[which(phylo$tip.label %in% names(trop_sav) == "FALSE")]
trop_phylo <- drop.tip(phylo, tip = missing_taxa)

missing_taxa <- phylo$tip.label[which(phylo$tip.label %in% names(temp_sav) == "FALSE")]
temp_phylo <- drop.tip(phylo, tip = missing_taxa)

kappa_trop <- physignal(trop_sav, phy = trop_phylo, iter = 999, print.progress = F)
kappa_trop # kappa = 1.5851; z = 6.5741; p = 0.001

kappa_temp <- physignal(temp_sav, phy = temp_phylo, iter = 999, print.progress = F)
kappa_temp # kappa = 0.5223; z = 6.0836; p = 0.001

source("Data_Analyses/CompareZ.r")
temp_vs_trop <- compare.Z(kappa_trop, kappa_bwl)
temp_vs_trop








