# Analyses

treesizes <-  2^(5:10)

ANOVA.data <- read.csv("Data_Analyses/Munged_Data/ANOVA.fulldataset.csv")
ANOVA.data <- as.data.frame(ANOVA.data)

Regression.data <- read.csv("Data_Analyses/Munged_Data/Regression.fulldataset.csv")
Regression.data <- as.data.frame(Regression.data)

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
