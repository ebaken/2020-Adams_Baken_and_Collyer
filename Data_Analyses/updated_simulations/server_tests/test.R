library(geomorph)

a<-3
b<-66
c <- a*b
d<-seq(a,b)

plot(d, (d*2 + rnorm(length(d))), pch = 1)

source("/work/LAS/dcadams-lab/baken/updated_simulations/test_func.R")
pth <- test_func(6,9)
write.csv(pth, "test_pth.csv")