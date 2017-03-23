# Joshua Burkhart
# 3/14/2017
# PHPH 618
# Problem Set 2

set.seed(88)

library(magrittr)

setwd("~/SoftwareProjects/PHPH618/src")

binding_data <- read.csv("../data/problem_set_2.csv")

# log plot
plot(log10(binding_data$ligand_uM),binding_data$fraction_bound)

# 1.2   2.5   5.2  10.5  21.1  42.1  73.2 111.3 185.5 256.9 325.7 392.0
binding_data$ligand_uM

# 0.00 0.01 0.02 0.03 0.04 0.15 0.42 0.80 0.97 1.00 1.01 0.99
binding_data$fraction_bound

L  <- binding_data$ligand_uM
Kd.init.guess <- 60
n.init.guess <- 1
Fb  <- binding_data$fraction_bound

# fitting equation
nonlinear_model2 <- nls(Fb ~ (1 / ((Kd/L)^n + 1)), start = c(Kd = Kd.init.guess,n = n.init.guess),model = TRUE,trace=TRUE)

#Formula: Fb ~ (1/((Kd/L)^n + 1))
#Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
#Kd  77.1000     1.7494   44.07 8.69e-13 ***
#n    3.4398     0.2516   13.67 8.50e-08 ***
nonlinear_model2 %>% summary()

# plot nonlinear model
lines(log10(binding_data$ligand_uM),predict(nonlinear_model2))

# half_max = 0.505
half_max <- max(binding_data$fraction_bound) / 2
half_max_r <- round(half_max,digits = 3)

# log10 Kd for plotting
Kd <- nonlinear_model2 %>% summary() %>% coefficients() %>% .[1]
Kd_r <- round(Kd,digits = 3)
log10Kd_r <- round(log10(Kd),3)

n <- nonlinear_model2 %>% summary() %>% coefficients() %>% .[2]

abline(h=half_max)
abline(v=log10Kd_r)

text(x=log10(Kd) - .6,y=half_max + .06,labels = paste("(",log10Kd_r,",",half_max_r,")",sep=""))
title(paste("Half-Max = ",half_max_r,", Kd = 10^",log10Kd_r," = ",Kd_r,sep=""))

# Kd = K0.5^n = K0.5^3.4398 = 77.1000
# K0.5 = Kd^(1/n) = 77.1000^(1/3.4398) = 3.536644
Kd^(1/n)
