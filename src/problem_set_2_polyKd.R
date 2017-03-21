set.seed(88)

library(magrittr)

setwd("~/SoftwareProjects/PHPH618/src")

binding_data <- read.csv("../data/problem_set_2.csv")

# 1.2   2.5   5.2  10.5  21.1  42.1  73.2 111.3 185.5 256.9 325.7 392.0
binding_data$ligand_uM

# 0.00 0.01 0.02 0.03 0.04 0.15 0.42 0.80 0.97 1.00 1.01 0.99
binding_data$fraction_bound

# log plot
plot(log10(binding_data$ligand_uM),binding_data$fraction_bound)

# create model of the form Y = (1)x + (2)x^2 + (3)x^3 + (4)x^4 + (Intercept)
nonlinear_model <- lm(binding_data$fraction_bound ~ poly(log10(binding_data$ligand_uM),4,raw=TRUE))

# (1) = 1.35028
# (2) = -2.81397
# (3) = 1.90015
# (4) = -0.36699
# (Intercept) = -0.11813
# Adjusted R-squared:  0.9744
nonlinear_model %>% summary()

# plot nonlinear model
lines(log10(binding_data$ligand_uM),predict(nonlinear_model))

# half_max = 
half_max <- max(binding_data$fraction_bound) / 2
half_max_r <- round(half_max,digits = 3)

# 0.505 = -0.1181305 + 1.3502812*x - 2.8139697*x^2 + 1.9001524*x^3 - 0.3669943*x^4 
library(polynom)
pnom <- polynomial(nonlinear_model$coefficients)
pnom

# 0.2012695-0.5215514i 0.2012695+0.5215514i 1.8704424+0.0000000i 2.9046256+0.0000000i
solve(pnom,b=half_max)

# unlog Kd
log10Kd_r <- (1.8704424+0i)
Kd <- 10^log10Kd_r
Kd_r <- round(Kd,digits = 3)

abline(h=half_max)
abline(v=log10Kd_r)

text(x=log10(Kd) + .6,y=half_max + .06,labels = paste("(",log10Kd_r,",",half_max_r,")",sep=""))
title(paste("Half-Max = ",half_max_r,", Kd = 10^",log10Kd_r," = ",Kd_r,sep=""))
