# Joshua Burkhart
# 3/14/2017
# PHPH 618
# Problem Set 1

library(magrittr)

setwd("~/SoftwareProjects/PHPH618/src")

binding_data <- read.csv("../data/problem_set_1.csv")

# log plot
plot(log2(binding_data$ligand_nM),binding_data$fraction_bound)

# create model of the form Y = (1)x + (2)x^2 + (3)x^3 + (Intercept)
nonlinear_model <- lm(binding_data$fraction_bound ~ poly(log2(binding_data$ligand_nM),3))

# (1) = 1.45977
# (2) = 0.29031
# (3) = -0.33380
# (Intercept) = 0.40073
# Adjusted R-squared:  0.9719
nonlinear_model %>% summary()

# plot nonlinear model
lines(log2(binding_data$ligand_nM),predict(nonlinear_model))

# half_max = 0.4925
half_max <- max(binding_data$fraction_bound) / 2
half_max_r <- round(half_max,digits = 3)

# 0.4925 = 1.45977x + 0.29031x^2 + -0.33380x^3 + 0.40073
library(polynom)
pnom <- polynomial(c(0.40073,1.45977,0.29031,-0.33380))
pnom

# -1.73781250  0.06215273  2.54537217
solve(pnom,b=half_max)

# unlog Kd
Kd <- 2^2.54537217
Kd_r <- round(Kd,digits = 3)

abline(h=half_max)
abline(v=Kd)

text(x=Kd - 3,y=half_max + .1,labels = paste("(",Kd_r,",",half_max_r,")",sep=""))
title(paste("Half-Max = ",half_max_r,", Kd = ",Kd_r,sep=""))
