# Joshua Burkhart
# 3/14/2017
# PHPH 618
# Problem Set 1

set.seed(88)

library(magrittr)

setwd("~/SoftwareProjects/PHPH618/src")

binding_data <- read.csv("../data/problem_set_1.csv")

# log plot
plot(log2(binding_data$ligand_nM),binding_data$fraction_bound)

# create model of the form Y = (1)x + (2)x^2 + (3)x^3 + (Intercept)
nonlinear_model <- lm(binding_data$fraction_bound ~ poly(log2(binding_data$ligand_nM),3,raw = TRUE))

# (1) = 0.0281644
# (2) = 0.0170429
# (3) = -0.0010036
# (Intercept) = -0.0561749
# Adjusted R-squared:  0.9719
nonlinear_model %>% summary()

# plot nonlinear model
lines(log2(binding_data$ligand_nM),predict(nonlinear_model))

# half_max = 0.4925
half_max <- max(binding_data$fraction_bound) / 2
half_max_r <- round(half_max,digits = 3)

# 0.4925 = 0.0281644x + 0.0170429x^2 + -0.0010036x^3 + -0.0561749
library(polynom)
pnom <- polynomial(nonlinear_model$coefficients)

pnom

# -5.583139  5.862778 16.701541
solve(pnom,b=half_max)

# unlog Kd
Kd <- 5.862778
Kd_r <- round(Kd,digits = 3)

abline(h=half_max)
abline(v=Kd)

text(x=Kd - 3,y=half_max + .1,labels = paste("(",Kd_r,",",half_max_r,")",sep=""))
title(paste("Half-Max = ",half_max_r,", Kd = ",Kd_r,sep=""))
