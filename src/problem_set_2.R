# Joshua Burkhart
# 3/14/2017
# PHPH 618
# Problem Set 2

set.seed(88)

library(magrittr)

setwd("~/SoftwareProjects/PHPH618/src")

binding_data <- read.csv("../data/problem_set_2.csv")

# log plot
plot(log2(binding_data$ligand_uM),binding_data$fraction_bound)

# create model of the form Y = (1)x + (2)x^2 + (3)x^3 + (4)x^4 + (Intercept)
nonlinear_model <- lm(binding_data$fraction_bound ~ poly(log2(binding_data$ligand_uM),4,raw=TRUE))

# (1) = 0.4064751
# (2) = -0.2549993
# (3) = 0.0518344
# (4) = -0.0030137
# (Intercept) = -0.1181305
# Adjusted R-squared:  0.9744
nonlinear_model %>% summary()

# plot nonlinear model
lines(log2(binding_data$ligand_uM),predict(nonlinear_model))

# half_max = 
half_max <- max(binding_data$fraction_bound) / 2
half_max_r <- round(half_max,digits = 3)

# 0.505 = -0.1181305 + 0.4064751*x - 0.2549993*x^2 + 0.05183436*x^3 - 0.003013689*x^4 
library(polynom)
pnom <- polynomial(nonlinear_model$coefficients)
pnom

# 0.668603-1.732556i 0.668603+1.732556i 6.213475+0i 9.648957+0i
solve(pnom,b=half_max)

# unlog Kd
Kd <- (6.213475+0i)
Kd_r <- round(Kd,digits = 3)

abline(h=half_max)
abline(v=Kd)

text(x=Kd - 3,y=half_max + .1,labels = paste("(",Kd_r,",",half_max_r,")",sep=""))
title(paste("Half-Max = ",half_max_r,", Kd = ",Kd_r,sep=""))

# Hill Equation (based on https://github.com/dritoshi/Fitting-Hill-equation/blob/master/bin/hill.r)
L  <- log2(binding_data$ligand_uM)
n  <- 2
Fb  <- binding_data$fraction_bound

# initial
n.init <- 1

# fitting Hill equation
nonlinear_hill_model <- nls(Fb ~ L^n / (Kd + L^n), start = c(n = n.init))

# n = 0.52436
nonlinear_hill_model %>% summary()

# K0.5^0.52436 = Kd
# K0.5 = Kd^(1/0.52436) = 32.5805
Kd^(1/0.52436)
