# Joshua Burkhart
# 3/14/2017
# 3/21/2017 (revised)
# PHPH 618
# Problem Set 1

set.seed(88)

library(magrittr)

setwd("~/SoftwareProjects/PHPH618/src")

binding_data <- read.csv("../data/problem_set_1.csv")

# log plot
plot(log10(binding_data$ligand_nM),binding_data$fraction_bound)

L  <- binding_data$ligand_nM
Kd.init.guess <- 1 # any number works...
Fb  <- binding_data$fraction_bound

# fitting equation
nonlinear_model2 <- nls(Fb ~ (1 / (Kd/L + 1)), start = c(Kd = Kd.init.guess),model = TRUE)

#Formula: Fb ~ (1/(Kd/L + 1))
#Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
#Kd    67.29       4.47   15.05 4.86e-10 ***
nonlinear_model2 %>% summary()

# plot nonlinear model
lines(log10(binding_data$ligand_nM),predict(nonlinear_model2))

# half_max = 0.4925
half_max <- max(binding_data$fraction_bound) / 2
half_max_r <- round(half_max,digits = 3)

# log10 Kd for plotting
Kd <- nonlinear_model2 %>% summary() %>% coefficients() %>% .[1]
Kd_r <- round(Kd,digits = 3)
log10Kd_r <- round(log10(Kd),3)

abline(h=half_max)
abline(v=log10Kd_r)

text(x=log10(Kd) - .6,y=half_max + .06,labels = paste("(",log10Kd_r,",",half_max_r,")",sep=""))
title(paste("Half-Max = ",half_max_r,", Kd = 10^",log10Kd_r," = ",Kd_r,sep=""))
