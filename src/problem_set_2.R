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

L  <- binding_data$ligand_uM
Kd.init.guess <- 1 # any number works...
Fb  <- binding_data$fraction_bound

# fitting equation
nonlinear_model2 <- nls(Fb ~ (1 / (Kd/L + 1)), start = c(Kd = Kd.init.guess),model = TRUE,trace=TRUE)

#Formula: Fb ~ (1/(Kd/L + 1))
#Parameters:
#   Estimate Std. Error t value Pr(>|t|)   
#Kd    62.82      18.22   3.447  0.00546 **
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

abline(h=half_max)
abline(v=log10Kd_r)

text(x=log10(Kd) - .6,y=half_max + .06,labels = paste("(",log10Kd_r,",",half_max_r,")",sep=""))
title(paste("Half-Max = ",half_max_r,", Kd = 10^",log10Kd_r," = ",Kd_r,sep=""))

# Hill Equation (based on https://github.com/dritoshi/Fitting-Hill-equation/blob/master/bin/hill.r)
L  <- binding_data$ligand_uM
n.init <- 1 # any number works...
Fb  <- binding_data$fraction_bound

# initial


# fitting Hill equation
nonlinear_hill_model <- nls(Fb ~ (L^n / (Kd + L^n)), start = c(n = n.init))

#Formula: Fb ~ (L^n/(Kd + L^n))
#Parameters:
#  Estimate Std. Error t value Pr(>|t|)    
#n  1.03272    0.06766   15.26 9.48e-09 ***
nonlinear_hill_model %>% summary()

n <- nonlinear_hill_model %>% summary() %>% coefficients() %>% .[1]

# K0.5^1.03272 = Kd
# K0.5 = Kd^(1/0.52436) = 32.5805
Kd^(1/n)
