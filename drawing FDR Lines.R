# As per "Uses and misuses of the fudge factor in quantitative discovery proteomics"
# https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/pmic.201600132

# Install & Load packages ----
if (!require("cp4p")) {
  install.packages("cp4p")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("multtest")
  BiocManager::install("qvalue")
  BiocManager::install("siggenes")
  BiocManager::install("limma")
}

library(cp4p) # For data
library(siggenes) # For fudge2()

# Define custom functions ----
source("functions.R")

# Get data ----
tabl <- LFQRatio2

# Columns A.R1, A.R2 and A.R3 correspond to the (numeric) abundance values of proteins in the three replicates of condition A.
# Columns B.R1, B.R2 and B.R3 correspond to the (numeric) abundance values of proteins in the three replicates of condition B.
# Column Welch.test.pval contains the p-values of the Welch t-test between condition A and condition B computed with the Perseus software.

# Calculate the Fold-Change of the averages of each group ----
n1 = 3
n2 = 3
m1 = apply(tabl[, 1:3], 1, mean)
m2 = apply(tabl[, 4:6], 1, mean)
FC = m1 - m2

# Get standard deviations for t-test
sd1 = apply(tabl[, 1:3], 1, sd)
sd2 = apply(tabl[, 4:6], 1, sd)
sdd = sqrt(((sd1^2)+(sd2^2))/(n1+n2))

# Compute the fudge factor with the method of Tusher, V., Tibshirani, R., and Chu, G. (2001)
ff0 = fudge2(FC, sdd)$s.zero

# Get degrees of freedom for significance testing
df = (n1 + n2) - 2

# Get t.stats and p.values 
t.stat <- FC/sdd
p.value <- 2 * (1 - pt(abs(t.stat), df = df))

####################################################
# Prepare displays
abs.inf = -3.5
abs.sup = 3.5

#################################
# Figure a
#################################
# Display the volcano plot on the interval [abs.inf,abs.sup]
plot(FC,-log10(p.value),xlim=c(abs.inf,abs.sup),col=4);

# Display the human proteins which should be identified as differentially abundant
points(FC[(tabl$Organism=="human")],-log10(p.value)[(tabl$Organism=="human")],col=2,pch=3);

confidence_level=0.975; ff=0.5;
ht <- -log10(1-confidence_level);
points(c(-ff,ff), c(10,10), type='h', col='green',lwd=2);
lines(c(abs.inf,abs.sup), c(ht,ht), col='green',lwd=2);

#################################
# Figure b
#################################
# Display the volcano plot on the interval [abs.inf,abs.sup]
plot(FC, -log10(p.value), xlim = c(abs.inf, abs.sup), col = 4)

# Display the human proteins which should be identified as differentially abundant
points(FC[(tabl$Organism == "human")], -log10(p.value)[(tabl$Organism == "human")], col = 2, pch = 3)

confidence_level = 0.975
ff = 0.5
smoothcurve = smooth.threshold(seq(abs.inf, abs.sup, by = 0.0001), ta = qt(confidence_level, df = df), s0 = ff, df = df)
lines(smoothcurve, col = 'green', lwd = 2)


## ggplot2 version:


library(tidyverse)
myData <- tibble(FC = FC, p.value.log = -log10(p.value))

myData$sig <- "black"
# myData$sig[myData$FC[tabl$Organism == "human"] & myData$p.value.log[tabl$Organism == "human"]] <- "red"
myData$sig[tabl$Organism == "human"] <- "red"



# Refactored code version 1 ----
confidence_level = 0.975 # 95% CI 
ff = 0.5
smoothcurve = smooth.threshold(x = seq(abs.inf, abs.sup, by = 0.0001), 
                               ta = qt(confidence_level, df = df),
                               s0 = ff,
                               df = df)

FDR_line <- tibble(FC = smoothcurve[,1],
                   p.value.log = smoothcurve[,2])

FDR_line$p.value.log[FDR_line$p.value.log >= 7.5] <- NA

ggplot(myData, aes(x = FC, y = p.value.log, color = sig)) +
  geom_point(alpha = 0.45, shape = 16) +
  geom_line(data = FDR_line, color = "pink") +
  coord_cartesian(xlim = c(abs.inf, abs.sup), ylim = c(0, 7.5), expand = 0) + # Here xlim ZOOMS-IN on the data
  scale_color_identity()  + # alternatively use scale_color_manual() and define sig as a logical vector (TRUE/FALSE)
  scale_x_continuous("Fold Change") + # This (, limits = c(abs.inf, abs.sup)) FILTERS the data
  scale_y_continuous("-log10(p.value)")
  
# Refactored code version 2 ----
# Using a function which returns a tibble - myVersion
confidence_level = 0.975 # 95% CI 
ff = fudge2(FC, sdd)$s.zero # From S0 and not just 0.5

FDR_line <- ggthreshold(x = seq(abs.inf, abs.sup, by = 0.0001), 
                        ta = qt(confidence_level, df = df),
                        s0 = ff,
                        df = df)

ggplot(myData, aes(FC, p.value.log, color = sig)) +
  geom_point(alpha = 0.45, shape = 16) +
  geom_line(aes(x = x, y = y), data = FDR_line, color = "green") +
  coord_cartesian(xlim = c(abs.inf, abs.sup), ylim = c(0, 7.5), expand = 0) +
  scale_color_identity() +
  NULL

# Refactored code version 3 ----
confidence_level = 0.975 # 95% CI 
ff = 0.5

# Using a function which returns a tibble
# FDR_line <- ggthr (x = seq(abs.inf, abs.sup, by = 0.0001), 
#                         ta = qt(confidence_level, df = df),
#                         s0 = ff,
#                         df = df)

# ggplot(myData, aes(FC, p.value.log, color = sig)) +
#   geom_point(alpha = 0.45, shape = 16) +
#   geom_function()
#   geom_line(aes(x = x, y = y), data = FDR_line, color = "green", inherit.aes = FALSE) +
#   coord_cartesian(xlim = c(abs.inf, abs.sup)) +
#   scale_color_identity() +
#   NULL

#################################
# Figure c
#################################
# Display the volcano plot on the interval [abs.inf,abs.sup]
plot(FC,-log10(p.value),xlim=c(abs.inf,abs.sup),col=4);
# Display the human proteins which should be identified as differentially abundant
points(FC[(tabl$Organism=="human")],-log10(p.value)[(tabl$Organism=="human")],col=2,pch=3);

confidence_level=0.975; ff=0.5; ht <- -log10(1-confidence_level)
points(c(-ff,ff), c(10,10), type='h', col='black',lwd=1)
lines(c(abs.inf,abs.sup), c(ht,ht), col='black',lwd=1)

confidence_level=0.9; ff=0.6;
smoothcurve=smooth.threshold(seq(abs.inf,abs.sup,by=0.0001),ta=qt(confidence_level,df=df),s0=ff,df=df);
lines(smoothcurve, col='gray',lwd=1)

confidence_level=0.75; ff=2;
smoothcurve=smooth.threshold(seq(abs.inf,abs.sup,by=0.0001),ta=qt(confidence_level,df=df),s0=ff,df=df);
lines(smoothcurve, col='gray',lwd=1)

confidence_level=0.95; ff=0.6;
smoothcurve=smooth.threshold(seq(abs.inf,abs.sup,by=0.0001),ta=qt(confidence_level,df=df),s0=ff,df=df);
lines(smoothcurve, col='gray',lwd=1)

confidence_level=0.97; ff=0.15;
smoothcurve=smooth.threshold(seq(abs.inf,abs.sup,by=0.0001),ta=qt(confidence_level,df=df),s0=ff,df=df);
lines(smoothcurve, col='green',lwd=2)

#################################
# Figure d
#################################
# Display the volcano plot on the interval [abs.inf,abs.sup]
plot(FC,-log10(p.value),xlim=c(abs.inf,abs.sup),col=4);
# Display the human proteins which should be identified as differentially abundant
points(FC[(tabl$Organism=="human")],-log10(p.value)[(tabl$Organism=="human")],col=2,pch=3);

confidence_level=0.975; ff=ff0;
smoothcurve=smooth.threshold(seq(abs.inf,abs.sup,by=0.0001),ta=qt(confidence_level,df=df),s0=ff,df=df);
lines(smoothcurve, col='green',lwd=2);
