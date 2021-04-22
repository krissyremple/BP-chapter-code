## To make Figure 2 Table of chem slopes at high and low tides
## By: Krissy Remple
## Last updated: 22 Apr 2021
########################################

#libraries
library(tidyverse)
library(ggplot2)
library(gt)
library(lme4)
library(lmerTest)
library(EnvStats)
library(lm.beta)

#load data
Hi_Res <- read.csv("data/BP_data_Hi_Res_working.csv")


#Prep data by ensuring our favorite columns are the class we want and the order we want

#Assign order to the levels 
Hi_Res$Site.2 <- factor(Hi_Res$Site.2, levels = c("Spring", "Transition", "Diffuse", "Reef"))


#Make day a factor for model 2
Hi_Res$Day <- as.factor(Hi_Res$Day)



#Set colors for plotting later
Site.cols <- c('#d1495b', '#edae49', '#66a182', '#173F5F')
########################################
#Make a dataset for the chemistry of interest

#Get the columns with chemistry data and remove any samples that are missing data

Chem_param <- Hi_Res[,c(1,3,6,13,15,18,20:23,25:32,35,37,38,40:46)]

Chem_complete <-Chem_param[complete.cases(Chem_param),]

#Make a high and low tide dataset for later analysis

Hi_Tide <- Chem_complete %>%
  filter(Tide == "High")


Low_Tide <- Chem_complete %>%
  filter(Tide== "Low")

################################################################################
######################################## ANALYSIS ##############################


# Model 1: fit a model for elements that track salinity and those that track the second degree polynomial of salinity


#Get the cols to run in model 1
Chems <- colnames(Chem_complete[,7:19])

#Run Model 1: 

Mod1_Out <- data.frame()

for (i in Chems){ 
  
  lm_poly = lm(Chem_complete[,i] ~ poly(Salinity, 2, raw = FALSE), data = Chem_complete) 
  Sum_Poly = summary(lm_poly)
  co = Sum_Poly$coefficients
  Beta_lm = lm.beta(lm_poly)
  
  Mod1_Out = rbind(Mod1_Out, data.frame("Nutrient" = i, 
                                        "Salt CoEf" = co[[2]],
                                        "Salt_poly_2 CoEf" = co[[3]],
                                        "Salt P-Value" = co[2,4],
                                        "Salt_poly_2 P-Value" = co[3,4],
                                        "Model RSE" = Sum_Poly$sigma,
                                        "Model R2" = Sum_Poly$r.squared,
                                        "Model Adj R2" = Sum_Poly$adj.r.squared,
                                        "Model F-statistic" = Sum_Poly$fstatistic[1],
                                        "Std Beta Salinity" = Beta_lm$standardized.coefficients[2],
                                        "Std Beta Salinity_poly_2" = Beta_lm$standardized.coefficients[3]))
  
  
}


# Perform for high and low tide

Mod1_Low <- data.frame()

for (i in Chems){ 
  
  lm_poly = lm(Low_Tide[,i] ~ poly(Salinity, 2, raw = FALSE), data = Low_Tide) 
  Sum_Poly = summary(lm_poly)
  co = Sum_Poly$coefficients
  Beta_lm = lm.beta(lm_poly)
  
  Mod1_Low = rbind(Mod1_Low, data.frame("Nutrient" = i, 
                                        "Salt CoEf" = co[[2]],
                                        "Salt_poly_2 CoEf" = co[[3]],
                                        "Salt P-Value" = co[2,4],
                                        "Salt_poly_2 P-Value" = co[3,4],
                                        "Model RSE" = Sum_Poly$sigma,
                                        "Model R2" = Sum_Poly$r.squared,
                                        "Model Adj R2" = Sum_Poly$adj.r.squared,
                                        "Model F-statistic" = Sum_Poly$fstatistic[1],
                                        "Std Beta Salinity" = Beta_lm$standardized.coefficients[2],
                                        "Std Beta Salinity_poly_2" = Beta_lm$standardized.coefficients[3]))
  
  
}

write.csv(Mod1_Low, file = "Mod1_LowTide.csv", row.names = FALSE)


Mod1_High <- data.frame()

for (i in Chems){ 
  
  lm_poly = lm(Hi_Tide[,i] ~ poly(Salinity, 2, raw = FALSE), data = Hi_Tide) 
  Sum_Poly = summary(lm_poly)
  co = Sum_Poly$coefficients
  Beta_lm = lm.beta(lm_poly)
  
  Mod1_High = rbind(Mod1_High, data.frame("Nutrient" = i, 
                                        "Salt CoEf" = co[[2]],
                                        "Salt_poly_2 CoEf" = co[[3]],
                                        "Salt P-Value" = co[2,4],
                                        "Salt_poly_2 P-Value" = co[3,4],
                                        "Model RSE" = Sum_Poly$sigma,
                                        "Model R2" = Sum_Poly$r.squared,
                                        "Model Adj R2" = Sum_Poly$adj.r.squared,
                                        "Model F-statistic" = Sum_Poly$fstatistic[1],
                                        "Std Beta Salinity" = Beta_lm$standardized.coefficients[2],
                                        "Std Beta Salinity_poly_2" = Beta_lm$standardized.coefficients[3]))
  
  
}


write.csv(Mod1_High, file = "Mod1_HighTide.csv", row.names = FALSE)
