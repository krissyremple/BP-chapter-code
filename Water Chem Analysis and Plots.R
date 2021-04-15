## Analysis of water chemistry
## By: Krissy Remple
## Last updated: 14 Apr 2021
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




################################################################################
####################### CLASSIFICATIONS AND PLOTS ##############################

#Select SGD Parameters
SGD_Chem <- Mod1_Out %>%
  filter(Salt.P.Value <= 0.05 & Salt.CoEf <0 & (abs(Std.Beta.Salinity) > abs(Std.Beta.Salinity_poly_2)))


#Make a table of summary stats 

#Get the SGD chems of interest
Chem_Name <- SGD_Chem$Nutrient
col.nos <- which(colnames(Chem_complete) %in% Chem_Name)
SGD_Chem_DF <- cbind(Chem_complete[,c(1:3,6)],Chem_complete[col.nos] )

SGD_Chem_Stats <- SGD_Chem_DF %>%
  # tidy data
  gather(Parameter, value, -(1:4)) %>%
  group_by(Parameter, Site.2, Salinity) %>%
  dplyr::summarise(
    Mean = mean(value, na.rm=TRUE),
    Min = min(value, na.rm=TRUE),
    Max = max(value, na.rm=TRUE),
    Geo_Mean = geoMean(value, na.rm = TRUE)
  )

  

#Make some plots
se <- function(x) sd(x)/sqrt(length(x))

#Get data in long format so we can facet later
SGD_Long <- SGD_Chem_DF %>%
  gather("Parameter", "value", -c(1:4))

#Convert the parameter column to a factor and then specify the order we want displayed
SGD_Long$Parameter <- as.factor(SGD_Long$Parameter)

SGD_Long$Parameter <- factor(SGD_Long$Parameter, levels = c("N.N.umol.L", "Phosphate.umol.L", "Silicate.umol.L", "DON.SFIA.uM", "CobleA","HIX", "M.C" ))

#Plot each element against salinity, and the lm trend line with 95% confidence level (is the se default)
SGD_Chem_plots <- ggplot(SGD_Long, aes(x=Salinity, y = value))+
  geom_point()+
  stat_smooth(method = "lm", formula = y~x, size = 1, aes(group = 1))+
  facet_wrap(~Parameter, scales = "free", strip.position = "top" , ncol = 4)+
  theme_light()

SGD_Chem_plots



#ggsave(filename = "Chapter_SGD dotplots.png", plot = SGD_Chem_plots, width = 7, height = 5, units = "in")


#Make a box plot to see  each parameter by site


#Convert the parameter column to a factor and then specify the order we want displayed
SGD_Chem_Stats$Parameter <- as.factor(SGD_Chem_Stats$Parameter)

SGD_Chem_Stats$Parameter <- factor(SGD_Chem_Stats$Parameter, levels = c("N.N.umol.L", "Phosphate.umol.L", "Silicate.umol.L", "DON.SFIA.uM", "CobleA","HIX", "M.C" ))

SGD_Box <- ggplot(SGD_Chem_Stats, aes(x = Site.2, y = Mean, fill = Site.2))+
  geom_boxplot()+
  scale_fill_manual(values = Site.cols)+
  facet_wrap(~Parameter, scales = "free", strip.position = "top", ncol = 4)+
  theme_light()

SGD_Box

#ggsave(filename = "Chapter_SGD Boxplot.png", plot = SGD_Box, width = 12, height = 8, units = "in")

#Make these plots for the other chems
#Select chem parameters from model 1 output that were not selected for SGD.
Reef_Chems <- Mod1_Out[!Mod1_Out$Nutrient %in% Chem_Name,]

#Clear these parameters to avoid accidents
rm(col.nos, Chem_Name)

## Get the Reef chems of interest in a dataframe
Chem_Name <- Reef_Chems$Nutrient
col.nos <- which(colnames(Chem_complete) %in% Chem_Name)
Reef_Chem_DF <- cbind(Chem_complete[,c(1:3,6)],Chem_complete[col.nos] )

#Get the summary statistics
Reef_Chem_Stats <- Reef_Chem_DF %>%
  # tidy data
  gather(Parameter, value, -(1:4)) %>%
  group_by(Parameter, Site.2, Salinity) %>%
  summarise(
    Mean = mean(value, na.rm=TRUE),
    Min = min(value, na.rm=TRUE),
    Max = max(value, na.rm=TRUE),
    Geo_Mean = geoMean(value, na.rm = TRUE)
    
  )



## Make long DF for regressions
Reef_Long <- Reef_Chem_DF %>%
  gather("Parameter", "value", -c(1:4))

#For plotting, make parameter a factor and specify the order
Reef_Long$Parameter <- as.factor(Reef_Long$Parameter)

Reef_Long$Parameter <- factor(Reef_Long$Parameter, levels = c("Bacterioplankton", "NH3.umol.L", "Eukaryotes", "CobleT","Chla.ug.L", "DOC..uM."))

#Make the plots
Reef_Chem_plots <- ggplot(Reef_Long, aes(x=Salinity, y = value))+
  geom_point()+
  stat_smooth(method = "lm", formula = y~x, size = 1, aes(group = 1))+
  stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, aes(group = 1), color = "red")+
  facet_wrap(~Parameter, scales = "free", strip.position = "top", ncol = 2)+
  theme_light()
Reef_Chem_plots


#ggsave(filename = "Chapter Reef Chem Dot Plots.png", plot = Reef_Chem_plots, width = 8, height = 10, units = "in")

Reef_Chem_Stats$Parameter <- as.factor(Reef_Chem_Stats$Parameter)

Reef_Chem_Stats$Parameter <- factor(Reef_Chem_Stats$Parameter, levels = c("Bacterioplankton", "NH3.umol.L", "Eukaryotes", "CobleT","Chla.ug.L", "DOC..uM."))

Reef_Box <- ggplot(Reef_Chem_Stats, aes(x=Site.2, y = Mean, fill = Site.2))+
  geom_boxplot()+
  scale_fill_manual(values = Site.cols)+
  theme_light()+
  facet_wrap(~Parameter, scales = "free", strip.position = "top", ncol = 2)
Reef_Box

#ggsave(filename = "Chapter Reef boxplot.png", plot = Reef_Box, width = 8, height = 10, units = "in")
