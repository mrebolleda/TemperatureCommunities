########################################################
######### Predict R/F ratio from model  ################
########################################################
## Author: Ariel Favier & Maria Rebolleda-Gomez 
## Contact: mreboll1@uci.edu
## Last modified: 2023-11-05
########################################################

#####################################################S
# Load libraries
library(data.table)
library(mgcv)
library(lme4)
library(tidyverse)
library(lubridate)
library(dplyr)

# Set parent directory
setwd("TemperatureCommunities/Pairwise")

# Load colony count data
data <- fread("Coassay_colony_counts_231029.csv")

## Column names ##
# Temperature <- temperature of pair incubation 
# Dilution <- streaked dilution
# Replicate <- experimental replicate (1:3)
# Blue <- blue colonies
# White <- white colonies
# Ratio <- White / Blue
# Ratio values were only calculated when there were colonies of both types on the plate

# We calculate the natural logarithm of the R/F ratio
data$logratio <- log(data$Ratio)

# We summarize data by pair and temperature, 
# leaving out evident outliers (there were orders of magnitude higher colony numbers)
meanRF_Pair_Temp <- data %>% filter(log(Ratio)<=1, log(Ratio)>=-3) %>% group_by(Pair,Temperature) %>% summarise(mean = mean(Ratio,na.rm=TRUE),median = median(Ratio,na.rm=TRUE),sd = sd(Ratio,na.rm=TRUE))
meanRF_Temp <- data %>% filter (log(Ratio)<=1, log(Ratio)>=-3) %>% group_by(Temperature) %>% summarise(mean = mean(Ratio,na.rm=TRUE),median = median(Ratio,na.rm=TRUE),sd = sd(Ratio,na.rm=TRUE))


# Fit a linear model between ln(R/F) and Temperature
model2 <- lm(log(mean) ~ Temperature * Pair, data = meanRF_Pair_Temp)
summary(model2)
coefficients2 <- coef(model2)
intercept2 <- coefficients2[1]
slope2 <- coefficients2[2]
###############################

# RF ratio vs temperature

ggplot(meanRF_Pair_Temp, aes(x = Temperature, y = mean,color=as.factor(Pair))) +
  geom_point(size=3) +   
  geom_line() + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.2) +
  labs(
    x = "Temperature",
    y = "R/F Ratio")+
  scale_color_discrete()+
  theme_classic()

# RF ratio in log scale vs temperature

ggplot(meanRF_Pair_Temp, aes(x = Temperature, y = log(mean),color=as.factor(Pair))) +
  geom_point(size=3) +   
  geom_line() + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_errorbar(aes(ymin = log(mean-sd), ymax = log(mean+sd)), width = 0.2) +
    labs(
    x = "Temperature",
    y = "Log(Ratio)")+
  scale_color_discrete()+
  theme_classic()
#################

# RF ratio vs temperature aggregating all pairs

ggplot(meanRF_Temp, aes(x = Temperature, y = log(mean))) +
  geom_point(size=3) +             # Add points for the scatter plot
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    x = "Temperature",
    y = "log(R/F)")+
  geom_errorbar(aes(ymin = log(mean-sd), ymax = log(mean+sd)), width = 0.2) +
  scale_fill_brewer()+
  theme_classic()

#######################
## Growth rate data for the isolates of each pair 

Pair_gr <- fread("C:\\Users\\Ariel\\OneDrive - personalmicrosoftsoftware.uci.edu\\Desktop\\UCI\\Paper temperatura\\CoAssay\\Pair_gr_coassay.csv")

## Column names ##
# Pair <- Pair ID
# F <- Fermenter ID
# R <- Respirator ID
# F_maxrate and R_maxrate <- Maximum growth rate (rmax) calculated 
# after fitting a thermal performance curve to maximum growth rates (r) for Fermenters and Respirators
# at 23,28 and 30C in M9 glucose (1/h)
# F_meanrate and R_meanrate <- Maximum growth rate (rmax) calculated 
# as the average of maximum growth rates (r) at 23,28 and 30C in M9 glucose
# for Fermenters and Respirators (1/h)
# R_F_maxrate <- R_maxrate - F_maxrate (1/h)
# R_F_meanrate <- R_meanrate - F_meanrate (1/h)

newdata <- left_join(Pair_gr,meanRF_Pair_Temp)

# RF ratio vs temp coloring by the difference in rates between respirators and fermenters

ggplot(newdata, aes(x = Temperature, y = log(mean),color=R_F_meanrate,group=Pair)) +
  geom_point(size=3) +   
  geom_line() + 
  geom_errorbar(aes(ymin = log(mean-sd), ymax = log(mean+sd)), width = 0.2) +
  #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    color = "rR-rF",
    x = "Temperature",
    y = "log(Ratio)")+
  ylim(-3,1)+
  scale_color_continuous(type = "viridis")+
  theme_classic()

# RF ratio vs the difference in mean growth rates between respirators and fermenters coloring by temperature

ggplot(newdata, aes(x = R_F_meanrate, y = log(mean),color= Temperature)) +
  geom_point(size=3) +   
  #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    x = "R-F (1/h)",
    y = "log(Ratio)")+
  ylim(-3,1)+
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  geom_errorbar(aes(ymin = log(mean-sd), ymax = log(mean+sd)), width = 0.02) +# Add a linear regression line
  scale_color_continuous(type = "viridis")+
  theme_classic()

# RF ratio vs the difference in mean growth rates between respirators 
# and fermenters coloring and faceting by temperature
newdata <- newdata %>% filter(!is.na(Temperature))

ggplot(newdata, aes(x = R_F_meanrate, y = log(mean),color= Temperature)) +
  geom_point(size=3) +   
  #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    x = "R-F (1/h)",
    y = "log(Ratio)")+
  ylim(-3,1)+
  facet_wrap(~factor(Temperature), ncol=3)+
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  geom_errorbar(aes(ymin = log(mean-sd), ymax = log(mean+sd)), width = 0.02) +# Add a linear regression line
  scale_color_continuous(type = "viridis")+
  theme_classic()

# RF ratio vs the difference in mean growth rates between respirators 
# and fermenters coloring and faceting by temperature

ggplot(newdata, aes(x = R_F_meanrate, y = mean,color= Temperature)) +
  geom_point(size=3) +   
  #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    x = "R-F (1/h)",
    y = "R/F Ratio")+
  ylim(0,2.5)+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.02) +# Add a linear regression line
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  scale_color_continuous(type = "viridis")+
  theme_classic()

# RF ratio vs mean growth rates between respirators 
# and fermenters coloring and faceting by temperature

ggplot(newdata, aes(x = F_meanrate, y = log(mean),color= Temperature)) +
  geom_point(size=3) +   
  #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    x = "Fermenter growth rate (1/h)",
    y = "log(Ratio)")+
  geom_errorbar(aes(ymin = log(mean-sd), ymax = log(mean+sd)), width = 0.02) +# Add a linear regression line
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  scale_color_continuous(type = "viridis")+
  theme_classic()


ggplot(newdata, aes(x = R_meanrate, y = log(mean),color= Temperature)) +
  geom_point(size=3) +   
  #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    x = "Respirator growth rate (1/h)",
    y = "log(Ratio)")+
  geom_errorbar(aes(ymin = log(mean-sd), ymax = log(mean+sd)), width = 0.02) +# Add a linear regression line
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  scale_color_continuous(type = "viridis")+
  theme_classic()

######################################

# Import fermenter secretion data

fermenter_secretions <- fread("C:\\Users\\Ariel\\OneDrive - personalmicrosoftsoftware.uci.edu\\Desktop\\UCI\\Paper temperatura\\coassay_fermenter_secretions.csv")

# Merge into full dataset
coassay_dataset <- left_join(newdata,fermenter_secretions)
coassay_dataset$Ace_gr_ratio <- coassay_dataset$Ace_secretions/coassay_dataset$F_meanrate

## Column names ##

# Ace_secretions <- Acetate concentration (mM) of the supernatant of a 24hs fermenter monoculture in M9 at each temperature 
# Ace_gr_ratio <- Ratio of secretions over fermenter monoculture growth rate in M9 at each temperature

# Correlation between mean fermenter groth rate and acetate concentration

ggplot(coassay_dataset, aes(x = F_meanrate, y = Ace_secretions,color= Temperature)) +
  geom_point(size=3) +   
  #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    x = "F_rate (1/h)",
    y = "Acetate mM")+
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  scale_color_continuous(type = "viridis")+
  theme_classic()

# We obtain only the fermenter data
Fermenter_dataset <- coassay_dataset %>% dplyr::select(F,F_meanrate,F_maxrate,Ace_secretions,Temperature)
Fermenter_dataset$Ace_maxr_ratio <- Fermenter_dataset$Ace_secretions/Fermenter_dataset$F_maxrate
Fermenter_dataset$Ace_meanr_ratio <- Fermenter_dataset$Ace_secretions/Fermenter_dataset$F_meanrate
# Calculate the mean secretion for each Isolate and store it in a new data frame
Fermenter_dataset2 <- Fermenter_dataset %>%
  group_by(F) %>%
  summarize(mean_Ace_secretions_gr = mean(Ace_secretions/F_meanrate))
# Sort the isolates based on the mean secretions
Fermenter_dataset2 <- Fermenter_dataset2[order(Fermenter_dataset2$mean_Ace_secretions_gr), ]

# Acetate secretions per unit of growth rate for the fermenters

ggplot(Fermenter_dataset, aes(x = F, y = Ace_secretions/F_meanrate,color= Temperature,group=Temperature)) +
  geom_point(size=3) +   
  #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    x = "ID",
    y = "Acetate mM/gr")+
  #geom_line()+
  scale_x_discrete(limits = Fermenter_dataset2$F)+
  scale_color_continuous(type = "viridis")+
  theme_classic()


# Acetate secretions of the fermenters vs growth rate

ggplot(Fermenter_dataset, aes(x = F_meanrate, y = Ace_secretions,color= Temperature,group=Temperature)) +
  geom_point(size=3) +   
  #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    x = "r(1/h)",
    y = "Acetate mM")+
  #geom_line()+
  scale_color_continuous(type = "viridis")+
  theme_classic()

# RF ratio in log scale against the acetate secretions 
# and colored by fermenter growth rate

ggplot(coassay_dataset, aes(x = Ace_secretions, y = log(mean) ,color= F_meanrate)) +
  geom_point(size=3) +   
  #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    x = "Acetate secretions",
    y = "log(Ratio)")+
  #geom_line()+
  ylim(-3,1)+
  geom_errorbar(aes(ymin = log(mean-sd), ymax = log(mean+sd)), width = 0.02) +# Add a linear regression line
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  scale_color_continuous(type = "viridis")+
  theme_classic()

# Mean RF ratio in log scale against the acetate secretions 
# and colored by fermenter growth rate

ggplot(coassay_dataset, aes(x = Ace_secretions, y = log(mean) ,color= F_meanrate)) +
  geom_point(size=3) +   
  #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    x = "Acetate secretions",
    y = "log(Ratio)")+
  #geom_line()+
  facet_wrap(~Temperature)+
  ylim(-3,1)+
  geom_errorbar(aes(ymin = log(mean-sd), ymax = log(mean+sd)), width = 0.02) +# Add a linear regression line
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  scale_color_continuous(type = "viridis")+
  theme_classic()

# RF ratio in log scale against the ratio of secretions per unit of growth rate
# colored by temperature

ggplot(coassay_dataset, aes(x = Ace_gr_ratio, y = log(mean) ,color= Temperature)) +
  geom_point(size=3) +   
  #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    x = "Acetate secretions/r",
    y = "log(Ratio)")+
  #geom_line()+
  ylim(-3,1)+
  geom_errorbar(aes(ymin = log(mean-sd), ymax = log(mean+sd)), width = 0.02) +# Add a linear regression line
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  scale_color_continuous(type = "viridis")+
  theme_classic()

# RF ratio in log scale against the ratio of secretions per unit of growth rate
# colored by fermenter growth rate

ggplot(coassay_dataset, aes(x = Ace_gr_ratio, y = log(mean) ,color= F_meanrate)) +
  geom_point(size=3) +   
  #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    x = "Acetate secretions/r",
    y = "log(Ratio)")+
  #geom_line()+
  ylim(-3,1)+
  geom_errorbar(aes(ymin = log(mean-sd), ymax = log(mean+sd)), width = 0.02) +# Add a linear regression line
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  scale_color_continuous(type = "viridis")+
  theme_classic()

# RF ratio in log scale against the ratio of secretions 
# per unit of growth rate faceted by temperature 
# and colored by fermenter growth rate

ggplot(coassay_dataset, aes(x = Ace_gr_ratio, y = log(mean) ,color= F_meanrate)) +
  geom_point(size=3) +   
  #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  labs(
    x = "Acetate secretions/r",
    y = "log(Ratio)")+
  #geom_line()+
  facet_wrap(~Temperature)+
  ylim(-3,1)+
  geom_errorbar(aes(ymin = log(mean-sd), ymax = log(mean+sd)), width = 0.02) +# Add a linear regression line
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add a linear regression line
  scale_color_continuous(type = "viridis")+
  theme_classic()

# Fit a linear model
model <- lm(log(mean) ~ R_F_meanrate + Ace_secretions + R_F_meanrate*Ace_secretions, data = coassay_dataset)
summary(model)


