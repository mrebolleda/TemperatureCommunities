## " R script for Figure 6"
## Author: Ariel Favier ("afavier@uci.edu")
#  Effect of temperature on the carbon use of fermenters after glucose depletion"

#------------------------------------------------------------------------------

# Load packages
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(IDPmisc)
library(RColorBrewer)
library(nlme)
library(cowplot)
library(data.table)
library(scico)


# In this case we are using data from different parts of the paper
# therefore we set the paent as the overall containing folder. 
# Change folder as needed. 

parent <- "/Users/mrebolleda/Dropbox/Projects/TemperatureCommunities/TC_final/"
setwd(parent)

# Read data table
r_data <- read.csv("GrowthAnalysis/Processed_data/2023-06-14_growth-parameters.csv")
tpc_data <- read.csv("GrowthAnalysis/Processed_data/2023-07-03_coefs_tpc.csv")
fermenter_flux_data <- read.csv("GlucoseEfficiency/Data/fermenter_flux_data.csv")

# Get growth and thermal performance data only for fermenters in glucose
r_fermenters <- r_data %>% filter(as.factor(carbon_source) == "glucose" & as.factor(Fermenter) == "Fermenter" ) %>% select(!temperature)
r_fermenters <- r_fermenters %>% select(!X)

tpc_data_fermenters <- tpc_data %>% filter(as.factor(carbon) == "glucose" & as.factor(fermenter) == "Fermenter" )


# Join fermenter flux data and growth
fermenter_r_flux <- left_join(fermenter_flux_data,r_fermenters)
fermenter_data <- left_join(fermenter_r_flux,tpc_data_fermenters, by = "well")

# Remove NA values
fermenter_data <- fermenter_data %>% NaRV.omit() 

fermenter_data$subopt <- as.numeric(fermenter_data$tmp) > as.numeric(fermenter_data$topt) 
fermenter_data$high_temp <- fermenter_data$tmp  == 42

# Make figure 5B
# First we have to put values of glucose consumption and secretions in a long format as carbon fluxes
CGlucose_consumption_per_cell_fermenters <- fermenter_data %>% pull(Cglucose_consumption_percell) 
Cacid_secretion_per_cell_fermenters <- fermenter_data %>% pull(Cacid_percell)
Exptemperature_fermenters <- fermenter_data %>% pull(tmp)
Exptemperature_fermenters <- append(Exptemperature_fermenters,Exptemperature_fermenters)
Exptemperature_fermenters <- as.numeric(Exptemperature_fermenters)
mM_C <- append(CGlucose_consumption_per_cell_fermenters,Cacid_secretion_per_cell_fermenters)
Flux<- c(rep("CGlucose_consumption",215),rep("Cacid_secretion",215))
Flux <- as.factor(Flux)
IDseq <- fermenter_data %>% pull(SangerID)
IDseq <- append(IDseq,IDseq)

# Write data 
# write.csv(fermenter_data,"TC_final/GlucoseEfficiency/Data/fermenter_tpc_and_secretions.csv")


# Make different versions of figure 6A

Fig6A <- fermenter_data %>% ggplot(aes(x = mean_r, y = Cacid))+
  geom_point(alpha=0.6, aes(color= factor(tmp)),size = 2 )+
  scale_colour_brewer(palette = "RdYlBu",direction=-1,name ="Experimental Temperature [�C]")+
  geom_smooth(method = "lm", se = FALSE, color="darkgray")+
  stat_smooth(method="lm", se = FALSE, aes(linetype=high_temp),color = "darkgray")+
  scale_linetype_manual(values=c("dashed", "dotted"))+
  #ylim(0,7e-11)+
  xlim(0,0.4)+
  xlab(expression(paste("r [hr"^"-1","]")))+
  ylab("OA secretion [mM-C]")+
  theme(axis.title=element_text(size=10),legend.position= c(0.5,0.88),legend.title = element_text(size=9),legend.direction = "horizontal",legend.background = element_rect(fill = "transparent", colour = NA),panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),legend.box.background = element_rect(fill = "transparent", colour = NA),panel.backgroun = element_blank(),axis.line = element_line())+
  guides(linetype="none")+
  guides(colour = guide_legend(title.position = "top"))

# Using total OA secretions and coloring by above thermal optimum and fitting linear models

Fig6A_2 <-fermenter_data %>% filter(!is.na(subopt)) %>% ggplot(aes(x = mean_r, y = Cacid))+
  geom_point(alpha=0.6, aes(color= as.factor(subopt)),size = 2 )+
  scale_color_grey(labels = c("Below topt", "Above topt"),start = 0.1, end = 0.6 )+
  #scale_colour_brewer(palette = "RdYlBu",direction=-1,name ="Experimental Temperature [�C]")+
  geom_smooth(method = "lm", se = FALSE, color="darkgray")+
  stat_smooth(method="lm", se = FALSE, aes(linetype=subopt),color = "darkgray")+
  scale_linetype_manual(values=c("dashed", "dotted"))+
  #ylim(0,7e-11)+
  #xlim(0,0.7)+
  xlab(expression(paste("r [hr"^"-1","]")))+
  ylab("OA secretion [mM-C]")+
  theme(axis.title=element_text(size=10),legend.position= c(0.7,0.88),legend.title = element_text(size=9),legend.direction = "horizontal",legend.background = element_rect(fill = "transparent", colour = NA),panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),legend.box.background = element_rect(fill = "transparent", colour = NA),panel.backgroun = element_blank(),axis.line = element_line())+
  guides(linetype="none")+
  guides(colour = guide_legend(title.position = "top",title=element_blank()))

# Using OA secretions per cell and coloring by above thermal optimum and fitting linear models

Fig6A_3 <- fermenter_data %>% filter(!is.na(subopt)) %>% ggplot(aes(x = mean_r, y = Cacid_percell))+
  geom_point(alpha=0.6, aes(color= as.factor(subopt)),size = 2 )+
  scale_color_grey(labels = c("Below topt", "Above topt"),start = 0.1, end = 0.6 )+
  #scale_colour_brewer(palette = "RdYlBu",direction=-1,name ="Experimental Temperature [�C]")+
  geom_smooth(method = "lm", se = FALSE, color="darkgray")+
  stat_smooth(method="lm", se = FALSE, aes(linetype=subopt),color = "darkgray")+
  scale_linetype_manual(values=c("dashed", "dotted"))+
  ylim(0,7e-11)+
  #xlim(0,0.7)+
  xlab(expression(paste("r [hr"^"-1","]")))+
  ylab("OA secretion per cell [mM-C]")+
  theme(axis.title=element_text(size=10),legend.position= c(0.7,0.88),legend.title = element_text(size=9),legend.direction = "horizontal",legend.background = element_rect(fill = "transparent", colour = NA),panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),legend.box.background = element_rect(fill = "transparent", colour = NA),panel.backgroun = element_blank(),axis.line = element_line())+
  guides(linetype="none")+
  guides(colour = guide_legend(title.position = "top",title=element_blank()))

# We add a column with the difference between experimental temperature and topt

fermenter_data$diff_topt <- fermenter_data$tmp - fermenter_data$topt

Fig6A_4 <- fermenter_data %>% filter(!is.na(subopt)) %>% ggplot(aes(x = mean_r, y = Cacid_percell))+
  geom_point(alpha=0.9, aes(fill= diff_topt),size = 3,shape=21 )+
  scale_fill_scico(palette = 'vik',begin =0 ,end = 1,midpoint = 0,aesthetics = "fill") +
  #scale_color_gradient2( low = "red",mid = "white",    high = "blue",    midpoint = 0,    space = "Lab",   na.value = "grey50")+
  #scale_colour_brewer(palette = "RdYlBu",direction=-1,name ="Experimental Temperature [�C]")+
  #geom_smooth(method = "lm", se = FALSE, color="darkgray")+
  #stat_smooth(method="lm", se = FALSE, aes(linetype=subopt),color = "darkgray")+
  #scale_linetype_manual(values=c("dashed", "dotted"))+
  ylim(0,7e-11)+
  #xlim(0,0.7)+
  xlab(expression(paste("r [hr"^"-1","]")))+
  ylab("OA secretion per cell [mM-C]")+
  theme(axis.title=element_text(size=10),legend.position= c(0.7,0.88),legend.title = element_text(size=9),legend.direction = "horizontal",legend.background = element_rect(fill = "transparent", colour = NA),panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),legend.box.background = element_rect(fill = "transparent", colour = NA),panel.backgroun = element_blank(),axis.line = element_line())+
  #guides(linetype="none")+
  guides(colour = guide_legend(title.position = "top",title=element_blank()))

## Dataframe for figure 6B -> Carbon fluxes
data_fluxes <- data_frame(IDseq,Flux,mM_C,Exptemperature_fermenters)

Fig6B <- data_fluxes %>% NaRV.omit() %>%
  ggplot(aes(x= as.factor(Exptemperature_fermenters), y= mM_C, color=Flux )) +
  geom_point(outlier.color = NA,position=position_dodge(0),alpha=0.5,show.legend = F)+
  ylim(0,7e-11)+
  xlab("Experimental Temperature [�C]")+
  ylab("Carbon [mM]")+
  #scale_x_continuous(breaks = c(12,22,30,37,42))+
  stat_summary(geom="point",fun="median",shape=20,fill="blue",size=6,position=position_dodge(0))+
  stat_summary(geom="line",fun="median",aes(group = Flux, colour = Flux),lwd=0.1,position=position_dodge(0))+
  theme(axis.title=element_text(size=10),legend.position = c(0.5,0.93) ,legend.title= element_blank(),legend.direction = "horizontal", panel.background = element_blank(),legend.key = element_rect(colour = "transparent", fill = "NA"),legend.background = element_rect(fill = "transparent", colour = NA),legend.box.background = element_rect(fill = "transparent", colour = NA),axis.line = element_line(),panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  scale_color_grey(labels = c("OA secretion", "Glucose consumption"),start = 0.1, end = 0.6 )

# Make figure 6c

Fig6C <- fermenter_data %>% ggplot(aes(x = as.factor(tmp), y = Cacid))+
  geom_point(alpha=0.5)+
  geom_boxplot(show.legend = FALSE,alpha = 0.5)+
  ylab("OA secretion [mM-C]")+
  xlab("Experimental temperature [�C]")+
  theme(axis.title=element_text(size=10),panel.background = element_blank(),axis.line = element_line(),panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))

# Make figure 6D and arrange

Fig6D <- fermenter_data %>% ggplot(aes(x = as.factor(tmp), y = (1-acid_glu_ratio)))+
  geom_point(alpha=0.5)+
  geom_boxplot(show.legend = FALSE,alpha = 0.5)+
  ylab(expression(paste("Glucose Efficiency (1-D"[oa/glu],")")))+
  xlab("Experimental temperature [�C]")+
  theme(axis.title=element_text(size=10),panel.background = element_blank(),axis.line = element_line(),panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))



##################################################
## Comparisons across temperature categories

# We run an ANOVA for glucose consumption per cell
res.aov.glucose_consumption_percell <- aov(Cglucose_consumption_percell ~ factor(tmp), data = fermenter_data )
summary(res.aov.glucose_consumption_percell)
TukeyHSD(res.aov.glucose_consumption_percell, conf.level=.95) 


# We run an ANOVA for organic acid secretion per cell
res.aov.oa_secretion_percell <- aov(Cacid_percell ~ factor(tmp), data = fermenter_data )
summary(res.aov.oa_secretion_percell)
TukeyHSD(res.aov.oa_secretion_percell, conf.level=.95) %>% plot


################
### Linear models

# linear regression for glucose consumption per cell and temperature (temp of origin as random)
lm_glucose_consumption_percell <- lme(Cglucose_consumption_percell~tmp, random = ~1|as.factor(temp_origin), data = fermenter_data) 
summary(lm_glucose_consumption_percell)


# linear regression for oa secretions per cell and temperature
lm_oa_secretion_percell <- lme(Cacid_percell~tmp, random = ~1|as.factor(temp_origin), data = fermenter_data) 
summary(lm_oa_secretion_percell) 

## Test effect of 42 in models
fermenter_data_sin42 <- fermenter_data %>% filter (!tmp == "42")

# linear regression for glucose consumption per cell and temperature excluding 42C
# linear regression for glucose consumption per cell and temperature (temp of origin as random)
lm_glucose_percell_no42 <- lme(Cglucose_consumption_percell~tmp, random = ~1|as.factor(temp_origin), data = fermenter_data_sin42) 
summary(lm_glucose_percell_no42)

lm_oa_secretion_percell_sin42 <- lme(Cacid_percell~tmp, random = ~1|as.factor(temp_origin), data = fermenter_data_sin42) 
summary(lm_oa_secretion_percell_sin42) 


# linear regression for total oa secretions and temperature
lm_oa_secretion_total <- lm(Cacid~tmp, data = fermenter_data) #Create the linear regression
summary(lm_oa_secretion_total) #Review the results


# linear regression for oa secretions per cell and temperature removing extreme outliers
lm_oa_secretion_percell <- lm(Cacid_percell~mean_r, data = fermenter_data_no_outliers) #Create the linear regression
summary(lm_oa_secretion_percell) #Review the results

par(mfrow=c(2,2))
plot(lm_oa_secretion_percell_r)

###################################################
# linear regression for oa secretions and r 
lm_oa_r <- lm(Cacid~mean_r, data = fermenter_data) #Create the linear regression
summary(lm_oa_r) 

par(mfrow=c(2,2))
plot(lm_oa_r)

# linear regression for oa secretions below topt
fermenter_data_below_topt <- fermenter_data_no_outliers %>% filter(subopt == TRUE)
lm_oa_r_below_topt <- lm(Cacid~mean_r,data = fermenter_data_below_topt) #Create the linear regression
summary(lm_oa_r_below_topt)

par(mfrow=c(2,2))
plot(lm_oa_r_below_topt)

# linear regression for oa secretions per cell and r 
lm_oa_cell_r <- lm(Cacid_percell~mean_r, data = fermenter_data) #Create the linear regression
summary(lm_oa_cell_r) #Review the results

par(mfrow=c(2,2))
plot(lm_oa_cell_r)
# Poor fit of linear model because of high oa secretion residuals. 

# linear regression for oa secretions per cell below topt
lm_oa_cell_r_below_topt <- lm(Cacid_percell~mean_r, data = fermenter_data_below_topt) #Create the linear regression
summary(lm_oa_cell_r_below_topt) #Review the results

par(mfrow=c(2,2))
plot(lm_oa_cell_r_below_topt)
# Fit is better. Both give the same message. 


#####################

