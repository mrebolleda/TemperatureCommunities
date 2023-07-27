########################################################
######### Analysis from Estrela 2021 ###################
########################################################
## Author: Maria Rebolleda-Gomez 
## Contact: mreboll1@uci.edu
## Last modified: 2023-07-22
########################################################

# Load packages: 
library("data.table")
library("tidyverse")

# Set parent diorectory 
parent <- "/Users/mrebolleda/Dropbox/Projects/TemperatureCommunities/TC_final/GR_OA_Estrela2021/"

# setwd
setwd(parent)

########################################################
# Open data for organic acids (oa) and growth rate (gr)
oa <- fread("Data/oa.csv")
gr <- fread("Data/gr.csv")

########################################################
# Get acetate and lactate data and get their sum at 16hrs
ace_lac <- oa[(variable=="acetate (mM)"|variable=="lactate (mM)") & time_hours==28,]

oa_new <- dcast(ace_lac, Family + SangerID  ~ variable, value.var = "value")
oa_new[, oa := `acetate (mM)`+ `lactate (mM)`]

# Subset growth rate data to get growth in glucose and merge with org. acids. 
gr_glu <- gr[cs=="glucose",]
gr_ace <- merge(gr_glu, oa_new)

pdf("Plots/oa_gr_estrela.pdf", 6, 4)
ggplot(gr_ace, aes(x=gr_max, y=oa, color=family))+
    geom_point()+
    stat_smooth(method="lm", se=F)+
    ylab("Organic acids at 28hrs [mM]")+
    xlab("Growth rate (/hrs)")+
    scale_color_manual(values=c("black","gray"))+
    theme_bw()
dev.off()

mod_En <- lm(oa~gr_max, data=gr_ace[family=="Enterobacteriaceae",])
summary(mod_En)

mod_P <- lm(oa~gr_max, data=gr_ace[family=="Pseudomonadaceae",])
summary(mod_P)
