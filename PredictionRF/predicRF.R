########################################################
######### Predict R/F ratio from model  ################
########################################################
## Author: Maria Rebolleda-Gomez 
## Contact: mreboll1@uci.edu
## Last modified: 2023-07-23
########################################################

############# Set working space ########################
# Load libraries 
library(data.table)
library(tidyverse)

# Parent directory
# In this case we are using data from different parts of the paper
# therefore we set the paent as the overall containing folder. 
# Change folder as needed. 

parent <- "/Users/mrebolleda/Dropbox/Projects/TemperatureCommunities/TC_final/"
setwd(parent)

#########################################################
# Respirator biomass

# To calculate respirator biomass we can use OD620 of communities in acetate
# these communities are dominated by respirators, but we can weight by relative abundance. 

# Open data OD
od_files <- list.files("PredictionRF/OD_data", full.names = TRUE) 
read_od <- function(x){
    dt <-fread(x)
    dt <- dt[, .(Plate, Well, Abs)]
    trf <- strsplit(x, split="_")[[1]][3] 
    dt[, transfer := substr(trf,2,3) %>% as.numeric]
    dt
}

od_data <- lapply(od_files, read_od) %>% rbindlist

# Open carbon sources
carbons=read.csv("PredictionRF/carbon_map.csv")

#merge data
od_carb <- merge(od_data, carbons)

#Remove bubbles and weird artifacts
od_carb=od_carb[od_carb$Abs<=0.5,]

# Plot OD over time data
pdf("PredictionRF/Plots/od_transfers.pdf", 8, 11)
ggplot(data=od_carb, aes(x=transfer, y=Abs,  color=Plate)) +
    facet_wrap(~Carbon, ncol=6)+
    stat_summary(fun.y=mean,  geom="line") +
    geom_point(alpha=0.5)+
    ggtitle("OD at 620nm")+
    ylim(0,0.5)+
    scale_color_manual(values=c("#3288bd","#fdae61","#f46d43","#d7191c","#abdda4"))+
    xlab(NULL)+
    ylab(NULL)+
    theme_bw() 

dev.off()


# Get median OD water to remove base. 
odwater <- od_carb[Carbon == "Water", .(median_od = median(Abs))]

# Get OD only for acetate. Get an average of the last transfers
mean_od <- od_carb[Carbon == "Acetate", .(mean_ab = mean(Abs)-odwater[[1]]), by=c("Plate", "Well", "Replicate")]


# Open data ESV abundance
ESVdata <- fread("CommmunityAssembly_16s/Processed_Data/ESV_fulldata_long.csv")
Acetate_esv <- ESVdata[Carbon=="Acetate",]

# Calculate relative abundace of respirators. 
resp <- fread("CommmunityAssembly_16s/Data/Fermentation_litrev.csv")
esv_resp_ace <- merge(Acetate_esv, resp, by= "Family")


resp_sum <- esv_resp_ace[, .(sum_f=sum(Abundance)),
  by = c("SampleID", "Temperature", "Carbon", "Fermentative")]

# Make wide format and calculate relative abundance of resp. 
resp_sum_wide <- dcast(resp_sum, 
    SampleID + Temperature + Carbon ~ Fermentative, value.var = "sum_f")

resp_relab <- resp_sum_wide[,.(resp_relab = F/(F+`NA`+T)), 
    by = c("SampleID", "Temperature")]

summary(resp_relab$resp_relab)
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.5383  0.7452  0.9731  0.8758  0.9880  1.0000 


# Get respirator biomass 
mean_od$Temperature <- rep(c(22,12.5, 30,37.7, 42), 4)
resp_relab$Replicate <- rep(c(1,2,3,4),5)

od_resp <- merge(mean_od, resp_relab)

# Get respirators 
od_resp[, biomass_r_final := mean_ab*resp_relab]

# Calculate biomass removing initial and assuming all 70mM Carbon was consumed
od_resp[, biomass_r := (biomass_r_final - (biomass_r_final/125))/70]

# Get mean at each temperature
#resp_data <- od_resp[,.(mean_rB=mean(biomass_r)/70),by=Temperature]

########################################################
# Fementer biomass
ferm_data <- fread("GlucoseEfficiency/Data/fermenter_tpc_and_secretions.csv")

ferm_data_sub <-ferm_data[,.(tmp, temperature, well, CorrectedOD,
    Cglucose_consumption, ace_glu_ratio = (Acetate*2)/Cglucose_consumption), acid_glu_ratio]

# Remove data not matching temperature of origin
ferm_data_sm <- ferm_data_sub[tmp==temperature,]

# Sample parameters to obtain predicted values for each temperature
od_resp[, temperature:=rep(c(12,22,30,37,42),4)]

# Get first D/ferm and then multiply by R 
#ferm_data_sm[, D_f := acid_glu_ratio/(CorrectedOD/Cglucose_consumption)]
ferm_data_sm[, D_f := ace_glu_ratio/(CorrectedOD/Cglucose_consumption)]

# Get all combinations of 4 replicate estimates of R and all estimates of D/ferm
temp <- c(12, 22, 30, 37, 42)
rf_list <- vector(mode="list",5)

for (i in 1:5){
    t <- temp[i]
    x <- ferm_data_sm[temperature==t,.(D_f)][[1]]
    y <- od_resp[temperature==t,.(biomass_r)][[1]]
    r_f <- (rep(x, length(y)))*(rep(y, length(x)))
    temperature <- rep(t, length(y)*length(x))
    rf_list[[i]] <- data.table(r_f,temperature)
}

dt_rf_pred <- rbindlist(rf_list)
lapply(rf_list, summary)

r_f_predic <- dt_rf_pred[,.(rf_pred=mean(r_f)), by=temperature]
r_f_predic$temperature <- as.factor(r_f_predic$temperature)

######################################################
# Predict sugar communities
# Open R/F data
rf_data_sugs <- fread("CommmunityAssembly_16s/Processed_Data/rf_data.csv")
rf_data_sugs[, temperature := as.factor(Temp)]

# Remove fucose and water
no_fuc <- rf_data_sugs[Carbon!="Fuc"&Carbon!="H2O",]

# Merge with predicted datas
merged_sugs <- merge(r_f_predic, no_fuc , by="temperature")
mean_carbon <- merged_sugs[, .(meanrf = mean(rf), sdrf=sd(rf), pred=rf_pred), 
    by=c("Temp", "Carbon", "temperature")]


pdf("PredictionRF/Plots/prediction_rf_mean.pdf",6,4)
ggplot(merged_sugs, aes(x=temperature, y=log(rf)))+
    geom_boxplot(fill="gray70")+
    geom_point(aes(y= log(rf_pred)),size=3.5,fill="gray30" , shape=23)+
    #geom_point(data=merged_sugs[Carbon=="Glu",],shape=21,fill="white", size=1.5)+
    theme_bw()
dev.off()


mean_rf <- merged_sugs[,.(mean_rf=mean(rf), sdrf = sd(rf), rf_pred), by=c("Temp")]
mean_rf_glu <- merged_sugs[Carbon=="Glu" & Temp!=12,.(mean_rf=mean(rf), sdrf = sd(rf), rf_pred), by=c("Temp")]

pdf("PredictionRF/Plots/predicted_observed_trend.pdf",4,3)
ggplot(mean_rf[Temp!=12,], aes(y=log(mean_rf), x=Temp))+
    geom_errorbar(aes(ymin=log(mean_rf)-log(sdrf),ymax=log(mean_rf)+log(sdrf)), width=0.3, color="gray30")+
    geom_point(color="gray30", size=1)+
    geom_line(color="gray30")+
    geom_point(data=mean_rf_glu, aes(y= log(mean_rf)),size=1,color="gray60")+
    geom_errorbar(data=mean_rf_glu,aes(ymin=log(mean_rf)-log(sdrf),ymax=log(mean_rf)+log(sdrf)), width=0.3,color="gray60")+
    geom_line(data=mean_rf_glu, aes(y=log(mean_rf)),color="gray60")+
    geom_point(aes(y= log(rf_pred)),size=2, shape=23)+
    geom_line(aes(y= log(rf_pred)),linetype=2)+
    theme_bw()
dev.off()