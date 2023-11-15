######################################################
##### Growth curve data pairs ########################
## Authors: Ariel Favier and Maria Rebolleda-Gomez 
## Contact: mreboll1@uci.edu
## Last updated: 2023-11-13
########################################################


#####################################################
############# Set working space 
# Load libraries
library(data.table)
library(mgcv)
library(lme4)
library(tidyverse)


# Parent directory (change accordingly)
parent <- ("/Users/mrebolleda/Dropbox/Projects/TemperatureCommunities/TC_final/TemperatureCommunities/")

# Set parent as working directory
setwd(parent)


# Load growth curve function
source("GrowthAnalysis/Code/growth_curve_fit_function.r")

# Open metadata (plate map)
map96_coassay <- fread("Pairwise/coassay_gc_platemap.csv")

# Open colony counts
counts <- fread("Pairwise/Coassay_colony_counts_231029.csv")

####################################################
######### Get growth curve data ####################
####################################################

######################################
## Read files and format
raw_files <- list.files("Pairwise/GrowthData/", full.names = TRUE)

addtemp <- function (x){
  n <- x %>% str_locate(pattern = "glu") #  finds the number of plate and converts it is variable  
  n <- n[2]+1 # adds one to the plate possition
  z <- fread(x) 
  z[, tmp:=substring(x,n,n+1) %>% as.numeric]
}

# apply function to open and format all files
data <- lapply(raw_files, addtemp) %>% rbindlist

# Change column names
setnames(data, c("V1","V2","V3","V4","V5","V6","V7","V8","V9","tmp"),
     c("Plate","well","Group","Type","Sample","Kindofod","time","abs","NAS","tmp")) 

data <- data[,c("Plate","well","Type","abs","time","tmp")]

# Convert time intervals to hours
data[, t:=time*0.16551724137]

# Convert from well to row and column
data[, col:=sub("[A-Z]", "" , well) %>% as.numeric] 
data[, row:=substr(well, 1,1) %>% as.factor %>% as.numeric]

# Remove empty wells
isolates <- data[Type=="Unknown",]


# Get log of optical density. There are negative values because blank removal.
# Values under a 0.01 OD are very noisy and before exponential growth, remove:

isolates_no0 <- isolates[abs>0.01,]
isolates_no0[, lOD:=log(abs)]


#################################
# Apply growth curve function 

y <- isolates_no0[,c("abs","t","col","row","tmp","lOD")]
# Make list for each 
y_l <- split(y,by = c("col","row","tmp")) 


### Fit growth curves
d.fits_pair <- lapply(y_l, compute.gam)
dt.fits <- do.call(rbind, lapply(d.fits_pair, function(x) x[[1]])) %>% as.data.table

#write.csv(dt.fits, "Pairwise/Output/dt.fits.pairs.csv")

### Get table with parameters
dt.prm <- do.call(rbind, lapply(d.fits_pair, function(x) x[[2]])) %>% as.data.table
#write.csv(dt.prm, "Pairwise/Output/dt.parms.pairs.csv")


#### Merge with plate map 
fits <- merge(dt.fits,map96_coassay)
parms <- merge(dt.prm,map96_coassay)

########################################
###### Plot growth curve data
pdf("Pairwise/Output/gc_all_derivfit.pdf", width = 7, height = 8)

ggplot(fits, aes(x=t,y=deriv.fit, color=type))+
    facet_wrap(Isolate~tmp)+
    geom_point(size=0.1)+
    xlab("time (hr)")+
    ylab("predicted OD")+
    scale_color_manual(values=c("black","gray40"))+
    theme_bw()

dev.off()

pdf("Pairwise/Output/gc_all_fit.pdf", width = 7, height = 8)

ggplot(fits, aes(x=t,y=fit, color=type))+
    facet_wrap(Isolate~tmp)+
    geom_point(size=0.1)+
    xlab("time (hr)")+
    ylab("predicted OD")+
    scale_color_manual(values=c("black","gray40"))+
    theme_bw()

dev.off()


##########################################################
###### Determine slow growing isolate ####################
###########################################################

########## Map pairs to growth rates of each isolate
# Open pair reference table
pairs <- fread("Pairwise/Pairs.tsv")

# Calculate mean growth rate
parms_mean <- parms[,.(mean_gr=mean(r),sd_gr=sd(r)), by=c("tmp", "Isolate", "type")]

# Make wide table with growth rates for each temp in a different column
Ferm <- parms_mean[type=="Fermenter",] %>%
 setnames("Isolate", "F") %>%
 dcast(F~tmp, value.var="mean_gr") %>%
 setnames(c("23","28","30"), c("F23","F28","F30"))

Resp <- parms_mean[type=="Respirator",] %>%
 setnames("Isolate", "R") %>%
 dcast(R~tmp, value.var="mean_gr") %>%
 setnames(c("23","28","30"), c("R23","R28","R30"))

# Merge with pairs
F_gr <- merge(Ferm, pairs[,1:2])
R_gr <- merge(Resp, pairs[,-2])
pairs_gr <- merge(F_gr, R_gr, by="Pair")

#Calculate differences in growth
difs <- pairs_gr[,.(R23_F23=R23-F23,R28_F28=R28-F28,R30_F30=R30-F30,Pair)]

# Is gr difference always negative or positive? 
difs$slow <- NA

for (i in 1:nrow(difs)){
    y = (difs$R23_F23[i] > 0) + (difs$R28_F28[i] > 0) + (difs$R30_F30[i] > 0)
    if (y == 3){difs$slow[i]<- "F"} 
    else if (y == 0){difs$slow[i] <- "R"}
    else {difs$slow[i] <- "Neither"}
}

##################################################
### Obtain relative freq. slow, R/F ratio and plot
##################################################
# Remove NAs
counts_nona <- counts[!is.na(Blue),]
counts_gr <- merge(difs, counts_nona)

# Calculate frequency of slow growing isolate 
counts_gr[,freq_slow := NA] # First ratio assuming R(White) is slow
counts_gr$freq_slow[counts_gr$slow=="F"] <- counts_gr$Blue[counts_gr$slow=="F"]/(counts_gr$Blue[counts_gr$slow=="F"]+counts_gr$White[counts_gr$slow=="F"])
counts_gr$freq_slow[counts_gr$slow=="R"] <- counts_gr$White[counts_gr$slow=="R"]/(counts_gr$Blue[counts_gr$slow=="R"]+counts_gr$White[counts_gr$slow=="R"])


# Summarize data by pair and temperature. 
meanFreq_Pair_Temp <- counts_gr[!is.na(freq_slow),
    .(mean = mean(freq_slow),median = median(freq_slow),sd = sd(freq_slow)), 
    by=c("Pair","Temperature", "slow")]

pdf("Pairwise/Output/freq_slow.pdf", width = 6, height = 4)
ggplot(meanFreq_Pair_Temp, aes(x=Temperature, y=mean, group=as.factor(Pair)))+
  facet_wrap(1~slow)+
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2)+
  ylab("Frequency of slow")+
  xlab("Temperature")+
  geom_line()+
  theme_bw()
dev.off()


# Calculate R/F ratio
counts_gr[, RF_ratio := White/Blue, ]

# Plot distributio of R/F adding interquartile range from Estrela., 2021
pdf("Pairwise/Output/dist_rf.pdf", width = 6, height = 4)
ggplot(counts_gr, aes(x=log(RF_ratio)))+
  facet_wrap(~Temperature)+
  geom_density()+
  geom_vline(xintercept = log(0.17), linetype=2)+
  geom_vline(xintercept = log(0.69), linetype=2)+
  theme_bw()
dev.off()


##########################################
#### Plot growth rates 

# There are two sets o f pairs -> dominance or coexixtence 
# from distribution of colony counts we can see that less than 10 colonies 
# of each type is within a different peak distribution. 

gr_rf_mean <- counts_gr[RF_ratio > 0.02 & RF_ratio < 7 ,
    .(mean_rf = mean(log(RF_ratio)), sd_rf= sd(log(RF_ratio))), 
    by =c("Pair", "slow", "Temperature", "mean_gdif")]

pdf("Pairwise/Output/grdif_rf.pdf", width = 6, height = 4)
ggplot(gr_rf_mean, aes(y=mean_rf, x=mean_gdif))+
  facet_wrap(~Temperature)+
  geom_point()+
  geom_errorbar(aes(ymin=mean_rf-sd_rf, ymax=mean_rf+sd_rf))+
  stat_smooth(method="lm", se=F)+
  xlab("mean(r(R)-r(F)) (/hr)")+
  ylab("log(R/F)")+
  geom_vline(xintercept = 0, linetype = 2)+
  theme_bw()
dev.off()

