
########################################################
## Maria Rebolleda-Gomez 
## Contact: mreboll1@uci.edu
## Last updated: 2023-07-18
########################################################

###########################################################
# Code to get growth at 30C in two different plate readers#

# Load libraries
library(data.table)
library(mgcv)
library(tidyverse)

# Parent directory (change accordingly)
parent <- ("/Users/mrebolleda/Dropbox/Projects/TemperatureCommunities/TC_final/GrowthAnalysis")

# Set parent as working directory
setwd(parent)

# Load useful functions
source("Code/2023-03-06_ProcessGrowthCurves.r")
source("Code/growth_curve_fit_function.r")

# Read file from the other plate reader - growth in glucose at 30C
t30_plateReader2 <- shape_gc("GrowthCurves/TC30/TC30_GrowthCurve_20220317-19.csv", nskip = 9, 
    c = c(1, seq(4, 194, by = 2)), r = 1:288)

# Read file with previously calculated growth rates.
prms <- fread("Processed_data/2023-06-14_growth-parameters.csv")

# Open metadata
map96 <- fread("platemap_gc.csv")

##########################################
# Format data to evaluate growth.
setnames(t30_plateReader2, "Time [s]", "t")

# Change time to hours
t30_plateReader2$t <- t30_plateReader2$t/3600

# Get well names
t30_plateReader2$well <- substr(t30_plateReader2$well,9,11)
t30_plateReader2[, lOD:=log(abs)]

#################################################
# Get growth parameters
# Make a list for each well to calculate growth

# Change well format of the map and the t30 data
map96$well <- sub("([A-Z])0(.*)", "\\1\\2", map96$well)
t30_plateReader2$well <- sub("([A-Z])0(.*)", "\\1\\2", t30_plateReader2$well)
t30_plateReader2[t<36,]

# Separate into list and apply function
t30.2_list <- split(t30_plateReader2, by = "well")
d.fits_tc <- lapply(t30.2_list, compute.gam)

##### Get fits.
dt.fits <- do.call(rbind, lapply(d.fits_tc, function(x) x[[1]])) %>% as.data.table

###### Get the parameters table.
dt.prm <- do.call(rbind, lapply(d.fits_tc, function(x) x[[2]])) %>% as.data.table

##### Merge with plate map.
dt.prm.map <- merge(dt.prm, map96)
dt.fits.map <- merge(dt.fits, map96)
no_empty.prm <- dt.prm.map[!is.na(Fermenter),] #remove empty wells
no_empty.fits<- dt.fits.map[!is.na(Fermenter),] #remove empty wells


################################################
############## Compare data ###################
dt.fits <- fread("Processed_data/dt.fits.csv")
glu30.fit <- dt.fits[carbon_source=="glucose" & tmp==30,]

glu30.fit.map <- merge(glu30.fit, map96)
glu30_noempty <- glu30.fit.map[!is.na(Fermenter),] #remove empty wells


# Merge fits 
no_empty.fits$replicate <- 4
data.all <- rbind(no_empty.fits[,c("fit", "t", "well", "replicate")],
        glu30_noempty[,c("fit", "t", "well", "replicate")])

pdf("Plots/TC30_comparison.pdf", width = 10, height = 8)
ggplot(data = data.all[data.all$t < 24], aes(x=t, y=fit, color=as.factor(replicate)))+
        facet_wrap(~well, nrow=8)+
        geom_point(size=0.2)+
        scale_color_manual(values=c("gray20", "gray40", "gray60", "#90a955"))+
        theme_bw()

dev.off()

