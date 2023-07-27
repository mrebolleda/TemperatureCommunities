
########################################################
## Jackie Folmar and Maria Rebolleda-Gomez 
## Contact: mreboll1@uci.edu
## Last updated: 2023-05-10
########################################################


#######################################################
############# Process growth curves 30C ################

########################################################
# Set parameters and load functions

########################################################
# Set parameters and load functions

# Load libraries
library(lubridate)
library(tidyverse)
library(data.table)

# Set color palette
col_rep <- c("#90a955","#4f772d", "#31572c")

# Define parent directory as reference (change as appropriate) 
# parent <- "path to folder /TC_final/GrowthAnalysis/"

# Set parent directory as working directory
setwd(parent)

# Load functions
source("Code/2023-03-06_ProcessGrowthCurves.r")

### Fructose 
fruc <- shape_gc("GrowthCurves/TC42/TC42_Fruc_Day4Plate-GrowthCurve_R1_Reader2_20210225.csv",
  nskip = 9, c = c(1,seq(4,195,by=2)), r = 1:288)

fruc$well <- substr(fruc$well, 9,11)
fruc <- fruc %>% mutate(carbon_source = "fructose") %>% rename(t = `Time [s]`) 
fruc$t <- as.numeric(fruc$t)


### Galactose 
gal <- shape_gc("GrowthCurves/TC42/TC42_Gal_R1_Growthcurve_Reader2_20210219.csv",
  nskip = 9, c = c(1,seq(4,195,by=2)), r = 1:288)

gal$well <- substr(gal$well, 9,11)
gal <- gal %>% mutate(carbon_source = "galactose") %>% rename(t = `Time [s]`) 
gal$t <- as.numeric(gal$t)


### Glucose
glu <- fread("GrowthCurves/TC42/TC42_Glu_R1_Growthcurve_20210219-21_xinsavedon20210222.txt")

glu <- glu[,.(Well,Abs,`Meas. Time [s]`)]
glu[, carbon_source := "glucose"]

setnames(glu, c("Well", "Abs", "Meas. Time [s]"), c("well", "abs", "t"))

### Ribose
rib <- shape_gc("GrowthCurves/TC42/TC42_Rib_R1_Reader4_GrowthCurve_20210219.txt",
  nskip = 32, c = 2, r = 1:289, t = "Time")

rib <- rib %>% mutate(carbon_source = "ribose") %>% rename(t = Time) 
rib$t <- hms(rib$t) %>% period_to_seconds


# colnames all the same so now need to combine to one big df
big_df <- rbind(glu, gal, fruc, rib)

# converting time different units
big_df <- big_df %>% 
  mutate(tmp = 42) 

# making well name model after A1 not A01
big_df$well <- sub("([A-Z])0(.*)", "\\1\\2", big_df$well)


# fwrite <- fwrite(big_df, "Processed_data/gc_TC42.csv")


# Plot data 
carbons <- big_df$carbon_source %>% unique
list_plots <- vector(mode = "list", length = length(carbons)) # making empty list with # of elements of table plots
#creating empty list for for loop to fill in
for (i in 1:length(carbons)) {
   c <- carbons[i]
  list_plots[[i]] <-  ggplot(data = big_df[carbon_source==c], aes(x=t, y=abs))+
        facet_wrap(~well, nrow=8)+
        geom_point(size=0.2, color= col_rep[1])+
        ggtitle(c)+
        theme_bw()
}

pdf("Plots/TC42_platePlots.pdf", width = 10, height = 8)
for (i in 1:length(carbons)) print(list_plots[[i]])

dev.off()


