
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
# parent <- "path to folder /GrowthAnalysis/"

# Set parent directory as working directory
setwd(parent)

# Load functions
source("Code/2023-03-06_ProcessGrowthCurves.r")

### Fructose files ######
# Get list of files
fruc <- shape_gc("GrowthCurves/TC37/TC37_Fruc_R1_GrowthCurve_Reader3_20210210.txt",
  nskip = 9, c = c(1,seq(4,195,by=2)), r = 1:289)

fruc$well <- substr(fruc$well, 9,11)
fruc <- fruc %>% mutate(carbon_source ="fructose") %>% rename(t = `Time [s]`) 
fruc$t <- as.numeric(fruc$t)

 ggplot(data = fruc, aes(x=t, y=abs))+
        facet_wrap(~well, nrow=8)+
        geom_point(size=0.2, color = col_rep[1])+
        theme_bw()

### Glucose files ####
glu <- shape_gc("GrowthCurves/TC37/TC37_growthcurves_202220413-15.txt", nskip = 20, c=2, r=1:289, t="Time")
setnames(glu, "Time", "t")

# convert time to seconds
glu$t <- hms(glu$t) %>% period_to_seconds
#Add plate number
glu[, carbon_source := 'glucose']


###ribose files####
rib <- shape_gc("GrowthCurves/TC37/TC37_Rib_R1_Reader4_GrowthCurve_20210210.txt",
  nskip = 33, c = 2, r = 1:270, t = "Time")

# convert time to seconds
setnames(rib, "Time", "t")
rib$t <- hms(rib$t) %>% period_to_seconds

#Add plate number
rib[, carbon_source := 'ribose']


####galactose files####
gal <- shape_gc("GrowthCurves/TC37/TC37_Gal_R1_growthcurve_reader2_20210210.csv",
  nskip = 9, c = c(1, seq(4,195,by=2)), r = 1:288)

gal$well <- substr(gal$well, 9,11)
gal <- gal %>%
  mutate(carbon_source = "galactose") %>%
  rename(t = `Time [s]`)

gal$t <- as.numeric(gal$t)

## stack results into one big df
big_df37 <- rbind(glu, gal, fruc, rib)

#### Format data.table for merging with other temperatures
# changing well names from something like "A01" to "A1"
big_df37$well <- sub("([A-Z])0(.*)", "\\1\\2", big_df37$well)

# Add temperature
big_df37[, tmp:=37]

carbons <- big_df37$carbon_source %>% unique

# Plot data 
list_plots <- vector(mode = "list", length = length(carbons)) # making empty list with # of elements of table plots
#creating empty list for for loop to fill in
for (i in 1:length(carbons)) {
   c <- carbons[i]
  list_plots[[i]] <-  ggplot(data = big_df37[carbon_source==c], aes(x=t, y=abs))+
        facet_wrap(~well, nrow=8)+
        geom_point(size=0.2, color= col_rep[1])+
        ggtitle(c)+
        theme_bw()
}

pdf("Plots/TC37_platePlots.pdf", width = 10, height = 8)
for (i in 1:length(carbons)) print(list_plots[[i]])

dev.off()

# fwrite(big_df37, "Processed_data/gc_TC37.csv")

