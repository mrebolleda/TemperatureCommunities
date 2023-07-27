########################################################
## Jackie Folmar and Maria Rebolleda-Gomez 
## Contact: mreboll1@uci.edu
## Last updated: 2023-05-02
########################################################


#######################################################
############# Process growth curves 22C ################

########################################################
# Set parameters and load functions

# Load libraries
library(lubridate)
library(tidyverse)
library(data.table)

# Set color palette
col_rep <- c("#90a955","#4f772d", "#31572c")

# Define parent directory as reference (change as is appropriate) 
# parent <- "path to folder/TC_final/GrowthAnalysis/"

# Set parent directory as working directory
setwd(parent)

# Load functions
source("Code/2023-03-06_ProcessGrowthCurves.r")

# Open data and transform to long format
gc_files_tpc <- list.files("GrowthCurves/TC22", pattern = "TC_620", full.names = T) # pattern different 
data_tpc <- lapply(gc_files_tpc, shape_gc_batch, n.skip = 17, cols = 3:98) %>% rbindlist # changed cols to 98 given a 96 well plate

# Change name temperature column
setnames(data_tpc,  "T\xb0 620", "tmp")

# Transform time into hours
data_tpc[, t := hms(Time) %>% period_to_seconds]
data_tpc$plate <- data_tpc$plate - 1

# Remove blank plates (1 and 14) and change plate number
data_tpc <- data_tpc[plate != 1 & plate != 14,]
data_tpc$plate <- data_tpc$plate - 1

data_tpc <- data_tpc %>%
  select(-tmp) %>%
  mutate(tmp = 22)

######################################################
# Merge with plate index data and save - unchanged as of now
index_tpc <- fread("TC620_TC22_carbon_source_key.csv")
data_id <- merge(data_tpc, index_tpc, by = "plate")

# Make new time column in seconds and remove old
data_id[, t:= hms(Time) %>% period_to_seconds]
data_id[, Time:= NULL]

# Remove plate 
data_id[, plate:=NULL]

# Write formated data
# fwrite <- fwrite(data_id, "Processed_data/gc_TC22.csv")


### Plot to check replicates
carbons <- data_id$carbon_source %>% unique

list_plots <- vector(mode = "list", length = length(carbons)) # making empty list with # of elements of table plots
#creating empty list for for loop to fill in
for (i in 1:length(carbons)) {
   c <- carbons[i]
  list_plots[[i]] <-  ggplot(data = data_id[carbon_source==c], aes(x=t, y=abs, color=as.factor(replicate)))+
        facet_wrap(~well, nrow=8)+
        geom_point(size=0.2)+
        ggtitle(c)+
        scale_color_manual(values=col_rep)+
        theme_bw()
}

pdf("Plots/TC22_platePlots.pdf", width = 10, height = 8)
for (i in 1:length(carbons)) print(list_plots[[i]])

dev.off()

