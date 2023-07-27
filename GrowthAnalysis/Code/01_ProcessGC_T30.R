########################################################
## Jackie Folmar and Maria Rebolleda-Gomez 
## Contact: mreboll1@uci.edu
## Last updated: 2023-05-02
########################################################


#######################################################
############# Process growth curves 30C ################

########################################################
# Set parameters and load functions

# Load libraries
library(lubridate)
library(tidyverse)
library(data.table)

# Set color palette
col_rep <- c("#90a955","#4f772d", "#31572c")

# Define parent directory as reference (change as appropriate) 
# parent <- "path to directory /TC_final/GrowthAnalysis/"

# Set parent directory as working directory
setwd(parent)

# Load functions
source("Code/2023-03-06_ProcessGrowthCurves.r")

# Get list of files
gc_files_tpc30 <- list.files("GrowthCurves/TC30", 
                             pattern = "TC30_Day4_20211119_48hr_Plate", 
                             full.names = T)

#apply function
data_tpc30 <- lapply(gc_files_tpc30, shape_gc_batch, r=1:161 , n.skip = 0, cols = 3:98) %>% rbindlist 

# Check wells 
well_vec <- paste(rep(LETTERS[1:8], each=12), rep(1:12,8), sep= "")

# pulling all the factors and renaming them - this should keep everything in order
levels(data_tpc30$well) <- well_vec
# convert from factor to  character
data_tpc30$well <- as.character(data_tpc30$well)

data_tpc30 <- data_tpc30 %>%
  select(V1, well, abs, plate) %>%
  rename(time = V1)
# getting time in a workable format - seconds
data_tpc30[,t :=hms(time) %>% period_to_seconds]

# remove blanks and change plate numbers
data_tpc30 <- data_tpc30[plate != 1 & plate != 14,]
data_tpc30$plate <- data_tpc30$plate - 1

data_tpc30 <- data_tpc30 %>%
  select(well,t,abs,plate)

# Merge with plate index data and save - unchanged as of now
index_tpc <- fread("TC620_TC22_carbon_source_key.csv")
data_id30 <- merge(data_tpc30, index_tpc, by = "plate")

##############################
# Finalize details to be able to merge with rest of data
# Add temperature column
data_id30[, tmp:=30]

# Remove old time column, and plate column
data_id30[, Time:= NULL]
data_id30[, plate:=NULL]

# fwrite <- fwrite(data_id30, "Processed_data/gc_TC30.csv")


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

pdf("Plots/TC30_platePlots.pdf", width = 10, height = 8)
for (i in 1:length(carbons)) print(list_plots[[i]])

dev.off()