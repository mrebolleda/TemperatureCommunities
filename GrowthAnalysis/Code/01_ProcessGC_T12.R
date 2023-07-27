########################################################
## Jackie Folmar and Maria Rebolleda-Gomez 
## Contact: mreboll1@uci.edu
## Last updated: 2023-05-02
########################################################


#######################################################
############# Process growth curves 12C ################

########################################################
# Set parameters and load functions

# Load libraries
library(lubridate)
library(tidyverse)
library(data.table)

# Set color palette
col_rep <- c("#90a955","#4f772d", "#31572c")

# Define parent directory as reference (change as appropriate)
# parent <- "Change file path /TC_final/GrowthAnalysis/"

# Set parent directory as working directory
setwd(parent)

# Load functions
source("Code/2023-03-06_ProcessGrowthCurves.r")

# Read in files
gc_files_tpc12 <- list.files("GrowthCurves/TC12", full.names = T)

# separating out the differently formatted data
gc_files_tpc12_galglu <- gc_files_tpc12[1:6]
gc_files_tpc12_frurib <- gc_files_tpc12[7:12]

# Open data from fructose and ribose and format
data_tpc12_fr <- lapply(gc_files_tpc12_frurib, shape_gc_batch, n.skip = 0, cols = 3:98) %>%
    rbindlist %>%
    rename(Time = V1) %>%
    select(Time, well, abs, plate)

# creating vector of well values in order A1-H1, A2-H2,....
well_vec <- paste(rep(LETTERS[1:8], each = 12), rep(1:12, 8), sep= "")

# pulling all the factors and renaming them - this should keep everything in order
levels(data_tpc12_fr$well) <- well_vec

# convert from factor to character
data_tpc12_fr$well <- as.character(data_tpc12_fr$well)

# Open dara from galactose and fructose 
data_tpc12_gg <- lapply(gc_files_tpc12_galglu, shape_gc_batch, n.skip = 17, cols = 3:98) %>%
  rbindlist %>%
  select(Time, well, abs, plate)

# convert from factor to character
data_tpc12_gg$well <- as.character(data_tpc12_gg$well)

# combine the two dfs into one df
data_tpc12 <- rbind(data_tpc12_gg, data_tpc12_fr)
data_tpc12[, tmp := 12]

# getting time in a workable format - seconds
data_tpc12[,t :=hms(Time) %>% period_to_seconds]

# take out time col
data_tpc12 <- data_tpc12[,2:ncol(data_tpc12)]

# add right carbon sources and replicates
data_tpc12$plate <- data_tpc12$plate-1
index_tpc12 <- fread("TC620_TC22_carbon_source_key.csv")
data_id12 <- merge(data_tpc12, index_tpc12, by = "plate")

### Plot to check replicates
carbons <- data_id12$carbon_source %>% unique

list_plots <- vector(mode = "list", length = length(carbons)) # making empty list with # of elements of table plots
#creating empty list for for loop to fill in
for (i in 1:length(carbons)) {
   c <- carbons[i]
  list_plots[[i]] <-  ggplot(data = data_id12[carbon_source==c], aes(x=t, y=abs, color=as.factor(replicate)))+
        facet_wrap(~well, nrow=8)+
        geom_point(size=0.2)+
        ggtitle(c)+
        scale_color_manual(values=col_rep)+
        theme_bw()
}

pdf("Plots/TC12_platePlots.pdf", width = 10, height = 8)
for (i in 1:length(carbons)) print(list_plots[[i]])

dev.off()

# Finalize data format to merge with others
# Change name carbon source and remove plate column
setnames(data_id12, "carbon source", "carbon_source")
data_id12[, plate:=NULL]

# Save data
# fwrite(data_id12, "Processed_data/gc_TC12.csv")
