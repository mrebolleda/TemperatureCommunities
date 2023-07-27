########################################################
######### ESV community analysis #######################
########################################################
## Author: Maria Rebolleda-Gomez 
## Contact: mreboll1@uci.edu
## Last modified: 2023-05-06
########################################################

############# Set working space ########################
# Load libraries 
library(data.table)
library(tidyverse)
library(vegan)
library(emmeans)


# Assign path to data (change to the right path to directory)
#parent = " directory to this folder /CommmunityAssembly_16s/"


# Set working directory
setwd(parent)

# Set color palettes
cols_sugars <- c("#50514F", '#9E6E40',"#D6B79A",
    "#FF9BC3", '#F7A2A1', '#F25F5C',
    '#F9A061', "#FFE066",'#A0B971',"#59A14F",
    "#3B8886", "#70C1B3","#8BCAE5", "#227396",'#594B73')

cols_fr_long <- c("#FFA62B", "#82C0CC","#16697a","gray40",
    "gray60","#aa856d","#f0dcc4","#e65d5d","gray80","gray90")

cols_tmp <- c(c("#2c7bb6","#abd9e9","#f4da8c","#fdae61","#d7191c"))

################ Open files and format data ###############
metadata <- fread("Data/TC1_amplicon_metadata.csv")

ESV_dt <- fread("Processed_Data/rarefied_data.csv") %>% setnames("V1", "SampleID")
taxa <- fread("Processed_Data/taxa_silva_v132.txt") %>% setnames("V1", "ESV")
# Gives warning message because column 1 does not have name. 
# Assigns V1 to the first column, and our code changes it to ESV


# Remove rare taxa (not in the ESV_dt)
taxa_rare <- taxa[ESV%in%colnames(ESV_dt),]
ESV_meta <- merge(metadata, ESV_dt)

# Make ESV long table
ESV_long <- melt(ESV_meta, id.vars = 1:ncol(metadata),
  measure.vars=(ncol(metadata)+1):ncol(ESV_dt),
  variable.name = "ESV", value.name = "Abundance")

ESV_all <- merge(ESV_long, taxa)

# Save ESV table in long format with all metadata
# fwrite(ESV_all, "Processed_Data/ESV_fulldata_long.csv")

#########################################################
########## Calculate relative abundance #################

# SANITY CHECK: Get total reads per sample, 
# because sample is rarefied it should be = to min_reads 
totals <- ESV_all[, sum(Abundance), by = SampleID]
min_reads <- 8064

# Calculate relative abundance - Divide all abundance values by min_reads 
ESV_all[, relab := Abundance/min_reads,]
ESV_all

# SANITY CHECK: Get sum of relative abundace per sample (it should be = 1)
ESV_all[, sum(relab), by = SampleID]


########## Relative abundance plots ######################
ESV_abundant <- ESV_all[relab>0.01,]
ESV_abundant$Carbon <- factor(ESV_abundant$Carbon, levels = c("Glucose", "Maltose",
  "Lactose", "Melibiose", "Galactose", "Sucrose", "Fructose", "Mannose", "Ribose",
  "L-Arabinose", "Mannitol", "Sorbitol", "Galacticol", "Glycerol", "Fucose",
  "Rhamnose", "Oxoglutarate", "Pyruvate", "Oxaloacetate", "Fumarate", "Succinate",
  "Citrate", "Acetate", "Water"))


# Subset sugars only (these are the only relevant carbon sources for this paper)
sugars <- c("Glu", "Mal", "Lac", "Mel", "Gal", "Suc", "Fru", "Man",
  "Rib", "Ara", "Mol", "Sol","Gol", "Gly", "Fuc", "Rha")

ESV_sugars <- ESV_abundant[ESV_abundant$CarbonID %in% c(sugars),]

pdf(file = "Plots/Relative_abundance_sugars.pdf", width = 11, height = 5) 

ggplot(ESV_sugars, aes(x=Replicate, y=relab, group=interaction(Family,Genus), fill=Family))+
  geom_bar(stat="identity", color="black")+
  facet_grid(Temperature~CarbonID)+
  scale_fill_manual(values = cols_fr_long )+
  theme_bw()

dev.off()

##################################################################
############# Calculate R/F ratio and plot #######################
# Rename taxa, create F and R column. 
# Open fermentation data based in literature 

ferm <- fread("Data/Fermentation_litrev.csv")
esv_sugars_h2o <- ESV_all[ESV_all$CarbonID %in% c(sugars, "H2O"),]
esv_ferm <- merge(esv_sugars_h2o, ferm, by= "Family")


ferm_sum <- esv_ferm[, .(sum_f=sum(Abundance)),
  by = c("SampleID", "Temperature", "Carbon", "Fermentative")]

#NA produces weird results
ferm_sum <- ferm_sum[!is.na(Fermentative),]

# Check that we have accounted for most diversity
check <- ferm_sum[,(sum(sum_f))/min_reads,by=SampleID]
check[V1<0.9,] # Only fucose samples are not well covered. 


rf <- function(x){
  if(!any(x$Fermentative=="F")){
    0
  }
  else if (any(x$Fermentative=="T") & any(x$Fermentative=="F")){
    x$sum_f[x$Fermentative=="F"]/x$sum_f[x$Fermentative=="T"]
  } else {
    Inf
  }
}


rf_ratios <- ferm_sum %>% split(by="SampleID") %>% lapply(rf) %>% unlist
data_rf <- data.table(SampleID = names(rf_ratios), rf = rf_ratios) %>% 
  separate(SampleID, c("Temp","Carbon", "Rep"), sep="_")

data_rf$Temp <- substr(data_rf$Temp,3,4) %>% as.numeric

# fwrite(data_rf,"Processed_Data/rf_data.csv")

ggplot(data_rf, aes(x=Temp, y= rf, colour=Carbon))+
  geom_point()+
  stat_smooth(method="lm", se=F)+
  ylab("R/F")+
  xlab("Temp")+
  #scale_color_manual(values = cols_sugars)+
  theme_classic()

# Fucose is an outlier (check plates) as well as T12
data_rf <- as.data.table(data_rf)
no12 <- data_rf[Temp > 12, ]
large_rf <- c("Fuc", "H2O")
subset <- no12[!Carbon %in% large_rf,]

#log transformation of rf makes it farily normaly distributed
subset[, log_rf := log(rf)]

pdf(file = "Plots/rf_16s.pdf", width = 5, height = 3) 
ggplot(subset, aes(x = Temp, y = log_rf, colour = Carbon))+
  geom_point(size=0.2)+
  stat_smooth(method="lm", se=F, size=0.4)+
  stat_smooth(method="lm", se=F, color="black", linetype=2, size=1)+
  ylab("R/F")+
  xlab("Temp")+
  scale_color_manual(values = cols_sugars)+
  theme_classic()

dev.off()


# There is one infinite after log transformation, 
# we need to remove to fit linear model
subset[is.infinite(subset$log_rf),]
subset <- subset[is.finite(subset$log_rf),]

model <- lm(log_rf ~ Temp * Carbon, data = subset)
summary(model)

#To get overall slope controlling for differences in intercepts.
model2 <- lm(log_rf ~ Temp + Carbon, data = subset)
summary(model2)

# Check assumptions of model
par(mfrow = c(2,2))
plot(model)
# Linear model looks fairly decent

# Get slopes coeficients for each carbon
no12[, log_rf:=log(rf)]
no12_noinf <- no12[is.finite(no12$log_rf),]
dt_slopes <- no12_noinf[, .(slope = coef(summary(lm(log_rf ~ Temp)))[2,1],
                    ci_lower = confint(lm(log_rf ~ Temp))[2,1],
                    ci_upper = confint(lm(log_rf ~ Temp))[2,2]),
                by = Carbon]

# Add a new column with the rank of the slope values
dt_slopes[, slope_rank := rank(slope)]
dt_slopes <- dt_slopes[order(slope_rank)]

# Order the carbon_source column by slope_rank
dt_slopes$Carbon <- factor(dt_slopes$Carbon, levels = dt_slopes$Carbon)


pdf(file = "Plots/rf_slope_16s.pdf", width = 2, height = 1.5)
ggplot(dt_slopes, aes(y = slope, x = Carbon))+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper,  width=0.2))+
  geom_point(size=0.5)+
  xlab("Carbon source")+
  ylab(expression(beta))+
  geom_hline(yintercept = 0, linetype = 2)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Get the R/F difference between 22 and 30
only_growthcs <- no12[Carbon %in% c("Glu", "Gal", "Fru", "Rib")]
only_growthcs[Temp==30, mean(rf)] 
only_growthcs[Temp==22, mean(rf)]

############################################################
############################################################
############ Calculate dispersion and beta div #############
# First we need to obtain only sugars

sugar_labels <- separate(ESV_dt[,1],SampleID, c("Temp","Carbon","Rep"))
sugar_rows <- sugar_labels$Carbon %in% sugars

ESV_table_sugars <- ESV_dt[sugar_rows,]

#Make matrix of community data
M_com <- as.matrix(ESV_table_sugars[,2:ncol(ESV_table_sugars)])
rownames(M_com) <- ESV_table_sugars$SampleID

# Calculate relative abundance
# Sample is rarefied so every sample has min_reads = 8064
rowSums(M_com) 
min_reads = 8064
M_relab <- M_com/min_reads

# Calculate bray_curtis pairwise
dist_bray <- vegdist(M_relab,method="bray",binary=F)
M <- as.matrix(dist_bray)

# Labels for pairwise differences 
meta <- sugar_labels[sugar_rows,]

# Get average distance between replicates 
temps <- meta$Temp %>% unique
cs <- meta$Carbon %>% unique

bray <- data.frame(Carbon=character(),
                 Temp=character(), 
                 Mean_bc=numeric(), 
                 stdev = numeric(),
                 stringsAsFactors=FALSE) 

for (i in 1:length(temps)){
  for (j in 1:length(cs)){
    c <- cs[j]
    t <- temps[i]
    sub <- M[meta$Temp==t & meta$Carbon==c, meta$Temp==t & meta$Carbon==c]
    sub<- sub[lower.tri(sub)]
    temp <- substring(t,3,4) %>% as.numeric
    bray[j+(length(cs)*(i-1)),1] <- c
    bray[j+(length(cs)*(i-1)),2:4] <- c(temp, mean(sub), sd(sub))
  }
}

summary(bray$Mean_bc)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.001054 0.096705 0.370629 0.340503 0.526197 0.847884 
summary(bray$Mean_bc[bray$Temp==12])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3202  0.4600  0.5013  0.5168  0.5403  0.8479 

pdf("Plots/pairwisebeta.pdf")
ggplot(bray, aes(y=Mean_bc, x=Temp))+
  geom_boxplot()+
  theme_classic()
dev.off()

bray.cs <- data.frame(Carbon=character(),
                 Temp=character(), 
                 Mean_bc=numeric(), 
                 stdev = numeric(),
                 stringsAsFactors=FALSE)  


# Get bray-curtis differences across carbon sources at different temperatures
for (i in 1:length(temps)){
  for (j in 1:length(cs)){
    c <- cs[j]
    t <- temps[i]
    sub <- M[meta$Temp==t & meta$Carbon==c, meta$Temp==t & meta$Carbon!=c]
    sub<- sub[lower.tri(sub)]
    temp <- substring(t,3,4) %>% as.numeric
    bray.cs[j+(length(cs)*(i-1)),1] <- c
    bray.cs[j+(length(cs)*(i-1)),2:4] <- c(temp, mean(sub), sd(sub))
  }
}

summary(bray.cs$Mean_bc)

ggplot(bray.cs, aes(y=Mean_bc, x=Temp))+
  geom_boxplot()+
  theme_classic()

######################################################################

