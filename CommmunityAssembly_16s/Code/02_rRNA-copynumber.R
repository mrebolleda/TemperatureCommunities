########################################################
## rRNA abundance-weighted mean copy number  (MCN) #####
########################################################
## Author: Maria Rebolleda-Gomez 
## Contact: mreboll1@uci.edu
## Last modified: 2023-26-02
########################################################

############# Set working space ########################
# Load libraries 
library(data.table)
library(tidyverse)

# Assign path to data (change to the right path to directory)
parent = "/Users/mrebolleda/Dropbox/Projects/TemperatureCommunities/TC_final/TemperatureCommunities/CommmunityAssembly_16s"

# Set working directory
setwd(parent)

############## Open and format ESV table ########################
metadata <- fread("Data/TC1_amplicon_metadata.csv")

ESV_dt <- fread("Processed_Data/rarefied_data.csv") %>% setnames("V1", "SampleID")
taxa <- fread("Processed_Data/taxa_rdp.txt") %>% setnames("V1", "ESV")
# Gives warning message because column 1 does not have name. 
# Assigns V1 to the first column, and our code changes it to ESV


# Remove rare taxa (not in the ESV_dt)
taxa_rare <- taxa[ESV%in%colnames(ESV_dt),]
ESV_meta <- merge(metadata, ESV_dt)

# Make ESV long table
ESV_long <- melt(ESV_meta, id.vars = 1:ncol(metadata),
  measure.vars=(ncol(metadata)+1):ncol(ESV_dt),
  variable.name = "ESV", value.name = "Abundance")

ESV_all <- merge(ESV_long, taxa) %>% setnames(c("Genus", "Family", "Order", "Class", "Phylum"), 
c("genus", "family", "order", "class", "phylum"))


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



############# Open and format rRNA copy number data-base
# Database of rRNA copy number (rrnDB pan-taxa statistics 5.8)
rRNA_copy <- fread("Data/rrnDB-5.8_pantaxa_stats_RDP.tsv", sep="\t")

# Remove unecessary columns and split by rank 
rRNA_copy_list <- rRNA_copy[,.(taxid,rank,name,mean)] %>%
    split(by="rank")


########## Match sequentially taxa ##################
setnames(rRNA_copy_list$genus, "name", "genus")
tmp <- ESV_all[genus %in% rRNA_copy_list$genus$genus,]
ESVgenus <- merge(tmp, rRNA_copy_list$genus, by="genus")

tmp1 <- ESV_all[!(genus %in% rRNA_copy_list$genus$genus),]
setnames(rRNA_copy_list$family, "name", "family")
tmp2 <- tmp1[family %in% rRNA_copy_list$family$family]
ESVfamily <- merge(tmp2, rRNA_copy_list$family, by="family")

tmp2.1 <- tmp1[!(family %in% rRNA_copy_list$family$family),]
setnames(rRNA_copy_list$order, "name", "order")
tmp3 <- tmp2.1[order %in% rRNA_copy_list$order$order]
ESVorder <- merge(tmp3, rRNA_copy_list$order, by="order")

tmp3.1 <- tmp2.1[!(order %in% rRNA_copy_list$order$order),]
setnames(rRNA_copy_list$class, "name", "class")
tmp4 <- tmp3.1[class %in% rRNA_copy_list$class$class]
ESVclass <- merge(tmp4, rRNA_copy_list$order, by="order")

# Check what is still missing
ESVs <- c(tmp$ESV, tmp2$ESV, tmp3$ESV, tmp4$ESV)
left <- ESV_all[!(ESV %in% ESVs),]
hist(left$relab)

left$class %>% unique


ggplot(left, aes(x=CarbonID, y=relab*100, color=Temperature))+
    geom_point()+
    theme_bw()
# Rare taxa and non-classified

# Make sure the sets are independent
nrow(tmp) + nrow(tmp2) + nrow(tmp3) + nrow(tmp4)+ nrow(left) == nrow(ESV_all)


######## Calculate MCN ##################
ESV_all_rrna <- rbind(ESVgenus, ESVfamily, ESVorder, ESVclass)
ESV_all_rrna[, weigh_rrna := relab * mean]

MCN <- ESV_all_rrna[ , .(MCN = sum(weigh_rrna)), by = c("SampleID", "Temperature", 'Carbon', "Replicate")]

MCN$Carbon <- factor(MCN$Carbon, levels = c("Glucose", "Maltose",
  "Lactose", "Melibiose", "Galactose", "Sucrose", "Fructose", "Mannose", "Ribose",
  "L-Arabinose", "Mannitol", "Sorbitol", "Galacticol", "Glycerol", "Fucose",
  "Rhamnose", "Oxoglutarate", "Pyruvate", "Oxaloacetate", "Fumarate", "Succinate",
  "Citrate", "Acetate", "Water"))


pdf("Plots/MCN_alltemps.pdf", 7, 10)
ggplot(MCN, aes(x=Temperature, y=MCN, color = as.factor(Replicate)))+
    facet_wrap(~Carbon)+
    geom_point()+
    geom_line()+
    scale_color_manual(values=c("#264653", "#2a9d8f", "#E9C46A", "#EC6442"))+
    theme_bw()
dev.off()


pdf("Plots/MCN_midtemps.pdf", 7, 10)
ggplot(MCN[Temperature > 13 & Temperature < 40,], aes(x=Temperature, y=MCN, color = as.factor(Replicate)))+
    facet_wrap(~Carbon)+
    geom_point()+
    geom_line()+
    scale_color_manual(values=c("#264653", "#2a9d8f", "#E9C46A", "#EC6442"))+
    theme_bw()
dev.off()


######## Calculate MCN - standarized ##################
ESV_all_rrna[, normal_ab :=  Abundance/ mean]

MCN_norm <- ESV_all_rrna[ , .(MCN = sum(Abundance)/sum(normal_ab)), by = c("SampleID", "Temperature", 'Carbon', "Replicate")]

MCN_norm$Carbon <- factor(MCN$Carbon, levels = c("Glucose", "Maltose",
  "Lactose", "Melibiose", "Galactose", "Sucrose", "Fructose", "Mannose", "Ribose",
  "L-Arabinose", "Mannitol", "Sorbitol", "Galacticol", "Glycerol", "Fucose",
  "Rhamnose", "Oxoglutarate", "Pyruvate", "Oxaloacetate", "Fumarate", "Succinate",
  "Citrate", "Acetate", "Water"))


pdf("Plots/MCNnorm_alltemps.pdf", 7, 10)
ggplot(MCN_norm, aes(x=Temperature, y=MCN, color = as.factor(Replicate)))+
    facet_wrap(~Carbon)+
    geom_point()+
    geom_line()+
    ylab("Mean copy number (MCN)")+
    scale_color_manual(values=c("#264653", "#2a9d8f", "#E9C46A", "#EC6442"))+
    theme_bw()
dev.off()


pdf("Plots/MCNnorm_midtemps.pdf", 7, 10)
ggplot(MCN_norm[Temperature > 13 & Temperature < 40,], aes(x=Temperature, y=MCN, color = as.factor(Replicate)))+
    facet_wrap(~Carbon)+
    geom_point()+
    geom_line()+
    ylab("Mean copy number (MCN)")+
    scale_color_manual(values=c("#264653", "#2a9d8f", "#E9C46A", "#EC6442"))+
    theme_bw()
dev.off()
