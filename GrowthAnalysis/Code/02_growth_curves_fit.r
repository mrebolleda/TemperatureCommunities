########################################################
## Jackie Folmar and Maria Rebolleda-Gomez 
## Contact: mreboll1@uci.edu
## Last updated: 2023-06-14
########################################################

############# Calculate growth rate ################

########################################################

#####################################################
############# Set working space 
# Load libraries
library(data.table)
library(mgcv)
library(lme4)
library(tidyverse)
library(latex2exp)
library(emmeans)

# Parent directory (change accordingly)
# parent <- ("path to folder/TC_final/GrowthAnalysis")

# Set parent as working directory
setwd(parent)

# Read files (files for each temperature were processed separately with "01_ProcessGC_TX.R")
t12 <- fread("Processed_data/gc_TC12.csv")
t22 <- fread("Processed_data/gc_TC22.csv")
t30 <- fread("Processed_data/gc_TC30.csv")
t37 <- fread("Processed_data/gc_TC37.csv")
t42 <- fread("Processed_data/gc_TC42.csv")

# Load growth curve function
source("Code/growth_curve_fit_function.r")

#######################################################
############# Format data and merge ################### 

# T37 and T42 do not have replicates. Add a replicate column just with number 1
t37$replicate <- 1
t42$replicate <- 1

# For T37 and T42 there is too much evaporation, biofilm formation and other problems
# after 30 or so hours so we removed these points from analyses. However at lower temperatures 
# it takes longer for isolates to reach their stationary phase.
t37 <- t37[t37$t <=(30*3600)]
t42 <- t42[t42$t <=(30*3600)]


TC_df <- rbind(t12,t22,t30,t37,t42)

# Change time to hours
TC_df$t <- TC_df$t/3600

# Ribose 42 controls are contaminated and so we removed the whole plate
TC_df_nocont <- TC_df[!(carbon_source=="ribose" & tmp==42),]
TC_df_nocont <- TC_df_nocont

#### Make list to fit each growth curve ####
# Calculate the ln of OD
TC_df_nocont <- TC_df_nocont %>% mutate(lOD = log(abs))

TC_df.list <- split(TC_df_nocont, by = c("tmp", "well", "replicate", "carbon_source"))

names_vec <- names(TC_df.list)

# Sanity check
# List should have (96 wells * 4 carbon * 3 temp * 3 reps)+(96 * 4 carb)+(96 * 3carb) 
# We removed one temperature of ribose and there is only one rep of 37 and 42
n = (96 * 4 * 3 * 3)+(96 * 4)+(96 * 3)
length(TC_df.list)== n


####### Fit growth curves
d.fits_tc <- lapply(TC_df.list, compute.gam)
dt.fits <- do.call(rbind, lapply(d.fits_tc, function(x) x[[1]])) %>% as.data.table

#write.csv(dt.fits, "Processed_data/dt.fits.csv")

#### Merge with plate map 
map96 <- fread("platemap_gc.csv")

# replace A01 pattern with A1 pattern so the join works correctly
map96$well <- sub("([A-Z])0(.*)", "\\1\\2", map96$well)
dt.fits.map <- merge(dt.fits, map96)

no_empty <- dt.fits.map[!is.na(Fermenter),] #remove empty wells

cols_keep <- c("well", "t", "tmp", "carbon_source", "carbon", "temperature", "Fermenter")
no_empty_mean <- no_empty[, .(mean_fit = mean(fit), mean_deriv_fit = mean(deriv.fit)), by = cols_keep]


###### Plot growth curve data
pdf("Plots/gc_all_bytemp_deriv.pdf", width = 7, height = 8)

ggplot(no_empty_mean[no_empty_mean$t>2 & no_empty_mean$t<40,], aes(x=t, y=mean_deriv_fit, color=as.factor(temperature)))+
    facet_wrap(tmp~Fermenter, ncol=2)+
    geom_point(alpha=0.1, size=0.2)+
    stat_smooth()+
    scale_color_manual(values=c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"))+
    theme_bw()

dev.off()

# Plot from supplementary material S4
pdf("Plots/gc_all_bytemp.pdf", width = 7, height = 8)

ggplot(no_empty_mean[no_empty_mean$t<30,], aes(x=t, y=mean_fit, color=as.factor(temperature)))+
    facet_wrap(tmp~Fermenter, ncol=2)+
    geom_point(alpha=0.1, size=0.2)+
    stat_smooth()+
    geom_vline(xintercept=24, linetype=2)+
    scale_color_manual(values=c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"))+
    theme_bw()

dev.off()


###### Get the parameters table. 
dt.prm <- do.call(rbind, lapply(d.fits_tc, function(x) x[[2]])) %>% as.data.table
write.csv(dt.prm, "Processed_data/dt.parms.csv")


###### Merge with plate map 
dt.prm.map <- merge(dt.prm, map96)

no_empty_prm <- dt.prm.map[!is.na(Fermenter),] #remove empty wells

cols_keep_prm <- c("well", "tmp", "carbon_source", "carbon", "temperature", "Fermenter")
no_empty_prm <- no_empty_prm[, .(mean_r = mean(r), mean_lag = mean(lag)), by = cols_keep_prm]


#Remove E.coli and Pseudomonas
no_out <- no_empty_prm[!is.na(no_empty_prm$temperature),]
id <- paste(no_empty_prm$well, no_empty_prm$tmp, no_empty_prm$carbon_source)


# write.csv(no_out, "Processed_data/2023-06-14_growth-parameters.csv")

pdf("Plots/growthrate_all_fit.pdf", width = 7, height = 8)

ggplot(no_out, aes(x=tmp, y=mean_r, color=Fermenter))+
    facet_wrap(temperature ~ carbon_source, ncol=4)+
    geom_point(shape=21, alpha=0.6)+
    stat_smooth()+
    scale_color_manual(values=c("#16697a","#aa856d"))+
    theme_bw()

dev.off()

###### Get difference in growth fermenters after accounting for carbon source and 
# temperature of origin. 

dif <- function(x, pair){
    mod <- lm(mean_r ~ Fermenter + carbon_source + as.factor(temperature), data=x)
    emeans <- emmeans(mod, "Fermenter")

    if (pair){
        dif <- pairs(emeans)
        dt <- confint(dif) %>% as.data.table
    } else {
        dt <- emeans %>% as.data.table
    }

    dt$tmp <- x$tmp %>% unique
    dt
}


prm_list <- split(no_out, by=c("tmp"))

mod <- lm(mean_r ~ Fermenter + carbon_source + as.factor(temperature), data=prm_list[[4]])
difs <- lapply(prm_list, dif, pair=T) %>% rbindlist
means <- lapply(prm_list, dif, pair=F) %>% rbindlist

pdf("Plots/difs_gr_rf.pdf", 3, 2)
ggplot(difs, aes(x=tmp, y=estimate))+
    geom_point()+
    geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL), width=0.5)+
    #geom_line()+
    ylab(TeX("$\\mu (T)_{Ferm}-\\mu (T)_{Resp}$"))+
    xlab("Temperature (T) [°C]")+ 
    geom_hline(yintercept=0, linetype=2)+
    theme_bw()
dev.off()


pdf("Plots/gr_rf.pdf", 4, 2)
ggplot(means, aes(x=tmp, y=emmean, color=Fermenter))+
    geom_point()+
    geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL), width=0.5)+
    geom_line()+
    ylab(paste("estimated marginal mean",TeX("$\\mu (T)$")))+
    xlab("Temperature (T) [°C]")+ 
    scale_color_manual(values=c("black","gray40"))+
    theme_bw()
dev.off()