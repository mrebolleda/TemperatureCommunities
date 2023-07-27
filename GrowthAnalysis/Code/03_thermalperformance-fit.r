########################################################
## Maria Rebolleda-Gomez 
## Contact: mreboll1@uci.edu
## Last updated: 2023-05-02
########################################################

########## Fit thermal performance curves #############

#####################################################
############# Set working space
# Load libraries
library(data.table)
#remotes::install_github("padpadpadpad/rTPC")
library(rTPC)
library(nls.multstart)
library(broom)
library(lme4)
library(tidyverse)

# Parent directory (change accordingly)
parent <- ("/Users/mrebolleda/Dropbox/Projects/TemperatureCommunities/TC_final/GrowthAnalysis")

# Set parent as working directory
setwd(parent)

# Open data
dt.parms <- fread("Processed_data/2023-06-14_growth-parameters.csv")


# Create ID to subset each growth curve
dt.parms$ID <- paste(dt.parms$well, dt.parms$carbon_source, dt.parms$temperature, dt.parms$Fermenter, sep="_")
id <- unique(dt.parms$ID)[1]

dt.list <- split(dt.parms, dt.parms$ID)


### Fit gaussian response temperature
gaussian_1987(temp, rmax, topt, a)
fit_list_gaussian <- function(d){
    mod <- 'gaussian_1987'
    start_vals <- get_start_vals(d$tmp, d$mean_r, model_name = mod)
    low_lims <- get_lower_lims(d$tmp, d$mean_r, model_name = mod)
    upper_lims <- get_upper_lims(d$tmp, d$mean_r, model_name = mod) 


    fit <- nls_multstart(mean_r~gaussian_1987(temp = tmp, rmax, topt, a),
                                                     data = d,
                                                     iter = 500,
                                                     start_lower = start_vals - 5,
                                                     start_upper = start_vals + 5,
                                                     lower = low_lims,
                                                     upper = upper_lims,
                                                     supp_errors = 'Y')
    fit
}

# Split list to fit curves
dt.list <- split(dt.parms, dt.parms$ID)


fit_list <- lapply(dt.list, fit_list_gaussian)
new_data <- data.frame(tmp = seq(12, 42, 0.5))
preds <- lapply(fit_list, augment, newdata = new_data) %>% rbindlist(idcol = "ID")

coef_dt <- function(x){
    data.table(rmax = coef(x)[1], topt = coef(x)[2], a = coef(x)[3], rssq = sum(residuals(x)^2))
    }

coefs <- lapply(fit_list, coef_dt) %>% rbindlist(idcol = "ID")

coefs[, c("well", "carbon", "temp_origin", "fermenter") := tstrsplit(ID, "_", fixed=TRUE)]

#write.csv(coefs, "Processed_data/2023-07-03_coefs_tpc.csv")

ggplot(coefs, aes(x = as.factor(temp_origin), y = topt, fill = fermenter))+
    geom_boxplot()+
    theme_bw()

preds[, c("well", "carbon_source", "temperature", "Fermenter") := tstrsplit(ID, "_", fixed=TRUE)]

ggplot(dt.parms, aes(x = tmp, y = mean_r, color = Fermenter, group = interaction(well,Fermenter)))+
    geom_point()+
    facet_wrap(carbon_source~temperature)+
    geom_line(aes(tmp, .fitted),preds)+
    theme_bw()


# Some curves with predicted values significantly outside of observed range. 
# Remove based on estimated rmax and topt


out_growth <- preds$ID[preds$.fitted > .65] %>% unique
# 2 curves removed, split across temps and carbon. One respirator, one fermenter 

preds_bound <- preds[!ID %in% out_growth,]
dt.parms.bound <- dt.parms[!ID %in% out_growth,]
coefs.bound <- coefs[!ID %in% out_growth,]

# Get mean fit for each functional group, temperature, and community of origin
mean_fit <- preds_bound[, .(mean_fit = mean(.fitted)), by = c("carbon_source", "temperature", "Fermenter", "tmp")]
mean_fit$well <- "A0"

pdf("Plots/Thermal_gaussian.pdf", 8, 6)
ggplot(dt.parms.bound, aes(x = tmp, y = mean_r, color = Fermenter, group = interaction(well,Fermenter)))+
    geom_point(size=0.2)+
    facet_wrap(carbon_source~temperature)+
    geom_line(aes(tmp, .fitted),preds_bound, alpha=0.8, linewidth= 0.1)+
    geom_line(aes(tmp, mean_fit), mean_fit, linewidth= 0.4)+
    scale_color_manual(values = c("black","gray65"))+
    theme_bw()
dev.off()

pdf("Plots/Topt_by_Tisolation.pdf", 6, 4)
ggplot(coefs.bound, aes(x = as.numeric(temp_origin), y = topt, color = fermenter))+
    geom_point(shape=21)+
    scale_color_manual(values = c("black","gray"))+
    xlab("Temperature of isolation [°C]")+
    ylab(expression(T[opt]*"[°C]"))+
    stat_summary(fun.y= mean, shape =19, size = 0.5)+
    stat_smooth(method="lm", se=F)+
    #stat_summary(fun.y= median, geom="line")+
    #stat_summary(fun.y= median, shape =22, size = 0.5)+
    geom_abline(slope = 1, intercept = 0, linetype=2)+
    theme_classic()
dev.off()

# Fit linear model 
coefs$temp_origin <- as.numeric(coefs$temp_origin)
mod_topt <- lm(temp_origin ~ fermenter*topt , data = coefs.bound)
summary(mod_topt)


par(mfrow=c(2,2))
plot(mod_topt)

confint(mod_topt)

# Summarize by group
coefs.bound[,.(mean(topt), sd(topt)),by="fermenter"]

a1 <- lm(rmax ~ temp_origin * carbon * fermenter, data = coefs.bound)
anova(a1)
coefs.bound[,.N, by=c("fermenter", "carbon", "temp_origin")]

# To evaluate effect of sample size inbalances remove 
sub <- coefs.bound[!(temp_origin=="22" & fermenter=="Respirator")&
    !(temp_origin=="30" & fermenter=="Respirator")]

sub[,.N, by=c("fermenter", "carbon", "temp_origin")]
a2 <- lm(rmax ~ temp_origin * carbon * fermenter, data = sub)
anova(a2)