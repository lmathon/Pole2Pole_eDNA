library(funk)
library(vegan)
library(ade4)
library(lme4)
library(DHARMa)
library(car)
library(visreg)
library(ecospat)
library(modEvA)
library(psych)
library(MASS)
library(lmeInfo)
library(pscl)
library(spdep)
library(tidyverse)
library(nlme)
library(MuMIn)
library(rcompanion)
library(ggpubr)
library(ggplot2)
library(effectsize)

load("Rdata/richness_station.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(rich_station) <- rich_station$station

data <- left_join(exp_var_num, rich_station[,c("MOTUs", "station")], by="station")
data <- data %>%
  dplyr::select(-c(station))

hist(data$MOTUs, main = "MOTUs", xlab ="MOTUs")
data$MOTUs <- log10(data$MOTUs +1)

hist(data$MOTUs, main = "log10(MOTUs+1)", xlab ="log10(MOTUs+1)")

# join longitude & latitude
meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(rich_station),]
identical(as.character(rownames(meta)), rownames(rich_station))
coor <- meta[, c("longitude_start", "latitude_start", "province")]
data <- cbind(data, coor)


#### GLS to account for spatial autocorrelation ####

mexp <- gls(MOTUs ~ . -latitude_start - longitude_start, correlation = corExp(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mgau <- gls(MOTUs ~ . -latitude_start - longitude_start, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

msph <- gls(MOTUs ~ . -latitude_start - longitude_start, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mlin <- gls(MOTUs ~ . -latitude_start - longitude_start, correlation = corLin(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mrat <- gls(MOTUs ~ . -latitude_start - longitude_start, correlation = corRatio(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
            

# Info full model
AIC(mexp, mgau, msph, mlin, mrat)

gls.full <- mgau

# remove colinear variables from VIF
gls.motus <- gls(MOTUs ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Voice_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

save(gls.motus, file="Rdata/gls_motus.rdata")


AIC(gls.motus)
summary(gls.motus)
anova(gls.motus, type = "marginal")

# R² for GLS
MOTU_pred <- predict(gls.motus)
fit <- lm(MOTU_pred ~ data$MOTUs)
RsquareAdj(fit)

hist(gls.motus$residuals)


fit.grav.motus <- visreg(gls.motus,"Gravity",scale="response")
save(fit.grav.motus, file="Rdata/fit.grav.motus.rdata")
fit.SST.motus <- visreg(gls.motus,"mean_SST_1year",scale="response")
save(fit.SST.motus, file="Rdata/fit.SST.motus.rdata")
fit.SSS.motus <- visreg(gls.motus,"mean_sss_1year",scale="response")
save(fit.SSS.motus, file="Rdata/fit.SSS.motus.rdata")
fit.MED.motus <- visreg(gls.motus,"MarineEcosystemDependency",scale="response")
save(fit.MED.motus, file="Rdata/fit.MED.motus.rdata")
fit.CT.motus <- visreg(gls.motus,"dist_to_CT",scale="response")
save(fit.CT.motus, file="Rdata/fit.CT.motus.rdata")
fit.coast.motus <- visreg(gls.motus,"distCoast",scale="response")
save(fit.coast.motus, file="Rdata/fit.coast.motus.rdata")
fit.vol.motus <- visreg(gls.motus,"volume",scale="response")
save(fit.vol.motus, file="Rdata/fit.vol.motus.rdata")
fit.samp.motus <- visreg(gls.motus,"depth_sampling",scale="response")
save(fit.samp.motus, file="Rdata/fit.samp.motus.rdata")
fit.DHW.motus <- visreg(gls.motus,"mean_DHW_1year",scale="response")
save(fit.DHW.motus, file="Rdata/fit.DHW.motus.rdata")


fit.grav_med.motus <- visreg2d(gls.motus, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(MOTUs richness +1)", xlab="log10(Gravity +1)", plot.type="gg")
save(fit.grav_med.motus, file="Rdata/fit.grav_med.motus.rdata")

#### Variation partitioning ####
env_var <- data[,c("mean_DHW_1year", "mean_sss_1year", "mean_SST_1year", "mean_npp_1year")]
geo_var <- data[, c("bathy", "dist_to_CT", "distCoast","depth_sampling")]
socio_var <- data[,c("HDI2019", "conflicts", "Voice_mean", "Gravity", "MarineEcosystemDependency")]
samp_var <- data[, c("volume")]

varpart <- varpart(gls.motus$fitted, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))


# boxplot partition per variable type

partition <- data.frame(environment=0.142+0.247+0.052+0.152+0.157+0.041+0.016, 
                        geography=0.087+0.247+0.013+0.152+0.019+0.064, 
                        socioeconomy=0.035+0.013+0.042+0.247+0.052+0.157+0.152, 
                        sampling=0.017+0.042+0.064+0.157+0.152+0.016+0.041)

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("environment", "socioeconomy", "geography", "sampling"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("cumulated variance explained")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())


#### effect size ####

motus_effectsize <- effectsize(gls.motus)
motus_effectsize <- motus_effectsize[-1,]
motus_effectsize$taxa <- "Richness - all MOTUs"
motus_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","socio","socio","geography","geography","geography","geography","sampling")

save(motus_effectsize, file = "Rdata/motu_effectsize.rdata")
