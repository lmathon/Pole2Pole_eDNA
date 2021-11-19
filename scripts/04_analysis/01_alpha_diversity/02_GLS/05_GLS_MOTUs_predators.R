library(funk)
library(vegan)
library(ade4)
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

data <- left_join(exp_var_num, rich_station[,c("predator_MOTUs", "station")], by="station")
data <- data %>%
  dplyr::select(-c(station))

hist(data$predator_MOTUs, main = "MOTUs_predators", xlab ="MOTUs_predators")

data$predator_MOTUs <- log10(data$predator_MOTUs+1)
hist(data$predator_MOTUs, main = "log10(MOTUs_predators+1)", xlab ="log10(MOTUs_predators+1)")

# join longitude & latitude
meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(rich_station),]
identical(as.character(rownames(meta)), rownames(rich_station))
coor <- meta[, c("longitude_start", "latitude_start")]
data <- cbind(data, coor)


#### GLS to account for spatial autocorrelation ####

mexp <- gls(predator_MOTUs ~ . -latitude_start - longitude_start, correlation = corExp(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mgau <- gls(predator_MOTUs ~ . -latitude_start - longitude_start, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

msph <- gls(predator_MOTUs ~ . -latitude_start - longitude_start, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mlin <- gls(predator_MOTUs ~ . -latitude_start - longitude_start, correlation = corLin(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mrat <- gls(predator_MOTUs ~ . -latitude_start - longitude_start, correlation = corRatio(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")


# Info full model
AIC(mexp, mgau, mlin, msph, mrat)

gls.full <- mgau

# remove colinear variables from VIF
gls.predator <- gls(predator_MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Corruption_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

summary(gls.predator)
anova(gls.predator, type = "marginal")
AIC(gls.predator)


# R² for GLS
MOTU_pred <- predict(gls.predator)
fit <- lm(MOTU_pred ~ data$predator_MOTUs)
RsquareAdj(fit)

shapiro.test(gls.predator$residuals)
hist(gls.predator$residuals)



fit.grav.predator <- visreg(gls.predator,"Gravity",scale="response")
save(fit.grav.predator, file="Rdata/fit.grav.predator.rdata")
fit.SST.predator <- visreg(gls.predator,"mean_SST_1year",scale="response")
save(fit.SST.predator, file="Rdata/fit.SST.predator.rdata")
fit.SSS.predator <- visreg(gls.predator,"mean_sss_1year",scale="response")
save(fit.SSS.predator, file="Rdata/fit.SSS.predator.rdata")
fit.MED.predator <- visreg(gls.predator,"MarineEcosystemDependency",scale="response")
save(fit.MED.predator, file="Rdata/fit.MED.predator.rdata")
fit.CT.predator <- visreg(gls.predator,"dist_to_CT",scale="response")
save(fit.CT.predator, file="Rdata/fit.CT.predator.rdata")
fit.coast.predator <- visreg(gls.predator,"distCoast",scale="response")
save(fit.coast.predator, file="Rdata/fit.coast.predator.rdata")
fit.vol.predator <- visreg(gls.predator,"volume",scale="response")
save(fit.vol.predator, file="Rdata/fit.vol.predator.rdata")
fit.samp.predator <- visreg(gls.predator,"depth_sampling",scale="response")
save(fit.samp.predator, file="Rdata/fit.samp.predator.rdata")
fit.DHW.predator <- visreg(gls.predator,"mean_DHW_5year",scale="response")
save(fit.DHW.predator, file="Rdata/fit.DHW.predator.rdata")

fit.grav_med.predator <- visreg2d(gls.predator, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(predator fish richness +1)", xlab="lg10(Gravity +1)")
save(fit.grav_med.predator, file="Rdata/fit.grav_med.predator.rdata")


#### Variation partitioning ####
env_var <- data[,c("mean_DHW_1year", "mean_DHW_5year", "mean_sss_1year", "mean_SST_1year", "mean_npp_1year")]
geo_var <- data[, c("bathy", "dist_to_CT", "distCoast","depth_sampling")]
socio_var <- data[,c("HDI2019", "conflicts", "Corruption_mean", "Gravity", "MarineEcosystemDependency")]
samp_var <- data[, c("volume")]

varpart <- varpart(gls.predator$fitted, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))


# boxplot partition per variable type

partition <- data.frame(environment=0.207+0.042+0.215+0.244+0.074+0.054+0.049+0.031, 
                        geography=0.015+0.042+0.042+0.215+0.054+0.049+0.004, 
                        socioeconomy=0.022+0.042+0.215+0.244+0.074+0.054, 
                        sampling=0.004+0.049+0.054+0.074+0.031)

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("environment", "socioeconomy", "geography", "sampling"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("cumulated variance explained")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())

#### effect size ####

predator_effectsize <- effectsize(gls.predator)
predator_effectsize <- predator_effectsize[-1,]
predator_effectsize$taxa <- "Richness - Predators"
predator_effectsize$vargroup <- c("environment","environment","environment","environment","environment","socio","socio","socio","socio","socio","geography","geography","geography","geography","sampling")

save(predator_effectsize, file = "Rdata/predator_effectsize.rdata")

