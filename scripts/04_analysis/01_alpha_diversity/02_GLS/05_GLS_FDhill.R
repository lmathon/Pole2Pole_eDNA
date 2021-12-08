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
library(greekLetters)

load("Rdata/FD_Hill_alpha.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(FD_Hill) <- FD_Hill$station

data <- left_join(exp_var_num, FD_Hill[,c("FD_q2", "station")], by="station")
data <- data %>%
  dplyr::select(-c(station))


# join longitude & latitude
meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(FD_Hill),]
identical(as.character(rownames(meta)), rownames(FD_Hill))
coor <- meta[, c("longitude_start", "latitude_start", "province")]
data <- cbind(data, coor)


#### GLS to account for spatial autocorrelation ####

mexp <- gls(FD_q2 ~ . -latitude_start - longitude_start, correlation = corExp(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mgau <- gls(FD_q2 ~ . -latitude_start - longitude_start, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

msph <- gls(FD_q2 ~ . -latitude_start - longitude_start, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mlin <- gls(FD_q2 ~ . -latitude_start - longitude_start, correlation = corLin(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mrat <- gls(FD_q2 ~ . -latitude_start - longitude_start, correlation = corRatio(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")


# Info full model
AIC(mexp, mgau, msph, mlin, mrat)

gls.full <- mgau

# remove colinear variables from VIF
gls.FDq2 <- gls(FD_q2 ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+Gravity+MarineEcosystemDependency+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

save(gls.FDq2, file="Rdata/gls_FDq2.rdata")

AIC(gls.FDq2)
summary(gls.FDq2)
anova(gls.FDq2, type = "marginal")

# R² for GLS
FDq2_pred <- predict(gls.FDq2)
fit <- lm(FDq2_pred ~ data$FD_q2)
RsquareAdj(fit)

hist(gls.FDq2$residuals)


fit.grav.FDq2 <- visreg(gls.FDq2,"Gravity",scale="response")
save(fit.grav.FDq2, file="Rdata/fit.grav.FDq2.rdata")
fit.SST.FDq2 <- visreg(gls.FDq2,"mean_SST_1year",scale="response")
save(fit.SST.FDq2, file="Rdata/fit.SST.FDq2.rdata")
fit.SSS.FDq2 <- visreg(gls.FDq2,"mean_sss_1year",scale="response")
save(fit.SSS.FDq2, file="Rdata/fit.SSS.FDq2.rdata")
fit.MED.FDq2 <- visreg(gls.FDq2,"MarineEcosystemDependency",scale="response")
save(fit.MED.FDq2, file="Rdata/fit.MED.FDq2.rdata")
fit.CT.FDq2 <- visreg(gls.FDq2,"dist_to_CT",scale="response")
save(fit.CT.FDq2, file="Rdata/fit.CT.FDq2.rdata")
fit.coast.FDq2 <- visreg(gls.FDq2,"distCoast",scale="response")
save(fit.coast.FDq2, file="Rdata/fit.coast.FDq2.rdata")
fit.vol.FDq2 <- visreg(gls.FDq2,"volume",scale="response")
save(fit.vol.FDq2, file="Rdata/fit.vol.FDq2.rdata")
fit.samp.FDq2 <- visreg(gls.FDq2,"depth_sampling",scale="response")
save(fit.samp.FDq2, file="Rdata/fit.samp.FDq2.rdata")
fit.DHW.FDq2 <- visreg(gls.FDq2,"mean_DHW_1year",scale="response")
save(fit.DHW.FDq2, file="Rdata/fit.DHW.FDq2.rdata")


fit.grav_med.FDq2 <- visreg2d(gls.FDq2, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(FDq2 richness +1)", xlab="log10(Gravity +1)")
save(fit.grav_med.FDq2, file="Rdata/fit.grav_med.FDq2.rdata")

#### Variation partitioning ####
env_var <- data[,c("mean_DHW_1year", "mean_sss_1year", "mean_SST_1year", "mean_npp_1year")]
geo_var <- data[, c("bathy", "dist_to_CT", "distCoast","depth_sampling")]
socio_var <- data[,c("HDI2019", "Gravity", "MarineEcosystemDependency")]
samp_var <- data[, c("volume")]

varpart <- varpart(gls.FDq2$fitted, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))


# boxplot partition per variable type

partition <- data.frame(environment=0.297+0.031+0.208+0.123+0.008+0.081+0.226, 
                        geography=0.095+0.031+0.208+0.008+0.081+0.007, 
                        socioeconomy=0.021+0.02+0.208+0.008+0.123, 
                        sampling=0.008+0.02+0.007+0.008+0.081+0.226)

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c( "environment", "geography", "socioeconomy",  "sampling"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("cumulated variance explained")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())


#### effect size ####

FDq2_effectsize <- effectsize(gls.FDq2)
FDq2_effectsize <- FDq2_effectsize[-1,]
FDq2_effectsize$taxa <- "Functional a-diversity"
FDq2_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","geography","geography","geography","geography","sampling")

save(FDq2_effectsize, file = "Rdata/FDq2_effectsize.rdata")
