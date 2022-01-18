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
library(performance)
library(relaimpo)

load("Rdata/FD_Hill_alpha.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(FD_Hill) <- FD_Hill$station

data <- left_join(exp_var, FD_Hill[,c("FD_q0", "station")], by="station")
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
data$sample_method2 <- as.factor(data$sample_method2)

#### GLS to account for spatial autocorrelation ####

mexp <- gls(FD_q0 ~ . -latitude_start - longitude_start, correlation = corExp(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mgau <- gls(FD_q0 ~ . -latitude_start - longitude_start, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

msph <- gls(FD_q0 ~ . -latitude_start - longitude_start, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mlin <- gls(FD_q0 ~ . -latitude_start - longitude_start, correlation = corLin(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mrat <- gls(FD_q0 ~ . -latitude_start - longitude_start, correlation = corRatio(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")


# Info full model
AIC(mexp, mgau, msph, mlin, mrat)

gls.full <- mgau

# remove colinear variables from VIF
gls.FDq0 <- gls(FD_q0 ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+Gravity+MarineEcosystemDependency+dist_to_CT+bathy+depth_sampling+distCoast+volume+sample_method2, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

save(gls.FDq0, file="Rdata/gls_FDq0.rdata")

AIC(gls.FDq0)
summary(gls.FDq0)
anova(gls.FDq0, type = "marginal")

# R² for GLS
r2(gls.FDq0)

hist(gls.FDq0$residuals)


fit.grav.FDq0 <- visreg(gls.FDq0,"Gravity",scale="response")
save(fit.grav.FDq0, file="Rdata/fit.grav.FDq0.rdata")
fit.SST.FDq0 <- visreg(gls.FDq0,"mean_SST_1year",scale="response")
save(fit.SST.FDq0, file="Rdata/fit.SST.FDq0.rdata")
fit.SSS.FDq0 <- visreg(gls.FDq0,"mean_sss_1year",scale="response")
save(fit.SSS.FDq0, file="Rdata/fit.SSS.FDq0.rdata")
fit.MED.FDq0 <- visreg(gls.FDq0,"MarineEcosystemDependency",scale="response")
save(fit.MED.FDq0, file="Rdata/fit.MED.FDq0.rdata")
fit.CT.FDq0 <- visreg(gls.FDq0,"dist_to_CT",scale="response")
save(fit.CT.FDq0, file="Rdata/fit.CT.FDq0.rdata")
fit.coast.FDq0 <- visreg(gls.FDq0,"distCoast",scale="response")
save(fit.coast.FDq0, file="Rdata/fit.coast.FDq0.rdata")
fit.vol.FDq0 <- visreg(gls.FDq0,"volume",scale="response")
save(fit.vol.FDq0, file="Rdata/fit.vol.FDq0.rdata")
fit.samp.FDq0 <- visreg(gls.FDq0,"depth_sampling",scale="response")
save(fit.samp.FDq0, file="Rdata/fit.samp.FDq0.rdata")
fit.DHW.FDq0 <- visreg(gls.FDq0,"mean_DHW_1year",scale="response")
save(fit.DHW.FDq0, file="Rdata/fit.DHW.FDq0.rdata")
fit.method.FDq0 <- visreg(gls.FDq0,"sample_method2",scale="response")
save(fit.method.FDq0, file="Rdata/fit.method.FDq0.rdata")

fit.grav_med.FDq0 <- visreg2d(gls.FDq0, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(FDq0 richness +1)", xlab="log10(Gravity +1)")
save(fit.grav_med.FDq0, file="Rdata/fit.grav_med.FDq0.rdata")


#### part R² ####
relimpo <- calc.relimp(FD_q0 ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+Gravity+MarineEcosystemDependency+dist_to_CT+bathy+depth_sampling+distCoast+volume+sample_method2,  
                       data, type = c("lmg", "last", "first"))

r2_FDq0 <- as.data.frame(relimpo$lmg)



# boxplot partition per variable type

partition <- data.frame(environment=sum(r2_FDq0[1:4,]), 
                        geography=sum(r2_FDq0[8:11,]), 
                        socioeconomy=sum(r2_FDq0[5:7,]), 
                        sampling=r2_FDq0[c(12,13),])

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c( "environment", "geography", "socioeconomy",  "sampling"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("partial R²")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())


#### effect size ####

FDq0_effectsize <- effectsize(gls.FDq0)
FDq0_effectsize <- FDq0_effectsize[-1,]
FDq0_effectsize$taxa <- "Functional a-diversity"
FDq0_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","geography","geography","geography","geography","sampling","sampling")

save(FDq0_effectsize, file = "Rdata/FDq0_effectsize.rdata")
