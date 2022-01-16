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

data <- left_join(exp_var, FD_Hill[,c("FD_q1", "station")], by="station")
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

mexp <- gls(FD_q1 ~ . -latitude_start - longitude_start, correlation = corExp(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mgau <- gls(FD_q1 ~ . -latitude_start - longitude_start, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

msph <- gls(FD_q1 ~ . -latitude_start - longitude_start, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mlin <- gls(FD_q1 ~ . -latitude_start - longitude_start, correlation = corLin(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mrat <- gls(FD_q1 ~ . -latitude_start - longitude_start, correlation = corRatio(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")


# Info full model
AIC(mexp, mgau, msph, mlin, mrat)

gls.full <- mgau

# remove colinear variables from VIF
gls.FDq1 <- gls(FD_q1 ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+Gravity+MarineEcosystemDependency+dist_to_CT+bathy+depth_sampling+distCoast+volume+sample_method2, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

save(gls.FDq1, file="Rdata/gls_FDq1.rdata")

AIC(gls.FDq1)
summary(gls.FDq1)
anova(gls.FDq1, type = "marginal")

# R² for GLS
r2(gls.FDq1)

hist(gls.FDq1$residuals)


fit.grav.FDq1 <- visreg(gls.FDq1,"Gravity",scale="response")
save(fit.grav.FDq1, file="Rdata/fit.grav.FDq1.rdata")
fit.SST.FDq1 <- visreg(gls.FDq1,"mean_SST_1year",scale="response")
save(fit.SST.FDq1, file="Rdata/fit.SST.FDq1.rdata")
fit.SSS.FDq1 <- visreg(gls.FDq1,"mean_sss_1year",scale="response")
save(fit.SSS.FDq1, file="Rdata/fit.SSS.FDq1.rdata")
fit.MED.FDq1 <- visreg(gls.FDq1,"MarineEcosystemDependency",scale="response")
save(fit.MED.FDq1, file="Rdata/fit.MED.FDq1.rdata")
fit.CT.FDq1 <- visreg(gls.FDq1,"dist_to_CT",scale="response")
save(fit.CT.FDq1, file="Rdata/fit.CT.FDq1.rdata")
fit.coast.FDq1 <- visreg(gls.FDq1,"distCoast",scale="response")
save(fit.coast.FDq1, file="Rdata/fit.coast.FDq1.rdata")
fit.vol.FDq1 <- visreg(gls.FDq1,"volume",scale="response")
save(fit.vol.FDq1, file="Rdata/fit.vol.FDq1.rdata")
fit.samp.FDq1 <- visreg(gls.FDq1,"depth_sampling",scale="response")
save(fit.samp.FDq1, file="Rdata/fit.samp.FDq1.rdata")
fit.DHW.FDq1 <- visreg(gls.FDq1,"mean_DHW_1year",scale="response")
save(fit.DHW.FDq1, file="Rdata/fit.DHW.FDq1.rdata")
fit.method.FDq1 <- visreg(gls.FDq1,"sample_method2",scale="response")
save(fit.method.FDq1, file="Rdata/fit.method.FDq1.rdata")

fit.grav_med.FDq1 <- visreg2d(gls.FDq1, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(FDq1 richness +1)", xlab="log10(Gravity +1)")
save(fit.grav_med.FDq1, file="Rdata/fit.grav_med.FDq1.rdata")


#### part R² ####
relimpo <- calc.relimp(FD_q1 ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+Gravity+MarineEcosystemDependency+dist_to_CT+bathy+depth_sampling+distCoast+volume+sample_method2,  
                       data, type = c("lmg", "last", "first"))

r2_FDq1 <- as.data.frame(relimpo$lmg)



# boxplot partition per variable type

partition <- data.frame(environment=sum(r2_FDq1[1:4,]), 
                        geography=sum(r2_FDq1[8:11,]), 
                        socioeconomy=sum(r2_FDq1[5:7,]), 
                        sampling=r2_FDq1[c(12,13),])

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c( "environment", "geography", "socioeconomy",  "sampling"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("partial R²")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())


#### effect size ####

FDq1_effectsize <- effectsize(gls.FDq1)
FDq1_effectsize <- FDq1_effectsize[-1,]
FDq1_effectsize$taxa <- "Functional a-diversity"
FDq1_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","geography","geography","geography","geography","sampling","sampling")

save(FDq1_effectsize, file = "Rdata/FDq1_effectsize.rdata")
