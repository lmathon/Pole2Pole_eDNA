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
library(performance)
library(relaimpo)


load("Rdata/richness_station.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(rich_station) <- rich_station$station

data <- left_join(exp_var, rich_station[,c("largefish_MOTUs", "station")], by="station")
data <- data %>%
  dplyr::select(-c(station))

hist(data$largefish_MOTUs, main = "MOTUs_largefish", xlab ="MOTUs_largefish")

data$largefish_MOTUs <- log10(data$largefish_MOTUs+1)
hist(data$largefish_MOTUs, main = "log10(MOTUs_largefish+1)", xlab ="log10(MOTUs_largefish+1)")

# join longitude & latitude
meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(rich_station),]
identical(as.character(rownames(meta)), rownames(rich_station))
coor <- meta[, c("longitude_start", "latitude_start")]
data <- cbind(data, coor)
data$sample_method2 <- as.factor(data$sample_method2)

#### GLS to account for spatial autocorrelation ####

mexp <- gls(largefish_MOTUs ~ . -latitude_start - longitude_start, correlation = corExp(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mgau <- gls(largefish_MOTUs ~ . -latitude_start - longitude_start, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

msph <- gls(largefish_MOTUs ~ . -latitude_start - longitude_start, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mlin <- gls(largefish_MOTUs ~ . -latitude_start - longitude_start, correlation = corLin(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mrat <- gls(largefish_MOTUs ~ . -latitude_start - longitude_start, correlation = corRatio(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")


# Info full model
AIC(mexp, mgau, mlin, msph, mrat)

gls.full <- mgau

# remove colinear variables from VIF
gls.largefish <- gls(largefish_MOTUs ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+Gravity+MarineEcosystemDependency+dist_to_CT+bathy+depth_sampling+distCoast+volume+sample_method2, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

save(gls.largefish, file="Rdata/gls_large.rdata")

summary(gls.largefish)
anova(gls.largefish, type = "marginal")
AIC(gls.largefish)


# R² for GLS
r2(gls.largefish)


shapiro.test(gls.largefish$residuals)
hist(gls.largefish$residuals)



fit.grav.large <- visreg(gls.largefish,"Gravity",scale="response")
save(fit.grav.large, file="Rdata/fit.grav.large.rdata")
fit.SST.large <- visreg(gls.largefish,"mean_SST_1year",scale="response")
save(fit.SST.large, file="Rdata/fit.SST.large.rdata")
fit.SSS.large <- visreg(gls.largefish,"mean_sss_1year",scale="response")
save(fit.SSS.large, file="Rdata/fit.SSS.large.rdata")
fit.MED.large <- visreg(gls.largefish,"MarineEcosystemDependency",scale="response")
save(fit.MED.large, file="Rdata/fit.MED.large.rdata")
fit.CT.large <- visreg(gls.largefish,"dist_to_CT",scale="response")
save(fit.CT.large, file="Rdata/fit.CT.large.rdata")
fit.coast.large <- visreg(gls.largefish,"distCoast",scale="response")
save(fit.coast.large, file="Rdata/fit.coast.large.rdata")
fit.vol.large <- visreg(gls.largefish,"volume",scale="response")
save(fit.vol.large, file="Rdata/fit.vol.large.rdata")
fit.samp.large <- visreg(gls.largefish,"depth_sampling",scale="response")
save(fit.samp.large, file="Rdata/fit.samp.large.rdata")
fit.DHW.large <- visreg(gls.largefish,"mean_DHW_1year",scale="response")
save(fit.DHW.large, file="Rdata/fit.DHW.large.rdata")
fit.method.large <- visreg(gls.largefish,"sample_method2",scale="response")
save(fit.method.large, file="Rdata/fit.method.large.rdata")


fit.grav_med.large <- visreg2d(gls.largefish, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", xlab="log10(Gravity +1)", zlab="Large fish\na-diversity", plot.type="gg", color=c("red", "white", "blue"))
save(fit.grav_med.large, file="Rdata/fit.grav_med.large.rdata")

#### part R² ####
relimpo <- calc.relimp(largefish_MOTUs ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+Gravity+MarineEcosystemDependency+dist_to_CT+bathy+depth_sampling+distCoast+volume+sample_method2,  
                       data, type = c("lmg", "last", "first"))

r2_large <- as.data.frame(relimpo$lmg)



# boxplot partition per variable type

partition <- data.frame(environment=sum(r2_large[1:4,]), 
                        geography=sum(r2_large[8:11,]), 
                        socioeconomy=sum(r2_large[5:7,]), 
                        sampling=r2_large[c(12,13),])

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("environment", "geography", "socioeconomy", "sampling"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("partial R²")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())

#### effect size ####

large_effectsize <- effectsize(gls.largefish)
large_effectsize <- large_effectsize[-1,]
large_effectsize$taxa <- "Richness - Large fish"
large_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","geography","geography","geography","geography","sampling","sampling")

save(large_effectsize, file = "Rdata/large_effectsize.rdata")

