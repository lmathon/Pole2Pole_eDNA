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

data <- left_join(exp_var, rich_station[,c("crypto_MOTUs", "station")], by="station")
data <- data %>%
  dplyr::select(-c(station))

hist(data$crypto_MOTUs, main = "MOTUs_crypto", xlab ="MOTUs_crypto")

data$crypto_MOTUs <- log10(data$crypto_MOTUs+1)
hist(data$crypto_MOTUs, main = "log10(MOTUs_crypto+1)", xlab ="log10(MOTUs_crypto+1)")

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

mexp <- gls(crypto_MOTUs ~ . -latitude_start - longitude_start, correlation = corExp(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mgau <- gls(crypto_MOTUs ~ . -latitude_start - longitude_start, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

msph <- gls(crypto_MOTUs ~ . -latitude_start - longitude_start, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mlin <- gls(crypto_MOTUs ~ . -latitude_start - longitude_start, correlation = corLin(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mrat <- gls(crypto_MOTUs ~ . -latitude_start - longitude_start, correlation = corRatio(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")


# Info full model
AIC(mexp, mgau, mlin, msph, mrat)

gls.full <- mgau

# remove colinear variables from VIF
gls.crypto <- gls(crypto_MOTUs ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+Gravity+MarineEcosystemDependency+dist_to_CT+bathy+depth_sampling+distCoast+volume+sample_method2, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

save(gls.crypto, file="Rdata/gls_crypto.rdata")

summary(gls.crypto)
anova(gls.crypto, type = "marginal")
AIC(gls.crypto)


# R² for GLS
r2(gls.crypto)

shapiro.test(gls.crypto$residuals)
hist(gls.crypto$residuals)



fit.grav.crypto <- visreg(gls.crypto,"Gravity",scale="response")
save(fit.grav.crypto, file="Rdata/fit.grav.crypto.rdata")
fit.SST.crypto <- visreg(gls.crypto,"mean_SST_1year",scale="response")
save(fit.SST.crypto, file="Rdata/fit.SST.crypto.rdata")
fit.SSS.crypto <- visreg(gls.crypto,"mean_sss_1year",scale="response")
save(fit.SSS.crypto, file="Rdata/fit.SSS.crypto.rdata")
fit.MED.crypto <- visreg(gls.crypto,"MarineEcosystemDependency",scale="response")
save(fit.MED.crypto, file="Rdata/fit.MED.crypto.rdata")
fit.CT.crypto <- visreg(gls.crypto,"dist_to_CT",scale="response")
save(fit.CT.crypto, file="Rdata/fit.CT.crypto.rdata")
fit.coast.crypto <- visreg(gls.crypto,"distCoast",scale="response")
save(fit.coast.crypto, file="Rdata/fit.coast.crypto.rdata")
fit.vol.crypto <- visreg(gls.crypto,"volume",scale="response")
save(fit.vol.crypto, file="Rdata/fit.vol.crypto.rdata")
fit.samp.crypto <- visreg(gls.crypto,"depth_sampling",scale="response")
save(fit.samp.crypto, file="Rdata/fit.samp.crypto.rdata")
fit.DHW.crypto <- visreg(gls.crypto,"mean_DHW_1year",scale="response")
save(fit.DHW.crypto, file="Rdata/fit.DHW.crypto.rdata")
fit.method.crypto <- visreg(gls.crypto,"sample_method2",scale="response", xlab="sample method", ylab="Cryptobenthic taxonomic diversity", gg=TRUE)
save(fit.method.crypto, file="Rdata/fit.method.crypto.rdata")

fit.grav_med.crypto <- visreg2d(gls.crypto, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", xlab="log10(Gravity +1)", zlab=expression(paste("Cryptobenthic ",alpha,"-diversity")), plot.type="gg", color=c("red", "white", "blue"))
save(fit.grav_med.crypto, file="Rdata/fit.grav_med.crypto.rdata")

#### part R² ####
relimpo <- calc.relimp(crypto_MOTUs ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+Gravity+MarineEcosystemDependency+dist_to_CT+bathy+depth_sampling+distCoast+volume+sample_method2, 
                       data, type = c("lmg", "last", "first"))

r2_crypto <- as.data.frame(relimpo$lmg)



# boxplot partition per variable type

partition <- data.frame(environment=sum(r2_crypto[1:4,]), 
                        geography=sum(r2_crypto[8:11,]), 
                        socioeconomy=sum(r2_crypto[5:7,]), 
                        sampling=r2_crypto[c(12,13),])

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("environment", "geography", "socioeconomy", "sampling"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("partial R²")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())

#### effect size ####

crypto_effectsize <- effectsize(gls.crypto)
crypto_effectsize <- crypto_effectsize[-1,]
crypto_effectsize$taxa <- "Richness - Cryptobenthics"
crypto_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","geography","geography","geography","geography","sampling","sampling")

save(crypto_effectsize, file = "Rdata/crypto_effectsize.rdata")

