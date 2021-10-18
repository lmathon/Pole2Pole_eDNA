library(vegan)
library(visreg)
library(lmeInfo)
library(tidyverse)
library(nlme)
library(rcompanion)
library(ggplot2)
library(effectsize)


load("Rdata/richness_station.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(rich_station) <- rich_station$station

data <- left_join(exp_var_num, rich_station[,c("largefish_MOTUs", "station")], by="station")
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
gls.largefish <- gls(largefish_MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Corruption_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

summary(gls.largefish)
anova(gls.largefish, type = "marginal")
AIC(gls.largefish)


# R² for GLS
MOTU_pred <- predict(gls.largefish)
fit <- lm(MOTU_pred ~ data$largefish_MOTUs)
RsquareAdj(fit)

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
fit.DHW.large <- visreg(gls.largefish,"mean_DHW_5year",scale="response")
save(fit.DHW.large, file="Rdata/fit.DHW.large.rdata")

fit.grav_med.large <- visreg2d(gls.largefish, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(Large fish richness +1)", xlab="lg10(Gravity +1)")
save(fit.grav_med.large, file="Rdata/fit.grav_med.large.rdata")


#### Variation partitioning ####
env_var <- data[,c("mean_DHW_1year", "mean_DHW_5year", "mean_sss_1year", "mean_SST_1year", "mean_npp_1year")]
geo_var <- data[, c("bathy", "dist_to_CT", "distCoast","depth_sampling")]
socio_var <- data[,c("HDI2019", "conflicts", "Corruption_mean", "Gravity", "MarineEcosystemDependency")]
samp_var <- data[, c("volume")]

varpart <- varpart(gls.largefish$fitted, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))


# boxplot partition per variable type

partition <- data.frame(environment=0.194+0.004+0.217+0.157+0.112+0.082+0.022+0.098, 
                        geography=0.076+0.004+0.004+0.217+0.082+0.022+0.004, 
                        socioeconomy=0.026+0.01+0.004+0.217+0.082+0.157+0.112, 
                        sampling=0.002+0.01+0.004+0.082+0.022+0.112+0.098)

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("environment", "socioeconomy", "geography", "sampling"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("cumulated variance explained")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())

#### effect size ####

large_effectsize <- effectsize(gls.largefish)
large_effectsize <- large_effectsize[-1,]
large_effectsize$taxa <- "Large fish"
large_effectsize$vargroup <- c("environment","environment","environment","environment","environment","socio","socio","socio","socio","socio","geography","geography","geography","geography","sampling")

save(large_effectsize, file = "Rdata/large_effectsize.rdata")

