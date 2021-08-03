library(agricolae)
library(ggpubr)


load("Rdata/all_explanatory_variables_categorical.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
load("Rdata/richness_station.rdata")

data <- left_join(exp_var_num, rich_station, by="station")

# Kruskal-Wallis with categorical variables
shapiro.test(rich_station$MOTUs) #  distribution pas normale
bartlett.test(rich_station$MOTUs, exp_var_cat$sample_method) # non homocedasticite des variances
kruskal.test(rich_station$MOTUs, exp_var_cat$sample_method) # differences
kruskal(rich_station$MOTUs, exp_var_cat$sample_method, group=TRUE, p.adj="bonferroni")$groups

bartlett.test(rich_station$MOTUs, exp_var_cat$sequencer) # non homocedasticite des variances
kruskal.test(rich_station$MOTUs, exp_var_cat$sequencer) #  differences
kruskal(rich_station$MOTUs, exp_var_cat$sequencer, group=TRUE, p.adj="bonferroni")$groups

bartlett.test(rich_station$MOTUs, exp_var_cat$province) # non homocedasticite des variances
kruskal.test(rich_station$MOTUs, exp_var_cat$province) #  differences
kruskal(rich_station$MOTUs, exp_var_cat$province, group=TRUE, p.adj="bonferroni")$groups


# linear models with numerical variables
lm_DHW1 <- lm(rich_station$MOTUs~exp_var_num$mean_DHW_1year)
summary(lm_DHW1)
plot(rich_station$MOTUs~exp_var_num$mean_DHW_1year)

lm_DHW5 <- lm(rich_station$MOTUs~exp_var_num$mean_DHW_5year)
summary(lm_DHW5)
plot(rich_station$MOTUs~exp_var_num$mean_DHW_5year)

lm_SSS <- lm(rich_station$MOTUs~exp_var_num$mean_sss_1year)
summary(lm_SSS)
sss <- ggplot(data, aes(mean_sss_1year, MOTUs))+
  geom_point()+
  geom_abline(slope = 13.494, intercept = -403.558, size=0.8)+
  xlim(0,40)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA))
  

lm_SST <- lm(rich_station$MOTUs~exp_var_num$mean_SST_1year)
summary(lm_SST)
sst <- ggplot(data, aes(mean_SST_1year, MOTUs))+
  geom_point()+
  geom_abline(slope = 2.6647, intercept = 16.1552, size=0.8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA))

lm_npp <- lm(rich_station$MOTUs~exp_var_num$mean_npp_1year)
summary(lm_npp)

lm_ph <- lm(rich_station$MOTUs~exp_var_num$pH_mean)
summary(lm_ph)

lm_hdi <- lm(rich_station$MOTUs~exp_var_num$HDI2019)
summary(lm_hdi)
plot(rich_station$MOTUs~exp_var_num$HDI2019)


lm_neartt <- lm(rich_station$MOTUs~exp_var_num$neartt)
summary(lm_neartt)

lm_gravity <- lm(rich_station$MOTUs~exp_var_num$Gravity)
summary(lm_gravity)
gravity <- ggplot(data, aes(Gravity, MOTUs))+
  geom_point()+
  geom_abline(slope = 0.08298, intercept = 59.53905, size=0.8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA))

lm_ngo <- lm(rich_station$MOTUs~exp_var_num$NGO)
summary(lm_ngo)
ngo <- ggplot(data, aes(NGO, MOTUs))+
  geom_point()+
  geom_abline(slope = 0.036098, intercept = 35.197857, size=0.8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA))

lm_MarEco <- lm(rich_station$MOTUs~exp_var_num$MarineEcosystemDependency)
summary(lm_MarEco)

lm_natur <- lm(rich_station$MOTUs~exp_var_num$Naturalresourcesrents)
summary(lm_natur)
natur <- ggplot(data, aes(Naturalresourcesrents, MOTUs))+
  geom_point()+
  geom_abline(slope = -9.261, intercept = 98.112, size=0.8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA))

lm_distCT <- lm(rich_station$MOTUs~exp_var_num$dist_to_CT)
summary(lm_distCT)
distCT <- ggplot(data, aes(dist_to_CT, MOTUs))+
  geom_point()+
  geom_abline(slope = -0.0024369, intercept = 88.9153977, size=0.8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA))

lm_bathy <- lm(rich_station$MOTUs~exp_var_num$bathy)
summary(lm_bathy)
bathy <- ggplot(data, aes(bathy, MOTUs))+
  geom_point()+
  geom_abline(slope = 0.013833, intercept = 67.219780, size=0.8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA))

lm_depth <- lm(rich_station$MOTUs~exp_var_num$depth_sampling)
summary(lm_depth)

lm_latitude <- lm(rich_station$MOTUs~exp_var_num$latitude)
summary(lm_latitude)
latitude <- ggplot(data, aes(latitude, MOTUs))+
  geom_point()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA))

lm_distCoast <- lm(rich_station$MOTUs~exp_var_num$distCoast)
summary(lm_distCoast)
distcoast <- ggplot(data, aes(distCoast, MOTUs))+
  geom_point()+
  geom_abline(slope = -0.0002168, intercept = 66.9570458, size=0.8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA))

lm_volume <- lm(rich_station$MOTUs~exp_var_num$volume)
summary(lm_volume)
volume <- ggplot(data, aes(volume, MOTUs))+
  geom_point()+
  geom_abline(slope = 0.8368, intercept = 29.7677, size=0.8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA))

ggarrange(sss, sst, gravity, ngo, natur, distCT, bathy, latitude, distcoast, volume, nrow=4, ncol=3)
