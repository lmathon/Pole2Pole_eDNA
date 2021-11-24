library(ggplot2)
library(grid)
library(cowplot)
library(tidyverse)

#### Gravity ####

load("Rdata/fit.grav.motus.rdata")
load("Rdata/fit.grav.large.rdata")
load("Rdata/fit.grav.crypto.rdata")
load("Rdata/fit.grav.MNTD.rdata")


ylim.prim <- c(0, 2)   
ylim.sec <- c(-1, 0.5) 

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

gravity <- ggplot() + 
  geom_line(data=fit.grav.motus$fit, aes(Gravity, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.motus$fit, aes(x=Gravity, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.grav.crypto$fit, aes(Gravity, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.crypto$fit, aes(x=Gravity, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.grav.MNTD$fit, aes(Gravity, (visregFit*b)+a), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.MNTD$fit, aes(x=Gravity, ymin=(visregLwr*b)+a, ymax=(visregUpr*b)+a), fill='orange', alpha=.2)+
  geom_line(data=fit.grav.large$fit, aes(Gravity, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.large$fit, aes(x=Gravity, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  labs(x="log10(Gravity+1)")+
  scale_y_continuous("log10(MOTUs richness+1)", sec.axis = sec_axis(~ (. - a)/b, name = "ses.MNTD")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

ggsave("outputs/GLS/visreg_gravity.png", height = 5, width = 5)  

#### SST ####

load("Rdata/fit.SST.motus.rdata")
load("Rdata/fit.SST.large.rdata")
load("Rdata/fit.SST.crypto.rdata")
load("Rdata/fit.SST.MNTD.rdata")

sst <- ggplot() + 
  geom_line(data=fit.SST.motus$fit, aes(mean_SST_1year, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.SST.motus$fit, aes(x=mean_SST_1year, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.SST.crypto$fit, aes(mean_SST_1year, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.SST.crypto$fit, aes(x=mean_SST_1year, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.SST.MNTD$fit, aes(mean_SST_1year, (visregFit*b)+a), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.SST.MNTD$fit, aes(x=mean_SST_1year, ymin=(visregLwr*b)+a, ymax=(visregUpr*b)+a), fill='orange', alpha=.2)+
  geom_line(data=fit.SST.large$fit, aes(mean_SST_1year, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.SST.large$fit, aes(x=mean_SST_1year, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  annotate(geom="text", x=0, y=2, label="all MOTUs", hjust=0, size=3, color="darkred", fontface = "bold")+
  annotate(geom="text", x=0, y=2.1, label="ses.MNTD", hjust=0, size=3, color="orange", fontface = "bold")+
  annotate(geom="text", x=0, y=1.9, label="Crypto", hjust=0, size=3, color="navy", fontface = "bold")+
  annotate(geom="text", x=0, y=1.8, label="Large fish", hjust=0, size=3, color="darkgreen", fontface = "bold")+
  labs(x="mean SST")+
  scale_y_continuous("log10(MOTUs richness+1)", sec.axis = sec_axis(~ (. - a)/b, name = "ses.MNTD")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

ggsave("outputs/GLS/visreg_SST.png", height = 5, width = 5)  

#### SSS ####

load("Rdata/fit.SSS.motus.rdata")
load("Rdata/fit.SSS.large.rdata")
load("Rdata/fit.SSS.crypto.rdata")
load("Rdata/fit.SSS.MNTD.rdata")

sss <- ggplot() + 
  geom_line(data=fit.SSS.motus$fit, aes(mean_sss_1year, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.SSS.motus$fit, aes(x=mean_sss_1year, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.SSS.crypto$fit, aes(mean_sss_1year, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.SSS.crypto$fit, aes(x=mean_sss_1year, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.SSS.MNTD$fit, aes(mean_sss_1year, (visregFit*b)+a), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.SSS.MNTD$fit, aes(x=mean_sss_1year, ymin=(visregLwr*b)+a, ymax=(visregUpr*b)+a), fill='orange', alpha=.2)+
  geom_line(data=fit.SSS.large$fit, aes(mean_sss_1year, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.SSS.large$fit, aes(x=mean_sss_1year, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  labs(x="mean SSS")+
  scale_y_continuous("log10(MOTUs richness+1)", sec.axis = sec_axis(~ (. - a)/b, name = "ses.MNTD")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

ggsave("outputs/GLS/visreg_SSS.png", height = 5, width = 5)  

#### MED ####

load("Rdata/fit.MED.motus.rdata")
load("Rdata/fit.MED.large.rdata")
load("Rdata/fit.MED.crypto.rdata")
load("Rdata/fit.MED.MNTD.rdata")

med <- ggplot() + 
  geom_line(data=fit.MED.motus$fit, aes(MarineEcosystemDependency, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.MED.motus$fit, aes(x=MarineEcosystemDependency, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.MED.crypto$fit, aes(MarineEcosystemDependency, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.MED.crypto$fit, aes(x=MarineEcosystemDependency, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.MED.MNTD$fit, aes(MarineEcosystemDependency, (visregFit*b)+a), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.MED.MNTD$fit, aes(x=MarineEcosystemDependency, ymin=(visregLwr*b)+a, ymax=(visregUpr*b)+a), fill='orange', alpha=.2)+
  geom_line(data=fit.MED.large$fit, aes(MarineEcosystemDependency, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.MED.large$fit, aes(x=MarineEcosystemDependency, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  labs(x="Marine Ecosystem Dependency")+
  scale_y_continuous("log10(MOTUs richness+1)", sec.axis = sec_axis(~ (. - a)/b, name = "ses.MNTD")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

ggsave("outputs/GLS/visreg_MED.png", height = 5, width = 5)  

#### Dist CT ####

load("Rdata/fit.CT.motus.rdata")
load("Rdata/fit.CT.large.rdata")
load("Rdata/fit.CT.crypto.rdata")
load("Rdata/fit.CT.MNTD.rdata")

distCT <- ggplot() + 
  geom_line(data=fit.CT.motus$fit, aes(dist_to_CT, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.CT.motus$fit, aes(x=dist_to_CT, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.CT.crypto$fit, aes(dist_to_CT, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.CT.crypto$fit, aes(x=dist_to_CT, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.CT.MNTD$fit, aes(dist_to_CT, (visregFit*b)+a), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.CT.MNTD$fit, aes(x=dist_to_CT, ymin=(visregLwr*b)+a, ymax=(visregUpr*b)+a), fill='orange', alpha=.2)+
  geom_line(data=fit.CT.large$fit, aes(dist_to_CT, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.CT.large$fit, aes(x=dist_to_CT, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  labs(x="log10(Distance to Coral Triangle (km) +1)")+
  scale_y_continuous("log10(MOTUs richness+1)", sec.axis = sec_axis(~ (. - a)/b, name = "ses.MNTD")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

ggsave("outputs/GLS/visreg_distCT.png", height = 5, width = 5)  

#### Dist Coast ####

load("Rdata/fit.coast.motus.rdata")
load("Rdata/fit.coast.large.rdata")
load("Rdata/fit.coast.crypto.rdata")
load("Rdata/fit.coast.MNTD.rdata")

distcoast <- ggplot() + 
  geom_line(data=fit.coast.motus$fit, aes(distCoast, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.coast.motus$fit, aes(x=distCoast, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.coast.crypto$fit, aes(distCoast, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.coast.crypto$fit, aes(x=distCoast, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.coast.MNTD$fit, aes(distCoast, (visregFit*b)+a), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.coast.MNTD$fit, aes(x=distCoast, ymin=(visregLwr*b)+a, ymax=(visregUpr*b)+a), fill='orange', alpha=.2)+
  geom_line(data=fit.coast.large$fit, aes(distCoast, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.coast.large$fit, aes(x=distCoast, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  labs(x="log10(Distance to coast (m) +1)")+
  scale_y_continuous("log10(MOTUs richness+1)", sec.axis = sec_axis(~ (. - a)/b, name = "ses.MNTD")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

ggsave("outputs/GLS/visreg_distCoast.png", height = 5, width = 5)  

#### volume ####

load("Rdata/fit.vol.motus.rdata")
load("Rdata/fit.vol.large.rdata")
load("Rdata/fit.vol.crypto.rdata")
load("Rdata/fit.vol.MNTD.rdata")

volume <- ggplot() + 
  geom_line(data=fit.vol.motus$fit, aes(volume, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.vol.motus$fit, aes(x=volume, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.vol.crypto$fit, aes(volume, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.vol.crypto$fit, aes(x=volume, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.vol.MNTD$fit, aes(volume, (visregFit*b)+a), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.vol.MNTD$fit, aes(x=volume, ymin=(visregLwr*b)+a, ymax=(visregUpr*b)+a), fill='orange', alpha=.2)+
  geom_line(data=fit.vol.large$fit, aes(volume, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.vol.large$fit, aes(x=volume, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  labs(x="log10(volume (L) +1)")+
  scale_y_continuous("log10(MOTUs richness+1)", sec.axis = sec_axis(~ (. - a)/b, name = "ses.MNTD")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

ggsave("outputs/GLS/visreg_volume.png", height = 5, width = 5)  

#### depth sampling ####

load("Rdata/fit.samp.motus.rdata")
load("Rdata/fit.samp.large.rdata")
load("Rdata/fit.samp.crypto.rdata")
load("Rdata/fit.samp.MNTD.rdata")

depthsamp <- ggplot() + 
  geom_line(data=fit.samp.motus$fit, aes(depth_sampling, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.samp.motus$fit, aes(x=depth_sampling, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.samp.crypto$fit, aes(depth_sampling, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.samp.crypto$fit, aes(x=depth_sampling, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.samp.MNTD$fit, aes(depth_sampling, (visregFit*b)+a), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.samp.MNTD$fit, aes(x=depth_sampling, ymin=(visregLwr*b)+a, ymax=(visregUpr*b)+a), fill='orange', alpha=.2)+
  geom_line(data=fit.samp.large$fit, aes(depth_sampling, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.samp.large$fit, aes(x=depth_sampling, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  labs(x="depth of sampling (m)")+
  scale_y_continuous("log10(MOTUs richness+1)", sec.axis = sec_axis(~ (. - a)/b, name = "ses.MNTD")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

ggsave("outputs/GLS/visreg_depthsampling.png", height = 5, width = 5)  


#### DHW1 ####

load("Rdata/fit.DHW.motus.rdata")
load("Rdata/fit.DHW.large.rdata")
load("Rdata/fit.DHW.crypto.rdata")
load("Rdata/fit.DHW.MNTD.rdata")

dhw <- ggplot() + 
  geom_line(data=fit.DHW.motus$fit, aes(mean_DHW_1year, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.DHW.motus$fit, aes(x=mean_DHW_1year, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.DHW.crypto$fit, aes(mean_DHW_1year, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.DHW.crypto$fit, aes(x=mean_DHW_1year, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.DHW.MNTD$fit, aes(mean_DHW_1year, (visregFit*b)+a), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.DHW.MNTD$fit, aes(x=mean_DHW_1year, ymin=(visregLwr*b)+a, ymax=(visregUpr*b)+a), fill='orange', alpha=.2)+
  geom_line(data=fit.DHW.large$fit, aes(mean_DHW_1year, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.DHW.large$fit, aes(x=mean_DHW_1year, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  labs(x="log10(mean DHW 1year +1)")+
  scale_y_continuous("log10(MOTUs richness+1)", sec.axis = sec_axis(~ (. - a)/b, name = "ses.MNTD")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

ggsave("outputs/GLS/visreg_DHW1.png", height = 5, width = 5)  


#### plot all together ####

ggarrange(sst, distCT, gravity, med, nrow=2, ncol=2, labels = c("A", "B", "C", "D"))
ggsave("outputs/GLS/visreg_multiple.png", width = 8, height = 7)


#### Plot visreg gravity + MED ####

load("Rdata/gls_motus.rdata")
load("Rdata/gls_mntd.rdata")
load("Rdata/gls_crypto.rdata")
load("Rdata/gls_large.rdata")

load("Rdata/richness_station.rdata")
load("Rdata/MNTD_station.rdata")

load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(rich_station) <- rich_station$station

data <- left_join(exp_var_num, rich_station, by="station")
data <- left_join(data, mntd_stations)
data <- data %>%
  dplyr::select(-c(station))


visreg2d(gls.motus, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(MOTUs richness +1)", xlab="log10(Gravity +1)", plot.type="gg")
visreg2d(gls.MNTD, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="ses.MNTD", xlab="log10(Gravity +1)", plot.type="gg")
visreg2d(gls.crypto, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(Crypto richness +1)", xlab="log10(Gravity +1)", plot.type="gg")
visreg2d(gls.largefish, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(Large fish richness +1)", xlab="log10(Gravity +1)", plot.type="gg")

