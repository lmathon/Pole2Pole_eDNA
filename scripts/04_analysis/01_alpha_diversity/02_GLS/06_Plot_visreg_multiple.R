library(ggplot2)
library(pdp)
library(tidyverse)

#### Gravity ####

load("Rdata/fit.grav.motus.rdata")
load("Rdata/fit.grav.chondri.rdata")
load("Rdata/fit.grav.large.rdata")
load("Rdata/fit.grav.crypto.rdata")

ggplot() + 
  geom_line(data=fit.grav.motus$fit, aes(Gravity, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.motus$fit, aes(x=Gravity, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.grav.crypto$fit, aes(Gravity, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.crypto$fit, aes(x=Gravity, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.grav.chondri$fit, aes(Gravity, visregFit), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.chondri$fit, aes(x=Gravity, ymin=visregLwr, ymax=visregUpr), fill='orange', alpha=.2)+
  geom_line(data=fit.grav.large$fit, aes(Gravity, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.large$fit, aes(x=Gravity, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  annotate(geom="text", x=8, y=4, label="all MOTUs", hjust=1, size=4, color="darkred", fontface = "bold")+
  annotate(geom="text", x=8, y=3.85, label="Chondri", hjust=1, size=4, color="orange", fontface = "bold")+
  annotate(geom="text", x=8, y=3.7, label="Crypto", hjust=1, size=4, color="navy", fontface = "bold")+
  annotate(geom="text", x=8, y=3.55, label="Large fish", hjust=1, size=4, color="darkgreen", fontface = "bold")+
  labs(x="log(Gravity)", y="log(MOTUs richness)")+
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
load("Rdata/fit.SST.chondri.rdata")
load("Rdata/fit.SST.large.rdata")
load("Rdata/fit.SST.crypto.rdata")

ggplot() + 
  geom_line(data=fit.SST.motus$fit, aes(mean_SST_1year, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.SST.motus$fit, aes(x=mean_SST_1year, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.SST.crypto$fit, aes(mean_SST_1year, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.SST.crypto$fit, aes(x=mean_SST_1year, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.SST.chondri$fit, aes(mean_SST_1year, visregFit), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.SST.chondri$fit, aes(x=mean_SST_1year, ymin=visregLwr, ymax=visregUpr), fill='orange', alpha=.2)+
  geom_line(data=fit.SST.large$fit, aes(mean_SST_1year, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.SST.large$fit, aes(x=mean_SST_1year, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  annotate(geom="text", x=0, y=4, label="all MOTUs", hjust=0, size=4, color="darkred", fontface = "bold")+
  annotate(geom="text", x=0, y=3.8, label="Chondri", hjust=0, size=4, color="orange", fontface = "bold")+
  annotate(geom="text", x=0, y=3.6, label="Crypto", hjust=0, size=4, color="navy", fontface = "bold")+
  annotate(geom="text", x=0, y=3.4, label="Large fish", hjust=0, size=4, color="darkgreen", fontface = "bold")+
  labs(x="mean SST", y="log(MOTUs richness)")+
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
load("Rdata/fit.SSS.chondri.rdata")
load("Rdata/fit.SSS.large.rdata")
load("Rdata/fit.SSS.crypto.rdata")

ggplot() + 
  geom_line(data=fit.SSS.motus$fit, aes(mean_sss_1year, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.SSS.motus$fit, aes(x=mean_sss_1year, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.SSS.crypto$fit, aes(mean_sss_1year, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.SSS.crypto$fit, aes(x=mean_sss_1year, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.SSS.chondri$fit, aes(mean_sss_1year, visregFit), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.SSS.chondri$fit, aes(x=mean_sss_1year, ymin=visregLwr, ymax=visregUpr), fill='orange', alpha=.2)+
  geom_line(data=fit.SSS.large$fit, aes(mean_sss_1year, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.SSS.large$fit, aes(x=mean_sss_1year, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  annotate(geom="text", x=31.25, y=4.2, label="all MOTUs", hjust=0, size=4, color="darkred", fontface = "bold")+
  annotate(geom="text", x=31.25, y=4, label="Chondri", hjust=0, size=4, color="orange", fontface = "bold")+
  annotate(geom="text", x=31.25, y=3.8, label="Crypto", hjust=0, size=4, color="navy", fontface = "bold")+
  annotate(geom="text", x=31.25, y=3.6, label="Large fish", hjust=0, size=4, color="darkgreen", fontface = "bold")+
  labs(x="mean SSS", y="log(MOTUs richness)")+
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
load("Rdata/fit.MED.chondri.rdata")
load("Rdata/fit.MED.large.rdata")
load("Rdata/fit.MED.crypto.rdata")

ggplot() + 
  geom_line(data=fit.MED.motus$fit, aes(MarineEcosystemDependency, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.MED.motus$fit, aes(x=MarineEcosystemDependency, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.MED.crypto$fit, aes(MarineEcosystemDependency, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.MED.crypto$fit, aes(x=MarineEcosystemDependency, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.MED.chondri$fit, aes(MarineEcosystemDependency, visregFit), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.MED.chondri$fit, aes(x=MarineEcosystemDependency, ymin=visregLwr, ymax=visregUpr), fill='orange', alpha=.2)+
  geom_line(data=fit.MED.large$fit, aes(MarineEcosystemDependency, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.MED.large$fit, aes(x=MarineEcosystemDependency, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  annotate(geom="text", x=0.4, y=4.2, label="all MOTUs", hjust=1, size=4, color="darkred", fontface = "bold")+
  annotate(geom="text", x=0.4, y=4, label="Chondri", hjust=1, size=4, color="orange", fontface = "bold")+
  annotate(geom="text", x=0.4, y=3.8, label="Crypto", hjust=1, size=4, color="navy", fontface = "bold")+
  annotate(geom="text", x=0.4, y=3.6, label="Large fish", hjust=1, size=4, color="darkgreen", fontface = "bold")+
  labs(x="Marine Ecosystem Dependency", y="log(MOTUs richness)")+
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
load("Rdata/fit.CT.chondri.rdata")
load("Rdata/fit.CT.large.rdata")
load("Rdata/fit.CT.crypto.rdata")

ggplot() + 
  geom_line(data=fit.CT.motus$fit, aes(dist_to_CT, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.CT.motus$fit, aes(x=dist_to_CT, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.CT.crypto$fit, aes(dist_to_CT, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.CT.crypto$fit, aes(x=dist_to_CT, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.CT.chondri$fit, aes(dist_to_CT, visregFit), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.CT.chondri$fit, aes(x=dist_to_CT, ymin=visregLwr, ymax=visregUpr), fill='orange', alpha=.2)+
  geom_line(data=fit.CT.large$fit, aes(dist_to_CT, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.CT.large$fit, aes(x=dist_to_CT, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  annotate(geom="text", x=9.6, y=5.2, label="all MOTUs", hjust=1, size=4, color="darkred", fontface = "bold")+
  annotate(geom="text", x=9.6, y=5, label="Chondri", hjust=1, size=4, color="orange", fontface = "bold")+
  annotate(geom="text", x=9.6, y=4.8, label="Crypto", hjust=1, size=4, color="navy", fontface = "bold")+
  annotate(geom="text", x=9.6, y=4.6, label="Large fish", hjust=1, size=4, color="darkgreen", fontface = "bold")+
  labs(x="log(Distance to Coral Triangle (km))", y="log(MOTUs richness)")+
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
load("Rdata/fit.coast.chondri.rdata")
load("Rdata/fit.coast.large.rdata")
load("Rdata/fit.coast.crypto.rdata")

ggplot() + 
  geom_line(data=fit.coast.motus$fit, aes(distCoast, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.coast.motus$fit, aes(x=distCoast, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.coast.crypto$fit, aes(distCoast, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.coast.crypto$fit, aes(x=distCoast, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.coast.chondri$fit, aes(distCoast, visregFit), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.coast.chondri$fit, aes(x=distCoast, ymin=visregLwr, ymax=visregUpr), fill='orange', alpha=.2)+
  geom_line(data=fit.coast.large$fit, aes(distCoast, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.coast.large$fit, aes(x=distCoast, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  annotate(geom="text", x=13, y=4.5, label="all MOTUs", hjust=1, size=4, color="darkred", fontface = "bold")+
  annotate(geom="text", x=13, y=4.3, label="Chondri", hjust=1, size=4, color="orange", fontface = "bold")+
  annotate(geom="text", x=13, y=4.1, label="Crypto", hjust=1, size=4, color="navy", fontface = "bold")+
  annotate(geom="text", x=13, y=3.9, label="Large fish", hjust=1, size=4, color="darkgreen", fontface = "bold")+
  labs(x="log(Distance to coast (m))", y="log(MOTUs richness)")+
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
load("Rdata/fit.vol.chondri.rdata")
load("Rdata/fit.vol.large.rdata")
load("Rdata/fit.vol.crypto.rdata")

ggplot() + 
  geom_line(data=fit.vol.motus$fit, aes(volume, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.vol.motus$fit, aes(x=volume, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.vol.crypto$fit, aes(volume, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.vol.crypto$fit, aes(x=volume, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.vol.chondri$fit, aes(volume, visregFit), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.vol.chondri$fit, aes(x=volume, ymin=visregLwr, ymax=visregUpr), fill='orange', alpha=.2)+
  geom_line(data=fit.vol.large$fit, aes(volume, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.vol.large$fit, aes(x=volume, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  annotate(geom="text", x=1.6, y=4.5, label="all MOTUs", hjust=0, size=4, color="darkred", fontface = "bold")+
  annotate(geom="text", x=1.6, y=4.3, label="Chondri", hjust=0, size=4, color="orange", fontface = "bold")+
  annotate(geom="text", x=1.6, y=4.1, label="Crypto", hjust=0, size=4, color="navy", fontface = "bold")+
  annotate(geom="text", x=1.6, y=3.9, label="Large fish", hjust=0, size=4, color="darkgreen", fontface = "bold")+
  labs(x="log(volume (L))", y="log(MOTUs richness)")+
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
load("Rdata/fit.samp.chondri.rdata")
load("Rdata/fit.samp.large.rdata")
load("Rdata/fit.samp.crypto.rdata")

ggplot() + 
  geom_line(data=fit.samp.motus$fit, aes(depth_sampling, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.samp.motus$fit, aes(x=depth_sampling, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.samp.crypto$fit, aes(depth_sampling, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.samp.crypto$fit, aes(x=depth_sampling, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.samp.chondri$fit, aes(depth_sampling, visregFit), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.samp.chondri$fit, aes(x=depth_sampling, ymin=visregLwr, ymax=visregUpr), fill='orange', alpha=.2)+
  geom_line(data=fit.samp.large$fit, aes(depth_sampling, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.samp.large$fit, aes(x=depth_sampling, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  annotate(geom="text", x=1.6, y=4.5, label="all MOTUs", hjust=0, size=4, color="darkred", fontface = "bold")+
  annotate(geom="text", x=1.6, y=4.3, label="Chondri", hjust=0, size=4, color="orange", fontface = "bold")+
  annotate(geom="text", x=1.6, y=4.1, label="Crypto", hjust=0, size=4, color="navy", fontface = "bold")+
  annotate(geom="text", x=1.6, y=3.9, label="Large fish", hjust=0, size=4, color="darkgreen", fontface = "bold")+
  labs(x="depth of sampling (m)", y="log(MOTUs richness)")+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

ggsave("outputs/GLS/visreg_depthsampling.png", height = 5, width = 5)  


#### DHW5 ####

load("Rdata/fit.DHW.motus.rdata")
load("Rdata/fit.DHW.chondri.rdata")
load("Rdata/fit.DHW.large.rdata")
load("Rdata/fit.DHW.crypto.rdata")

ggplot() + 
  geom_line(data=fit.DHW.motus$fit, aes(mean_DHW_5year, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.DHW.motus$fit, aes(x=mean_DHW_5year, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.DHW.crypto$fit, aes(mean_DHW_5year, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.DHW.crypto$fit, aes(x=mean_DHW_5year, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.DHW.chondri$fit, aes(mean_DHW_5year, visregFit), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.DHW.chondri$fit, aes(x=mean_DHW_5year, ymin=visregLwr, ymax=visregUpr), fill='orange', alpha=.2)+
  geom_line(data=fit.DHW.large$fit, aes(mean_DHW_5year, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.DHW.large$fit, aes(x=mean_DHW_5year, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  annotate(geom="text", x=0, y=4.5, label="all MOTUs", hjust=0, size=4, color="darkred", fontface = "bold")+
  annotate(geom="text", x=0, y=4.3, label="Chondri", hjust=0, size=4, color="orange", fontface = "bold")+
  annotate(geom="text", x=0, y=4.1, label="Crypto", hjust=0, size=4, color="navy", fontface = "bold")+
  annotate(geom="text", x=0, y=3.9, label="Large fish", hjust=0, size=4, color="darkgreen", fontface = "bold")+
  labs(x="log(mean DHW 5year)", y="log(MOTUs richness)")+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

ggsave("outputs/GLS/visreg_DHW5.png", height = 5, width = 5)  
