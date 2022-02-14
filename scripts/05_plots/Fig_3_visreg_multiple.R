library(ggplot2)
library(grid)
library(cowplot)
library(tidyverse)

#### Gravity ####

load("Rdata/fit.grav.motus.rdata")
load("Rdata/fit.grav.large.rdata")
load("Rdata/fit.grav.crypto.rdata")
load("Rdata/fit.grav.FDq0.rdata")


ylim.prim <- c(0, 2)   
ylim.sec <- c(0, 15) 

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

gravity <- ggplot() + 
  geom_line(data=fit.grav.motus$fit, aes(Gravity, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.motus$fit, aes(x=Gravity, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.grav.crypto$fit, aes(Gravity, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.crypto$fit, aes(x=Gravity, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.grav.FDq0$fit, aes(Gravity, (visregFit*b)+a), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.FDq0$fit, aes(x=Gravity, ymin=(visregLwr*b)+a, ymax=(visregUpr*b)+a), fill='orange', alpha=.2)+
  geom_line(data=fit.grav.large$fit, aes(Gravity, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.large$fit, aes(x=Gravity, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  labs(x="log10(Gravity+1)")+
  scale_y_continuous("log(taxonomic a-diversity)", sec.axis = sec_axis(~ (. - a)/b, name = "sequence a-diversity")) +
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
load("Rdata/fit.SST.FDq0.rdata")

sst <- ggplot() + 
  geom_line(data=fit.SST.motus$fit, aes(mean_SST_1year, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.SST.motus$fit, aes(x=mean_SST_1year, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.SST.crypto$fit, aes(mean_SST_1year, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.SST.crypto$fit, aes(x=mean_SST_1year, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.SST.FDq0$fit, aes(mean_SST_1year, (visregFit*b)+a), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.SST.FDq0$fit, aes(x=mean_SST_1year, ymin=(visregLwr*b)+a, ymax=(visregUpr*b)+a), fill='orange', alpha=.2)+
  geom_line(data=fit.SST.large$fit, aes(mean_SST_1year, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.SST.large$fit, aes(x=mean_SST_1year, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  annotate(geom="text", x=0, y=2, label="All fish", hjust=0, size=3, color="darkred", fontface = "bold")+
  annotate(geom="text", x=0, y=2.1, label="Sequence", hjust=0, size=3, color="orange", fontface = "bold")+
  annotate(geom="text", x=0, y=1.9, label="Cryptobenthic", hjust=0, size=3, color="navy", fontface = "bold")+
  annotate(geom="text", x=0, y=1.8, label="Large fish", hjust=0, size=3, color="darkgreen", fontface = "bold")+
  labs(x="mean SST")+
  scale_y_continuous("log(taxonomic a-diversity)", sec.axis = sec_axis(~ (. - a)/b, name = "sequence a-diversity")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

ggsave("outputs/GLS/visreg_SST.png", height = 5, width = 5)  


#### MED ####

load("Rdata/fit.MED.motus.rdata")
load("Rdata/fit.MED.large.rdata")
load("Rdata/fit.MED.crypto.rdata")
load("Rdata/fit.MED.FDq0.rdata")

med <- ggplot() + 
  geom_line(data=fit.MED.motus$fit, aes(MarineEcosystemDependency, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.MED.motus$fit, aes(x=MarineEcosystemDependency, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.MED.crypto$fit, aes(MarineEcosystemDependency, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.MED.crypto$fit, aes(x=MarineEcosystemDependency, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.MED.FDq0$fit, aes(MarineEcosystemDependency, (visregFit*b)+a), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.MED.FDq0$fit, aes(x=MarineEcosystemDependency, ymin=(visregLwr*b)+a, ymax=(visregUpr*b)+a), fill='orange', alpha=.2)+
  geom_line(data=fit.MED.large$fit, aes(MarineEcosystemDependency, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.MED.large$fit, aes(x=MarineEcosystemDependency, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  labs(x="Marine Ecosystem Dependency")+
  scale_y_continuous("log(taxonomic a-diversity)", sec.axis = sec_axis(~ (. - a)/b, name = "sequence a-diversity")) +
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
load("Rdata/fit.CT.FDq0.rdata")

distCT <- ggplot() + 
  geom_line(data=fit.CT.motus$fit, aes(dist_to_CT, visregFit), colour='darkred', size=1, show.legend = T)+
  geom_ribbon(data=fit.CT.motus$fit, aes(x=dist_to_CT, ymin=visregLwr, ymax=visregUpr), fill='darkred', alpha=.2)+
  geom_line(data=fit.CT.crypto$fit, aes(dist_to_CT, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.CT.crypto$fit, aes(x=dist_to_CT, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.CT.FDq0$fit, aes(dist_to_CT, (visregFit*b)+a), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.CT.FDq0$fit, aes(x=dist_to_CT, ymin=(visregLwr*b)+a, ymax=(visregUpr*b)+a), fill='orange', alpha=.2)+
  geom_line(data=fit.CT.large$fit, aes(dist_to_CT, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.CT.large$fit, aes(x=dist_to_CT, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  labs(x="log10(Distance to Coral Triangle (km) +1)")+
  scale_y_continuous("log(taxonomic a-diversity)", sec.axis = sec_axis(~ (. - a)/b, name = "sequence a-diversity")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

ggsave("outputs/GLS/visreg_distCT.png", height = 5, width = 5)  



#### plot all together ####

ggarrange(sst, distCT, gravity, med, nrow=2, ncol=2, labels = c("A", "B", "C", "D"))
ggsave("outputs/GLS/visreg_multiple_FDq0.png", width = 8, height = 7)
ggsave("outputs/Figures_papier/Fig3.pdf", width = 8, height = 7)


