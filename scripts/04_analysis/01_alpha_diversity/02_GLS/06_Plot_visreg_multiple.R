library(ggplot2)
library(pdp)
library(tidyverse)

#### Gravity ####

load("Rdata/fit.grav.motus.rdata")
load("Rdata/fit.grav.chondri.rdata")
load("Rdata/fit.grav.large.rdata")
load("Rdata/fit.grav.crypto.rdata")

ggplot() + 
  geom_line(data=fit.grav.motus$fit, aes(Gravity, visregFit), colour='navy', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.motus$fit, aes(x=Gravity, ymin=visregLwr, ymax=visregUpr), fill='navy', alpha=.2)+
  geom_line(data=fit.grav.crypto$fit, aes(Gravity, visregFit), colour='tomato', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.crypto$fit, aes(x=Gravity, ymin=visregLwr, ymax=visregUpr), fill='tomato', alpha=.2)+
  geom_line(data=fit.grav.chondri$fit, aes(Gravity, visregFit), colour='orange', size=1, show.legend = T)+
  geom_ribbon(data=fit.grav.chondri$fit, aes(x=Gravity, ymin=visregLwr, ymax=visregUpr), fill='orange', alpha=.2)+
  geom_line(data=fit.gav.large$fit, aes(Gravity, visregFit), colour='darkgreen', size=1, show.legend = T)+
  geom_ribbon(data=fit.gav.large$fit, aes(x=Gravity, ymin=visregLwr, ymax=visregUpr), fill='darkgreen', alpha=.2)+
  annotate(geom="text", x=7, y=4, label="all MOTUs", hjust=1, size=4, color="navy", fontface = "bold")+
  annotate(geom="text", x=7, y=3.85, label="Chondri", hjust=1, size=4, color="orange", fontface = "bold")+
  annotate(geom="text", x=7, y=3.7, label="Crypto", hjust=1, size=4, color="tomato", fontface = "bold")+
  annotate(geom="text", x=7, y=3.55, label="Large fish", hjust=1, size=4, color="darkgreen", fontface = "bold")+
  labs(x="Gravity", y="log(MOTUs richness)")+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

ggsave("outputs/GLS/visreg_gravity.png", height = 5, width = 5)  
