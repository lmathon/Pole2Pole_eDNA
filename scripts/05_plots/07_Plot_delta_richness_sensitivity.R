library(ggpubr)
library(ggplot2)

load("Rdata/delta_richness_sensitivity.rdata")

color <- c("red", "black", "black", "red", "red", "red", "red", "red", "red",
           "red", "black", "black", "red", "red", "black", "red", "red", "red",
           "red", "black", "red", "red", "red", "red", "red", "red", "red",
           "red", "black", "red", "red", "red", "red", "red", "red", "red",
           "red", "black", "red", "red", "red", "red", "red", "red", "red",
           "red", "black", "red", "red", "red", "red", "red", "red", "red",
           "red", "black", "red", "red", "red", "red", "red", "red", "red",
           "red", "black", "black", "red", "red", "black", "red", "red", "red",
           "red", "black", "red", "red", "red", "red", "red", "red", "red",
           "red", "black", "black", "red", "red", "red", "red", "red", "red")

ggplot(data = delta_rich_fin, 
       aes(x = variable, y = delta)) +
  geom_hline(aes(yintercept = 0), colour = "black") + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), 
                colour = "black", size = 0.5, width = 0) +
  geom_point(size = 2, col=color) +
  coord_flip() + 
  ylim(-110, 150)+
  theme_sleek(base_size = 24) + 
  facet_grid(. ~ taxa, scales = "free_y", space = "free_y", switch = "y") + 
  #scale_y_continuous(breaks = c(-100,-50,0,50,100)) + 
  ylab(expression(paste(Delta,"% MOTU richness")))+
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12),
        strip.text.x = element_text(size=10,face="bold"),
        strip.text.y = element_text(size=10,face="bold"),
        strip.placement = "outside")

ggsave("outputs/GLS/delta_socio_sensitivity.png", width=6.5, height = 2)