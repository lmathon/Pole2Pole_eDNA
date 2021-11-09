library(ggplot2)
library(ggpubr)


load("Rdata/dbrda_MOTUs_MED.rdata")
load("Rdata/dbrda_MOTUs_SST.rdata")
load("Rdata/dbrda_MPD_MED.rdata")
load("Rdata/dbrda_MPD_SST.rdata")

dbrda_MOTUs_SST <- dbrda_MOTUs_SST +
  ggtitle("A. Jaccard MOTUs composition")+
  theme(plot.title.position = "plot",
        plot.title = element_text(size = 12, color = "black", face = "bold"))

ggarrange(dbrda_MOTUs_SST, dbrda_MPD_SST, dbrda_MOTUs_MED, dbrda_MPD_MED, nrow=2, ncol=2, 
          labels = c("", "C. MPD", "B", "D"), 
          font.label = list(size = 12, color = "black"),
          label.x = 0)


ggsave("outputs/dbRDA/dbrda_jaccard_MPD.png", width = 10, height = 11)
