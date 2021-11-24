library(ggplot2)
library(ggpubr)

load("Rdata/dbrda_MOTUs_province.rdata")
load("Rdata/dbrda_MOTUs_MED.rdata")
load("Rdata/dbrda_MOTUs_SST.rdata")
load("Rdata/dbrda_MNTD_province.rdata")
load("Rdata/dbrda_MNTD_MED.rdata")
load("Rdata/dbrda_MNTD_SST.rdata")


dbrda_MOTUs_province <- dbrda_MOTUs_province +
  ggtitle("A. Jaccard MOTUs composition")+
  theme(plot.title.position = "panel",
        plot.title = element_text(size = 12, color = "black", face = "bold"))

dbrda_province <- ggarrange(dbrda_MOTUs_province, dbrda_MNTD_province, 
                            common.legend = T, legend = c("bottom"), 
                            labels = c("", "D. MNTD"), font.label = list(size = 12, color = "black"),
                            label.x = 0)

dbrda_variable <- ggarrange(dbrda_MOTUs_SST, dbrda_MNTD_SST, dbrda_MOTUs_MED, dbrda_MNTD_MED, nrow=2, ncol=2, 
          labels = c("B", "E", "C", "F"), 
          font.label = list(size = 12, color = "black"),
          label.x = 0)
ggarrange(dbrda_province, dbrda_variable, nrow=2, ncol=1, heights = c(1,2))

ggsave("outputs/dbRDA/dbrda_jaccard_MNTD.png", width = 10, height = 12)
