library(ggplot2)
library(ggpubr)

load("Rdata/dbrda_MOTUs_province.rdata")
load("Rdata/dbrda_MOTUs_MED.rdata")
load("Rdata/dbrda_MOTUs_SST.rdata")
load("Rdata/dbrda_FD_province.rdata")
load("Rdata/dbrda_FD_MED.rdata")
load("Rdata/dbrda_FD_SST.rdata")


dbrda_MOTUs_province <- dbrda_MOTUs_prov +
  ggtitle(expression(paste("MOTU ", beta,"-diversity")))+
  labs(color="Region")+
  theme(plot.title.position = "panel",
        plot.title = element_text(size = 14, color = "black", face = "bold"),
        legend.text=element_text(size=11))

dbrda_FD_province <- dbrda_FD_prov +
  ggtitle(expression(paste("Sequence ", beta,"-diversity")))+  
  theme(plot.title.position = "panel",
        plot.title = element_text(size = 14, color = "black", face = "bold"))

dbrda_province <- ggarrange(dbrda_MOTUs_province, dbrda_FD_province, labels = c("(a)", "(b)"), 
                            common.legend = T, legend = c("bottom")) 
                            
ggsave(dbrda_province, file="outputs/dbRDA/dbrda_jaccard_FD.png", width = 10.7, height = 6.2)
ggsave(dbrda_province, file="outputs/Figures_papier/Fig4.pdf", width = 10.7, height = 6.2)
ggsave(dbrda_province, file="outputs/Figures_papier/Fig4.png", width = 10.7, height = 6.2)





dbrda_variable <- ggarrange(dbrda_MOTUs_SST, dbrda_FD_SST, dbrda_MOTUs_MED, dbrda_FD_MED, nrow=2, ncol=2, 
          labels = c("(a)", "(b)", "(c)", "(d)"), 
          font.label = list(size = 12, color = "black"),
          label.x = 0)

ggsave("outputs/dbRDA/dbrda_jaccard_FD_sup.png", width = 10, height = 8)
