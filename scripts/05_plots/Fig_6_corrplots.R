library(ggpubr)

load("Rdata/plot_alpha_phylo_gen.rdata")
load("Rdata/plot_alpha_trait_gen.rdata")
load("Rdata/plot_beta_phylo_gen.rdata")
load("Rdata/plot_beta_trait_gen.rdata")

plot_alpha_trait_gen <- plot_alpha_trait_gen +
  ylim(0,22)+
  annotate(geom="text", x=1, y=21, label="Pearson cor=0.91 \n p<0.001", hjust=0, size=6, fontface = "bold")

ggarrange(plot_alpha_phylo_gen, plot_alpha_trait_gen, plot_beta_phylo_gen, plot_beta_trait_gen,
          nrow=2, ncol=2, labels = c("A", "B", "C", "D"), font.label = list(size = 24))

ggsave("outputs/Figures_papier/Fig6.pdf", width = 13, height = 14, dpi = 600)
ggsave("outputs/Figures_papier/Fig6.png", width = 13, height = 14, dpi = 600)
