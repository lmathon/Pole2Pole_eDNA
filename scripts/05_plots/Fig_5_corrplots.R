library(ggpubr)

load("Rdata/plot_alpha_phylo_gen.rdata")
load("Rdata/plot_alpha_trait_gen.rdata")
load("Rdata/plot_beta_phylo_gen.rdata")
load("Rdata/plot_beta_trait_gen.rdata")


ggarrange(plot_alpha_phylo_gen, plot_alpha_trait_gen, plot_beta_phylo_gen, plot_beta_trait_gen,
          nrow=2, ncol=2, labels = c("a", "b", "c", "d"), font.label = list(size = 24))

ggsave("outputs/Figures_papier/Fig5.png", width = 13, height = 14, dpi = 600)
