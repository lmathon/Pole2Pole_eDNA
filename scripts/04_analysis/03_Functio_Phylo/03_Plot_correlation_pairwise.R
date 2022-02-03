library(ggplot2)
library(ggpubr)

load("Rdata/plot_phylo_gen.rdata")
load("Rdata/plot_functio_gen.rdata")


ggarrange(plot_phylo_gen, plot_functio_gen, ncol=2, labels = c("a", "b"))

ggsave("outputs/Fig_sup_correlation_gen_phylo_functio.png", width = 12, height = 6)
