library(dplyr)
library(rfishbase)

# get species taxo and length from fishbase
x <- load_taxa()
species <- collect(x)

length <- species(fields = c("Species", "Length"))

species_length <- left_join(species, length, by="Species")
species_length <- species_length %>%
  distinct(Species, .keep_all=T)

# get mean and 5-95 interval of size per family
family_length <- data.frame(Family=as.character(), mean_size=as.numeric(), decile5_size=as.numeric(), decile95_size=as.numeric())

family <- unique(species_length$Family)

for (i in 1:length(family)) {
  sub <- species_length %>%
    subset(Family==family[i]) %>%
    filter(!is.na(Length))
  mean <- mean(sub$Length)
  quant <- quantile(sub$Length, c(.95, .05))
  family_length[i, "Family"] <- family[i]
  family_length[i, "mean_size"] <- mean
  family_length[i, "decile5_size"] <- quant[2]
  family_length[i, "decile95_size"] <- quant[1]
}

large_families <- family_length %>%
  subset(decile5_size > 20) 
large_families <- as.character(large_families$Family)

save(large_families, file = "Rdata/large_families.rdata")

# get mean and 5-95 interval of size per order
order_length <- data.frame(Order=as.character(), mean_size=as.numeric(), decile5_size=as.numeric(), decile95_size=as.numeric())

order <- unique(species_length$Order)

for (i in 1:length(order)) {
  sub <- species_length %>%
    subset(Order==order[i]) %>%
    filter(!is.na(Length))
  mean <- mean(sub$Length)
  quant <- quantile(sub$Length, c(.95, .05))
  order_length[i, "Order"] <- order[i]
  order_length[i, "mean_size"] <- mean
  order_length[i, "decile5_size"] <- quant[2]
  order_length[i, "decile95_size"] <- quant[1]
}

large_orders <- order_length %>%
  subset(decile5_size > 20) 
large_orders <- as.character(large_orders$Order)

save(large_orders, file = "Rdata/large_orders.rdata")
