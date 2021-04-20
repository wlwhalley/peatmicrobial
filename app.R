library(microViz)
library(phyloseq)

setwd("C:/VMshared/bacteria/shiny")
bac2 <- readRDS("bac2.rds")

#check sequencing depth
ps2 <- prune_samples(sample_sums(bac2) >1500, bac2) #Removes samples with less than 40 reads i.e. WF4
sample_sums(ps2)
ps3 <- subset_taxa(ps2, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) # remove ASVs not determined to bacteria at PHYLUM level
ps4<-rarefy_even_depth(ps3, sample.size = min(sample_sums(ps3)),
                       rngseed = 42, replace = FALSE, trimOTUs = TRUE, verbose = TRUE) # normalise the data to even the sequencing length - this is Phils choice because
metadata<-as.data.frame(as.matrix(sample_data(ps4))) 
#### put it in and open the app
phy <- tax_fix(ps4)
ord1 <- phy %>%
  tax_transform("identity", rank = "Genus") %>% 
  dist_calc("bray") %>% # bray curtis
  ord_calc() #
ord_explore(data = ord1, auto_caption = NA)

