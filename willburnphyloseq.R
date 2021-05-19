# this is written by Will Burn wlb501@york.ac.uk or wlwhalley@gmail.com
################################# THI BIT OF CODE INSTALLS PACKAGES #######################
# the analysis requires dada2, phyloseq, DECIPHER, remotes, speedyseq, taRifix, and vegan as essentials
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
library(phyloseq)
BiocManager::install("DECIPHER")
library(DECIPHER)
install.packages("remotes")
remotes::install_github("mikemc/speedyseq")
#get all the normal packages installed - delete the ones you don't need/use/want
packages <- c("plyr", "Hmisc", "ggplot2", "dplyr", "tidyr", "vegan", "DT", "reshape2", "knitr", "lubridate", "pwr", "psy", "car", "doBy", "corrplot", "RcmdrMisc", "questionr", "vcd", "multcomp", "KappaGUI", "rcompanion", "gridExtra", "factoextra", "corrplot", "FSA", "MASS", "scales", "nlme", "psych", "ordinal", "lmtest", "ggpubr", "dslabs", "stringr", "assist", "ggstatsplot", "forcats", "styler", "remedy", "snakecaser", "addinslist", "esquisse", "here", "summarytools", "magrittr", "tidyverse", "funModeling", "pander", "cluster", "abind")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])} # checks whether the packages in "packages" are installed and if not, installs them - only really works on CRAN packages, the ones on bioconductor or github normally need doing individually
# load em packages
invisible(lapply(packages, library, character.only = TRUE))
#import the data

################################## THIS BIT MAKES A PHYLOSEQ OBJECT ########################
seqtab<-readRDS("C:/VMshared/Bacteria/seqtab_bac.rds") # outputs from dada2 OR IT WILL TAKE .CSV FILES OF THESE THINGS
tax<-readRDS("C:/VMshared/Bacteria/tax_bac.rds") # outputs from dada2 OR A CSV OF YOUR TAXA LIST
metadata <- read.csv("C:/VMshared/Bacteria/master/metadata.csv", row.names=1)
tax<-as.matrix(tax)  ## <- this is required
# make a phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
                sample_data(metadata), 
                tax_table(tax))
ps


################################## THIS BIT CLUSTERS BY SEQUENCES SIMILARITY FOR 97% CLUSTERED OTUS ################
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps))) # this bit here changes your tax names to "ASV" rather than the sequences put out by dada2
dna <- refseq(ps)
nproc <- 2 # Increase to use multiple processors if you have a fancy PC - or are using the cluster!
aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
clusters <- DECIPHER::IdClusters(d, method = "complete",cutoff = 0.03, # corresponds to 97% OTUs
                                 processors = nproc)
ps2 <- speedyseq::merge_taxa_vec(ps, group = clusters$cluster, tax_adjust = 2)
ntaxa(ps) # see the OTUs before clustering
ntaxa(ps2) # see the OTUs post clustering

################################# THIS IS WHERE I FILTER MY TAXA ############################
ps3 <- prune_samples(sample_sums(ps2) >40, ps2) #Removes samples with less than 40 reads, check reads using sample_sums(ps2)
ps4 <- subset_taxa(ps3, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) # remove ASVs not determined to bacteria at PHYLUM level. You can also use this command to remove specific taxa groups just replace "uncharacterised"
ps5<-rarefy_even_depth(ps4, sample.size = min(sample_sums(ps4)),
                       rngseed = 42, replace = FALSE, trimOTUs = TRUE, verbose = TRUE) # normalise the data to even the sequencing length
write.csv(otu_table(ps5), "C:/VMshared/bacteria/tax4fun/rawotu.csv") # writes an OTU table of your new phyloseq
write.csv(tax_table(ps5), "C:/VMshared/bacteria/tax4fun/tax.csv") # writes a list of taxa for your new phyloseq

######## other filtering things for phyloseq
# lots more here: https://joey711.github.io/phyloseq/preprocess.html
# and some more here: https://ucdavis-bioinformatics-training.github.io/2017-September-Microbial-Community-Analysis-Workshop/friday/MCA_Workshop_R/phyloseq.html
# the second link there has code for filtering out low abundance taxa
species <- subset_taxa(bac, !is.na(Species) & !Species %in% c("", "NA")) # remove specific species if you want, or remove OTUs not identified to species by leaving it blank

################################ Transforming pyloseq ########################################
ps6<-transform_sample_counts(ps5, function(OTU) OTU/sum(OTU) ) #Relative abundance transformation
metadata<-as.data.frame(as.matrix(sample_data(ps5))) # edit metadata file to match newest version of phyloseq - important if you have removed low read count samples
## bray-curtis dissimilarity matrix
dist<-phyloseq::distance(ps5, method = "bray") # bray curtis distance matrix for later analysis in e.g PERMANOVA
dist2<-phyloseq::distance(ps5, method="wunifrac") # same as above but weighted unifrac (phylogenetically imformed) distance
## calculate diversity scores
e_alpha <- estimate_richness(ps4) # gives alha diversity scores for every measuure you need. Must be calculated on *INTEGERS therefore dont use abundance transformed data!!
write.csv(e_alpha, "yourfilestring.csv") # write a csv of your new diversity measures


################################ EDITING METADATA FOR ANALYSIS ###############################
# frequently the metadata is imported in stupid classes especially ecological data
#check the data
str(metadata)
# normally my data has been imported as a character rather than a numeric because of missing/na data. Also need to convert some characters
# to factors
metadata$treatment <- as.factor(metadata$treatment) # converts a column to a factor
metadata$soiltemp <- as.numeric(as.character(metadata$soiltemp)) # converts a column to numeric
library(taRifx)
metadata <- japply(metadata, which(sapply(metadata, class)=="character"), as.numeric) # takes any column that is character and converts it to numeric

#### DO A PERMANOVA
dist<-phyloseq::distance(ps5, method = "bray") #taxonomy based distance
# Test beta dispersal between factors
#This can also be thought of as a test for 'normality' between treatments
disp<-betadisper(dist, group = metadata$yourtreatmenthere)
disp <- permutest(disp, permutations = 999, pairwise = FALSE)
disp # PERMANOVA expects a non-significant p value here, if it is significant, within-treatment community differences might be contributing to a significant permanova
# however, if you have a balanced study design, this is not so much a problem 
# see Anderson & Walsh 2013 http://doi.org/10.1890/12-2010.1
# PERMANOVA analysis
vegan::adonis(dist ~ metadata$yourtreatmenthere, permutation = 9999, na.rm=TRUE) #Block Effect see -Bray-Curtis- taxonomic distance
# pairwise permanovas
install.packages("remotes")
remotes::install_github("leffj/mctoolsr")
library(mctoolsr)
metadata$natcode <- as.factor(metadata$natcode)
mctoolsr::calc_pairwise_permanovas(dist, metadata, "yourtreatmenthere", n_perm = 9999)

