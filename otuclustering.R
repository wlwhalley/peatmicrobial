#####################################################
# this was written by Will Burn williamlburn@gmail.com or lb501@york.av.uk
#####################################################
### Clutering by OTU
# packages required
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "phyloseq", "DECIPHER")
install.packages("remotes")
remotes::install_github("mikemc/speedyseq")
Packages <- c("Biostrings", "phyloseq", "DECIPHER", "speedyseq")
lapply(Packages, library, character.only = TRUE)
##### import your data and make a phyloseq
#import the data
seqtab<-readRDS("C:/VMshared/Bacteria/seqtab_bac.rds") # if your phyloseq inputs are in csv format use read.csv instead
tax<-readRDS("C:/VMshared/Bacteria/tax_bac.rds")
metadata_bac <- read.csv("C:/VMshared/Bacteria/master/bmasterph.csv", row.names=1)
tax<-as.matrix(tax)  ## this is required
#View what makes up a phyloseq object
bac <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
                sample_data(metadata_bac), 
                tax_table(tax))
bac
#### now cluster by 97% sequence similarity
dna <- Biostrings::DNAStringSet(taxa_names(bac))
names(dna) <- taxa_names(bac)
bac <- merge_phyloseq(bac, dna)
taxa_names(bac) <- paste0("ASV", seq(ntaxa(bac))) # choose whether or not to use this line. It will replace the 
# sequence outputs from DADA2 with "ASV1" for species names. 
dna <- refseq(bac)
nproc <- 2 # Increase to use multiple processors, e.g if you have 4 or several hundred if youre using the cluster :) 
aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
clusters <- DECIPHER::IdClusters(d, method = "complete",cutoff = 0.03, # corresponds to 97% OTUs, change for different cut off
                                 processors = nproc)
bac2 <- speedyseq::merge_taxa_vec(bac, group = clusters$cluster, tax_adjust = 2)
# inspect the changes in number of species
ntaxa(bac) #7158
ntaxa(bac2) #2367
sample_sums(bac)