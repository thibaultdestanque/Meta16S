#!/usr/bin/Rscript

session_info = utils::sessionInfo()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


###########################################
#######  Metagenomic with phyloseq  #######  
###########################################

# Clean environnement
ls()
rm(list=ls())
ls()


###############
##  Library  ##
###############

library(vroom)
library(plyr); library(dplyr) # plyr need to be loaded before dplyr
library(stringr)
library(phyloseq)
library(tibble)
library(vegan)
library(ggplot2)
library(tidyr)
library(devtools)
library(microbiome)
library(upstartr)
library("easystats")
library(kableExtra)
library("network")
library("hrbrthemes")
library("readxl")




###########################
## Set working directory ##
###########################

setwd("C:/Users/t.destanque/Documents/04_Metagenomic_16S")



###########################################
##  Import abondance table and taxonomy  ##
###########################################

Genoscreen_Table_ASV <- read_excel("Genoscreen_Table_ASV_97431.xlsx",skip = 1)
colnames(Genoscreen_Table_ASV)
colnames(Genoscreen_Table_ASV)[1] = "ASV_ID"



#####################
## Import metadata ##
#####################

metadata    = data.frame(read.table("metadata_methodo.txt",  sep = "\t", header=T, stringsAsFactors = F, check.names=FALSE))
head(metadata)


###########################
## Format taxonomy table ##
###########################

# import column of interest
taxo_table = Genoscreen_Table_ASV[,c("ASV_ID", "Taxon")]

# Remoove useless info : d__Bacteria => Bacteria;
taxo_table$Taxon = gsub(" *.__", "", taxo_table$Taxon)

taxo_table_split = cbind(taxo_table$ASV_ID, data.frame(str_split_fixed(taxo_table$Taxon, pattern = ";", n = 7)))
colnames(taxo_table_split) = c("ASV_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Specie")

taxo_table = taxo_table_split
head(taxo_table)


######################
## Format ASV table ##
######################

ASV_table = Genoscreen_Table_ASV
ASV_table$Taxon = NULL
colnames(ASV_table)

# On recupére les colonnes de notre projet (MARISA)
#ASV_table = ASV_table[,c("ASV_ID", "10-V3L1", "11-V13L2", "12-V6L3", "7-V4L1", "8-V2L2", "9-V8L3")]

# On renomme les noms de samples pour qu'ils matchent avec ceux du fichier metadata (on enleve le "[0-9]-")
#colnames(ASV_table) = gsub(".*-", "", colnames(ASV_table))
#head(ASV_table)

# On recupére les colonnes de notre projet (AGNESE)
ASV_table = ASV_table[,c("ASV_ID", "1-vortex", "2-vortex-centri", "3-BB3c", "4-BB3c-centri", "5-BB9c", "6-BB9c-centri")]

# On renomme les noms de samples pour qu'ils matchent avec ceux du fichier metadata (on enleve le "[0-9]-")
colnames(ASV_table) = gsub("[0-9]-", "", colnames(ASV_table))
#colnames(ASV_table) = gsub("-", "_", colnames(ASV_table))
head(ASV_table)


###########################
## Build phyloseq object ##
###########################

otu_mat <- ASV_table %>%
  tibble::column_to_rownames("ASV_ID") %>%
  as.matrix()
head(otu_mat)

tax_mat <- taxo_table %>% 
  tibble::column_to_rownames("ASV_ID") %>%
  as.matrix()
head(tax_mat)

samples_df <- metadata %>% 
  tibble::column_to_rownames("Sample") 
head(samples_df)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

# Transform to phyloseq objects
pseq <- phyloseq(OTU, TAX, samples)
pseq
pseq@otu_table

sample_names(pseq)
rank_names(pseq)
sample_variables(pseq)

pseq_raw = pseq

## TEST ##
# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(pseq))
standf = function(x, t=total) round(t * (x / sum(x)))
pseq_norm1 = transform_sample_counts(pseq, standf)
## END TEST ##

# Quelques infos utiles avant de se lancer dans les figures
ntaxa(pseq)
nsamples(pseq)
sample_names(pseq)[1:5] 
rank_names(pseq)  
sample_variables(pseq)  
otu_table(pseq)[1:5, 1:5]  
tax_table(pseq)[1:5, 1:4]


# Rarefy the phyloseq object to even depth prior various analysis
set.seed(1)
physeq_rarefy <- rarefy_even_depth(pseq, rngseed=1, sample.size=0.9*min(sample_sums(pseq)), replace=F)
plot_taxa_prevalence(pseq, "Phylum")
plot_taxa_prevalence(physeq_rarefy, "Phylum")



###########
# barplot #
###########

# Composition plot
plot_bar(physeq_rarefy, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")




###########
# Heatmap #
###########

plot_heatmap(pseq, method = "NMDS", distance = "bray")

pseq_abund <- filter_taxa(pseq, function(x) sum(x > total*0.01) > 0, TRUE)
pseq_abund

# Filter
plot_heatmap(pseq_abund, method = "NMDS", distance = "bray")

# Filter with taxa names
plot_heatmap(pseq_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Class", taxa.order = "Class")



##############
# Rare curve #
##############

tab <- otu_table(physeq_rarefy)
class(tab) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have
tab <- t(tab) # transpose observations to rows
raremax = min(rowSums(tab))
rare <- rarecurve(tab, step=100, lwd=2, ylab="Richness", sample = raremax, col = "blue", label=F)



###################
# Alpha diversity #
###################

plot_richness(pseq, measures=c("observed", "Chao1", "Shannon"))



##############
# Ordination #
##############

pseq.ord <- ordinate(pseq_abund, "NMDS", "bray")

plot_ordination(pseq, pseq.ord, type="taxa", color="Phylum", shape= "Phylum", 
                title="OTUs")

plot_ordination(physeq_rarefy, pseq.ord, type="taxa", color="Phylum", shape= "Phylum", 
                title="OTUs")
plot_ordination(physeq_rarefy, pseq.ord, type="taxa", color="Genus", shape= "Genus", 
                title="OTUs")













