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

metadata    = data.frame(read.table("metadata.txt",  sep = "\t", header=T, stringsAsFactors = F, check.names=FALSE))






#####################
## Format taxonomy ##
#####################

taxo_table = Genoscreen_Table_ASV[,c("ASV_ID", "Taxon")]
taxo_table$Taxon = gsub(" *.__", "", taxo_table$Taxon)

taxo_table_split = cbind(taxo_table$ASV_ID, data.frame(str_split_fixed(taxo_table$Taxon, pattern = ";", n = 7)))
colnames(taxo_table_split) = c("ASV_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Specie")





