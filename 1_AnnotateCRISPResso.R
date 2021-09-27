#!/usr/bin/env Rscript
###CRISPRESSO Alignment Project
###By Darryl Nousome
library(optparse)
library(tidyverse)
library(biomaRt)
library(Biostrings)
library(GenomicRanges)
library(httr)
library(parallel)

set_config(config(ssl_verifypeer = 0L))


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path or files with all CRISPRESSOFiles with Sample Directories/Alleles_frequency_table.zip", metavar="character"),
  make_option(c("-g", "--gene"), type="character", default="BRCA2", 
              help="Gene name to grab coordinates for", metavar="character"),
  make_option(c("-p", "--pamsite"), type="character", default=NULL, 
              help="PAMSites: If more than one add comma to indicate [default= %default]", metavar="character"),
  make_option(c("-a", "--pamsiteallele"), type="character", default=NULL, 
              help="PAMSitesAlleles: If more than one add comma to indicate [default= %default]", metavar="character"),
  make_option(c("-s", "--start"), type="integer", default=NULL, 
              help="Start Filter Range", metavar="integer"),
  make_option(c("-e", "--end"), type="integer", default=NULL, 
              help="End Filter Range", metavar="integer"),
  make_option(c("-o", "--out"), type="character", default="output", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Path or File must be specified", call.=FALSE)
}


#setwd=("~/Documents/projects/ccbr1046/combined_flowcells/")
#batch_dir="~/Documents/projects/ccbr1046/combined_flowcells/CRISPRessoBatch_on_batch"

##Check the file type to ensure that reading in works correctly
##Parse either Excel or Zipped files direct from CRISPRESSO

##Load the Gene Name and all in the function listed
##Read in the data
##Source the functions
source("0_FXN_update.R")


##Load or get gene coords/sequence from Biomart
getGeneInfo(opt$gene)

##Get File Extension
ext=tail(unlist(strsplit(opt$input,"\\.")),n=1)

##
out_tab=allele_freq_tab(opt$input,file_ext = ext)    




numCores=as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
if (is.na(numCores)){
  numCores=detectCores()-2
}



out=list()
for (i in 1:length(out_tab)){
  #out[[i]] <- mclapply(out_tab[[i]][1:20], align_crispresso, mc.cores = numCores)
  out[[i]] <- mclapply(out_tab[[i]], align_crispresso_p1, mc.cores = numCores)
}


names(out)=names(out_tab)
##Save output temporarily 
saveRDS(out,sprintf("temp_%s.rds",opt$out))
