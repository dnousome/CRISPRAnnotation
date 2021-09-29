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
library(rslurm)


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


##Prepare empty lists
out=list()



for (i in 1:length(out_tab)){
  #out_temp[[i]] <- mclapply(out_tab[[i]]$Aligned_Sequence, align_crispresso_p1, mc.cores = numCores)
  
  sopt1 <- list(time = '08:00:00',mem='24g')
  sjob <- slurm_apply(align_crispresso, out_tab[[i]], jobname = sprintf("%s_slurm",opt$out),
                      nodes = 12, cpus_per_node = 6, slurm_options=sopt1,global_objects = c("gene_sequence","gene_coords"),
                      submit = TRUE)
  out[[i]] <- get_slurm_out(sjob, outtype = 'raw', wait = TRUE)
  
  
}

  






names(out)=names(out_tab)
##Save output temporarily 
saveRDS(out,sprintf("temp_%s.rds",opt$out))




##For annovar Keep on the REF/ALT
vt=lapply(out,function(x){
  bind_rows(x) %>% arrange(Start) %>% distinct() 
})

vt_annovar=lapply(out,function(x){
  bind_rows(x) %>% arrange(Start) %>%  dplyr::select(chr,Start,End,REF,ALT) %>%
    distinct() 
  
})

##Load this VT_annovar object into ANNOVAR
lapply(names(vt_annovar),function(x){
  write_tsv(vt_annovar[[x]],paste0(x,"_toanno.avinput"),col_names = F)
})



###Run ANNOVAR output
##Use 2_runanno.sh
avin=paste0(names(out),"_toanno.avinput")
system("chmod 755 2_runannovar.sh")
lapply(avin,function(x)system(sprintf("./2_runannovar.sh --annovarin %s",x)))



#####Parse after output
annos=paste0("annovar/",names(out),".hg38_multianno.txt")
#annos=annos[order(sapply(strsplit(annos,"[/_]"),'[',2))]
#annos=annos[order(match(names(myfiles),names(out)))]


##Prep for Table output
afch=lapply(out,function(x){
  s=lapply(x,function(y)y[1,])
  bind_rows(s) %>% 
    dplyr::select(Aligned,Reference,n_deleted,n_inserted,n_mutated,Reads_n,Reads_prop) %>%
    mutate(af.id=1:nrow(.)) %>%
    #mutate(Reference=)
    dplyr::select(af.id,everything())
})


vt_full=mapply(function(x,y){
  anno_in=read_tsv(y,guess_max=20000)  %>%
    distinct() %>%
    dplyr::select(Chr,Start,End,Ref,Alt,Func.refGene,Gene.refGene,AAChange.refGene,
                  gnomad=`AF...132`,gnomad_non_topmed=non_topmed_AF_popmax,
                  gnomad_female=`AF_female...135`,gnomad_noncancer=non_cancer_AF_popmax,
                  SIFT_score,Polyphen2_HDIV_score,
                  Polyphen2_HVAR_score,CADD_raw=`CADD_raw...105`,phyloP30way_mammalian,
                  CLNSIG)
  
  names(x)=1:length(x)
  fin_dt=lapply(names(x),function(x1){
    if(!is.null(x[[x1]])){
      x[[x1]] %>% 
        mutate(af.id=as.numeric(x1)) 
    }
  })
  bind_rows(fin_dt) %>% 
    left_join(.,anno_in,by=c('chr'="Chr",'Start','End','REF'="Ref",'ALT'="Alt")) %>%
    dplyr::select(-Aligned,-Reference,-n_deleted,-n_inserted,-n_mutated,-Reads_n,-Reads_prop) %>%
    relocate(af.id,.after=last_col())
  
  
},out,annos,SIMPLIFY = F)

saveRDS(afch,sprintf("afch_%s.rds",opt$out))
saveRDS(vt_full,sprintf("vt_full_%s.rds",opt$out))



