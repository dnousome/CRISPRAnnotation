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



##Load or get gene coords/seqence from Biomart
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
  out[[i]] <- mclapply(out_tab[[i]], align_crispresso, mc.cores = numCores)
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
avin=paste0(names(myfiles),"_toanno.avinput")
system("chmod 755 2_runannovar.sh")
lapply(avin,function(x)system(sprintf("./2_runannovar.sh --annovarin %s",x)))



#####Parse after output
annos=paste0("annovar/",names(myfiles),".hg38_multianno.txt")
#annos=annos[order(sapply(strsplit(annos,"[/_]"),'[',2))]
annos=annos[order(match(names(myfiles),names(out)))]


##ADD THE FLOSS
#flossies=read_csv("~/Downloads/whi_ENSG00000139618_2021_06_21_15_13_26.csv")
#dleft_join(anno_out[[1]],flossies,
#          by=c("chr"="Chrom",'Start'='Position','REF'='Reference','ALT'='Alternate'))

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
                  gnomad=AF,gnomad_non_topmed=non_topmed_AF_popmax,
                  gnomad_female=AF_female,gnomad_noncancer=non_cancer_AF_popmax,
                  SIFT_score,Polyphen2_HDIV_score,
                  Polyphen2_HVAR_score,CADD_raw,phyloP30way_mammalian,
                  CLNSIG)
  
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




af_file=readRDS(sprintf("afch_%s.rds",opt$out))
vt_file=readRDS(sprintf("vt_full_%s.rds",opt$out))


lapply(1:length(af_file),function(i){
  title= names(af_file)[i]

  startrange=as.numeric(opt$start)
  endrange=as.numeric(opt$end)
  
  ##Option for PAMSITE
  if (!is.null(opt$pamsite)){
    pamsites=unlist(strsplit(as.character(opt$pamsite),","))
    pamsitesalt=unlist(strsplit(as.character(opt$pamsiteallele),","))
  }
  
  newdt=left_join(af_file[[i]],vt_file[[i]],by="af.id") %>%
    mutate(Biallelic=case_when(
      nchar(ALT)==1 & ALT!="-" ~ "Biallelic",
      TRUE~"Indel/Multiallelic"))  %>%
    mutate(WithinRange=case_when(
      (Start <= startrange|Start >= endrange)~"OutsideRange",
      (End <= startrange|End >= endrange)~"OutsideRange",
      TRUE~ "WithinRange"))  %>%
    mutate(AlignError=case_when(
      grepl("-",Aligned) & (n_deleted==0)~"CRISPRessoError",
      TRUE~"NoError")) %>%
    group_by(af.id) %>%
    mutate(Reference=ifelse(duplicated(Reference),"",Reference)) %>%
    mutate(TotalPAMSites = sum((Start %in% pamsites[1] & ALT==pamsitesalt[1] & REF %in% c("A","T","G","C")),
                               (Start %in% pamsites[2] & ALT==pamsitesalt[2] & REF %in% c("A","T","G","C")))) %>%
    mutate(PAMsite = case_when(
      TotalPAMSites==length(pamsites) ~ paste(as.character(pamsites),collapse="|"),
      TotalPAMSites==1 & any(Start %in% pamsites[1] & ALT==pamsitesalt[1] & REF %in% c("A","T","G","C")) ~as.character(pamsites[1]),
      TotalPAMSites==1 & any(Start %in% pamsites[2] & ALT==pamsitesalt[2] & REF %in% c("A","T","G","C")) ~as.character(pamsites[2])
    )) %>%
    mutate(true_n_del=sum((ALT=="-")),
           true_n_ins=sum((REF=="-")), 
           true_n_mutated=sum((ALT %in% c("A","T","G","C","N"))) 
    ) %>%
    ungroup() %>%
    relocate(Aligned,.after=ALT) %>%
    relocate(Reference,.after=Aligned) %>% distinct()
  
  
  newdt1=filter(newdt,
                PAMsite %in% c(paste(as.character(pamsites),collapse="|"),pamsites) &
                  AlignError=='NoError' &   WithinRange=="WithinRange" &  Biallelic=="Biallelic" &
                  true_n_del==0 & true_n_ins==0 & Start!=startrange & Start!=endrange) %>% 
    filter((true_n_mutated==2 & TotalPAMSites== 1)|(true_n_mutated==3 & TotalPAMSites==2)) %>%
    filter(!Start %in% pamsites) 
  
  newdt2=newdt1 %>% 
    group_by(chr,Start,REF,ALT) %>%
    summarize(Reads_total=sum(Reads_n),Reads_prop_total=sum(Reads_prop))
  
  newdt2=newdt1 %>% dplyr::select(chr,Start,End,REF,ALT,Func.refGene:CLNSIG) %>% distinct() %>%
    left_join(newdt2,.,by=c('chr','Start','REF','ALT'))
  
  dt=list(AllReads=newdt,FinalReads=newdt1,SummarizedFinalReads=newdt2)
  openxlsx::write.xlsx(dt,sprintf("%s_annotated.xlsx",title),overwrite = T)
})