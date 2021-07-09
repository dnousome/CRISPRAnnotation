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
  make_option(c("-i", "--inputpath"), type="character", default=NULL, 
              help="Path to All CRISPRESSOFiles with Sample Directories/Alleles_frequency_table.zip", metavar="character"),
  make_option(c("-g", "--gene"), type="character", default="BRCA2", 
              help="Gene name to grab coordinates for", metavar="character"),
  make_option(c("-p", "--pamsite"), type="character", default=NULL, 
              help="PAMSites: If more than one add comma to indicate [default= %default]", metavar="character"),
  make_option(c("-s", "--start"), type="integer", default=NULL, 
              help="Start Filter Range", metavar="integer"),
  make_option(c("-e", "--end"), type="integer", default=NULL, 
              help="End Filter Range", metavar="integer"),
  make_option(c("-o", "--out"), type="character", default="output", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$inputpath)){
  print_help(opt_parser)
  stop("Path Must be specified", call.=FALSE)
}


#setwd=("~/Documents/projects/ccbr1046/combined_flowcells/")

#batch_dir="~/Documents/projects/ccbr1046/combined_flowcells/CRISPRessoBatch_on_batch"

file_dir=opt$inputpath
crispresso_file_name="Alleles_frequency_table.zip"
myfiles = list.files(path=file_dir, pattern = crispresso_file_name, recursive = T, full.names = T)
names(myfiles) = gsub("^CRISPResso_on_","",basename(dirname(myfiles)))


##Load the Gene Name
if(file.exists(sprintf("%s_gene_coords.rds",opt$gene)) & file.exists(sprintf("%s_gene_sequence.rds",opt$gene))){
  gene_coords=readRDS(sprintf("%s_gene_coords.rds",opt$gene))
  gene_sequence=readRDS(sprintf("%s_gene_sequence.rds",opt$gene))
  
  }else{
mart = useMart('ensembl', dataset="hsapiens_gene_ensembl")
#mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl",host = "asia.ensembl.org")


my_attrs = c(gene="external_gene_name",chr="chromosome_name",start="start_position",
             end="end_position",strand="strand",sequence="gene_exon_intron")

gene_sequence_info = getBM(attributes = my_attrs,
                           filters = "external_gene_name", 
                           values = opt$gene, mart = mart, verbose=T)


gene_sequence_info = gene_sequence_info[,match(my_attrs, colnames(gene_sequence_info))]
colnames(gene_sequence_info) = names(my_attrs)


#colnames(gene_sequence_info)
#nchar(gene_sequence_info$sequence)



gene_coords = gene_sequence_info[gene_sequence_info$gene==opt$gene, 
                                 c("chr","start","end","strand")]
gene_sequence = gene_sequence_info$sequence[gene_sequence_info$gene==opt$gene]

saveRDS(gene_coords,sprintf("%s_gene_coords.rds",opt$gene))
saveRDS(gene_sequence,sprintf("%s_gene_sequence.rds",opt$gene))
}

##Read in the data
##Source the functions
#source("0_Functions.R")
source("0_FXN_update.R")

out_tab<-lapply(myfiles,allele_freq_tab)

out=list()


numCores=as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
#numCores=detectCores()-2

for (i in 1:length(out_tab)){
  out[[i]] <- mclapply(out_tab[[i]], align_crispresso, mc.cores = numCores)
}


names(out)=names(myfiles)
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
annos=paste0(names(myfiles),"_toanno.avinput.hg38_multianno.txt")
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


