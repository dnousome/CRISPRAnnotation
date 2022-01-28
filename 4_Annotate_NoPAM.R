########################
#
#Exon 20
#
##################################

setwd("/data/CCBR/rawdata/nousome/1046/Exon20/x20_annotate/CRISPRAnnotation")
source("0_FXN_update.R")
pacman::p_load("biomaRt")
getGeneInfo("BRCA2")

d1=Sys.glob("../CRISPResso*")
dout=sapply(strsplit(basename(d1),"_"),'[',6)
###Add dupliate statement

dout=make.unique(dout)


c1=sprintf('./CRISPR_Annotate.sh -i "%s/%s" -g BRCA2 -s 0 -e 0 -p 0 -a A -o %s',
           dirname(getwd()),basename(d1),dout)

write.table(c1,"all.sh",row.names=F,quote=F,col.names=F)
sapply(c1,system)



#####After annotation rerun with the updated codes
pacman::p_load(tidyverse,biomaRt,Biostrings,GenomicRanges,httr,readxl)
source("0_FXN_update.R")
opt=list()
opt$gene="BRCA2"

files=sub(".rds","",sub("afch_","",Sys.glob("afch*.rds")))
getGeneInfo(opt$gene)

##Add filter to the location to 
####Process all inital reads to calculate the N
process_output=function(filein){
  print(filein)
  af_file=readRDS(sprintf("afch_%s.rds",filein))
  vt_file=readRDS(sprintf("vt_full_%s.rds",filein)) 
  
  newdt=left_join(af_file[[1]],vt_file[[1]],by="af.id") %>%
    group_by(af.id) %>%
    
    mutate(true_n_del=max(str_count(Aligned,"-"),sum(ALT=="-")),
           true_n_ins=max(str_count(Reference,"-"),sum(REF=="-")),
           true_n_mutated=sum(ALT %in% c("A","T","G","C","N"))) %>% 
    ungroup() %>%
    relocate(Aligned,.after=ALT) %>%
    relocate(Reference,.after=Aligned) %>% distinct()
  
  ##Function for creating new
  out=list()
  newdt1=newdt %>% 
    filter(n_mutated==1 & n_deleted==0 & n_inserted==0)
  out[[1]]=newdt1 %>%
    dplyr::select(Reference,Aligned,AAChange.refGene,Reads_n,Reads_prop,REF,ALT,everything())
  newdt1=newdt %>% 
    filter(true_n_mutated==1 & true_n_del==0 & true_n_ins==0)
  out[[2]]=newdt1 %>% 
    dplyr::select(Reference,Aligned,AAChange.refGene,Reads_n,Reads_prop,REF,ALT,everything())
  return(out)
}

out_final=lapply(files,process_output)

names(out_final)=files
#names(out_final1)=files

openxlsx::write.xlsx(sapply(out_final,'[',2),sprintf("Annotated_filtered_%s_%s.xlsx",basename(dirname(getwd())),"Exon20"),overwrite=T)





##Try to split PAIRS by CIS/HAT,M15,OLA if they exists
pairs=lapply(c("(?i)CIS","(?i)HAT","(?i)M15","(?i)OLA"),function(x){
  grep(x,names(out_final),value=T)
})
names(pairs)=c("CIS","HAT","M15","OLA")


##Pivot FXN
pivot_all=function(pairs,inputfile){
  
  lapply(pairs,function(x){
    all_in=lapply(x,function(z){
      read_xlsx(inputfile,sheet=z)
    })
    names(all_in)=x
    all_in %>% bind_rows(.,.id = "id" ) %>% dplyr::select(-af.id) %>% #distinct() %>%
      #select(-Reference,-Aligned,-AAChange.refGene) %>%
      pivot_longer(cols=c(Reads_n,Reads_prop)) %>% 
      pivot_wider(#id_cols = c(REF,ALT,chr,Start,End),
        names_from=c(id,name),values_from=value) %>%
      dplyr::select(REF,ALT,everything())
  })
}
inputfile=sprintf("Annotated_filtered_%s_%s.xlsx",basename(dirname(getwd())),"Exon20")
out1=pivot_all(pairs,inputfile)
names(out1)=names(pairs)


openxlsx::write.xlsx(out1,sprintf("Annotated_bygroup_%s_%s.xlsx",basename(dirname(getwd())),"Exon20"),overwrite=T)




##Final Combination
final_combined=lapply(excel_sheets(sprintf("Annotated_bygroup_%s_%s.xlsx",basename(dirname(getwd())),"Exon20")),function(x){
  read_xlsx(sprintf("Annotated_bygroup_%s_%s.xlsx",basename(dirname(getwd())),"Exon20"),x)
}) %>% purrr::reduce(full_join, by = c('REF','ALT','Reference','Aligned','AAChange.refGene','n_deleted','n_inserted','n_mutated','chr',
                                       'Start','End','Func.refGene','Gene.refGene','gnomad','gnomad_non_topmed','gnomad_female','gnomad_noncancer',
                                       'SIFT_score','Polyphen2_HDIV_score','Polyphen2_HVAR_score','CADD_raw','phyloP30way_mammalian','CLNSIG','true_n_del',
                                       'true_n_ins','true_n_mutated')) %>% arrange(Start,REF)

openxlsx::write.xlsx(final_combined,sprintf("Annotated_combined_%s_%s.xlsx",basename(dirname(getwd())),"Exon20"),overwrite=T)















##ORiginal
out_final=lapply(files,function(z){
  print(z)
  af_file=readRDS(sprintf("afch_%s.rds",z))
  vt_file=readRDS(sprintf("vt_full_%s.rds",z)) 
  
  ##Filter out the vt file since there are a 
  temp=vt_file[[1]] #%>% filter(chr==13,Start>=32376646,End<=32376816)
  
  newdt=left_join(af_file[[1]],temp,by="af.id") %>%
    group_by(af.id) %>%
    
    mutate(true_n_del=max(str_count(Aligned,"-"),sum(ALT=="-")),
           #true_n_ins=max(str_count(Reference,"-"),sum(REF=="-")),
           true_n_mutated=sum(ALT %in% c("A","T","G","C","N"))) %>% 
    ungroup() %>%
    relocate(Aligned,.after=ALT) %>%
    relocate(Reference,.after=Aligned) %>% distinct()
  
  
  ##Function for creating new ref seq
  newdt1=newdt %>% 
    filter(n_mutated==1 & n_deleted==0 & n_inserted==0)
  newdt1 %>% dplyr::select(Reference,Aligned,AAChange.refGene,Reads_n,Reads_prop
                           ,REF,ALT,everything())
  
})


out_final1=lapply(files,function(z){
  print(z)
  af_file=readRDS(sprintf("afch_%s.rds",z))
  vt_file=readRDS(sprintf("vt_full_%s.rds",z)) 
  
  temp=vt_file[[1]] %>% filter(chr==13,Start>=32376646,End<=32376816)
  
  newdt=left_join(af_file[[1]],temp,by="af.id") %>%
    group_by(af.id) %>%
    
    
    mutate(true_n_del=max(str_count(Aligned,"-"),sum(ALT=="-")),
           #  true_n_ins=max(str_count(Reference,"-"),sum(REF=="-")),
           true_n_mutated=sum(ALT %in% c("A","T","G","C","N"))) %>% 
    ungroup() %>%
    relocate(Aligned,.after=ALT) %>%
    relocate(Reference,.after=Aligned) %>% distinct()
  
  
  newdt1=newdt %>% 
    filter(true_n_mutated==1 & true_n_del==0)
  newdt1 %>% dplyr::select(Reference,Aligned,AAChange.refGene,Reads_n,Reads_prop
                           ,REF,ALT,everything())
})

sapply(out_final,nrow)
sapply(out_final1,nrow)


names(out_final)=files
names(out_final1)=files
openxlsx::write.xlsx(out_final,sprintf("Annotated_Exon21.xlsx"),overwrite=T)
openxlsx::write.xlsx(out_final1,sprintf("Annotated_filtered_%s.xlsx",opt$out),overwrite=T)


openxlsx::write.xlsx(out_final,sprintf("Annotated_Exon21_test.xlsx"),overwrite=T)
openxlsx::write.xlsx(out_final1,sprintf("Annotated_filtered_%s_test.xlsx","Exon21"),overwrite=T)





##Try to split PAIRS by CIS/HAT,M15,OLA if they exists
pairs=lapply(c("(?i)CIS","(?i)HAT","(?i)M15","(?i)OLA"),function(x){
  grep(x,names(out_final),value=T)
})
names(pairs)=c("CIS","HAT","M15","OLA")

out1=lapply(pairs[[3]],function(x){
  all_in=lapply(x,function(z){
    read_xlsx("Annotated_filtered_Exon21_test.xlsx",sheet=z,
              col_types=c('text','text','text','numeric','numeric',
                          'text','text',
                          'numeric','numeric','numeric','numeric','numeric','numeric','numeric',
                          'text','text','numeric','numeric','numeric',
                          'numeric','numeric','numeric','numeric','numeric','numeric',
                          'text','numeric','numeric'))
  }) 
  names(all_in)=x
  all_in2=lapply(all_in,function(x){
    filter(x,!is.na(REF),!is.na(ALT))
  })
  #all_in1=lapply(all_in,function(x)
  #  classes(x)class)
  all_in2 %>% bind_rows(.,.id = "id" ) %>% dplyr::select(-af.id) %>% #distinct() %>%
    #select(-Reference,-Aligned,-AAChange.refGene) %>%
    pivot_longer(cols=c(Reads_n,Reads_prop)) %>% 
    pivot_wider(#id_cols = c(REF,ALT,chr,Start,End),
      names_from=c(id,name),values_from=value) %>%
    dplyr::select(REF,ALT,everything()) #%>%
  #filter(if_all(ends_with("Reads_prop"),~.>0.002))
  
  
})
names(out1)=c("M15","HAT","Cis","Ola")
openxlsx::write.xlsx(out1,"Annotated_filtered_bygroup_Exon21_test.xlsx",overwrite=T)




