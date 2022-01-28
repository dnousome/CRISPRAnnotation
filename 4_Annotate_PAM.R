###########
#
#
###############

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