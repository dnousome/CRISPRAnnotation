###Read in annotated ANNOVAR!
annos=Sys.glob("out/*")
annos=annos[order(sapply(strsplit(annos,"[/_]"),'[',2))]

anno_out=mapply(function(x,y){
  anno_in=read_tsv(x)  %>%
    dplyr::select(Chr,Start,End,Ref,Alt,Func.refGene,Gene.refGene,SIFT_score,
                  Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,CLNSIG)
  
  
  y %>% #mutate(chr=as.character(chr)) %>%
    left_join(.,anno_in,by=c('chr'="Chr",'Start','End','REF'="Ref",'ALT'="Alt")) %>%
    filter(Func.refGene=="exonic") %>% 
    dplyr::select(-Func.refGene)
  
  
  
},annos,vt,SIMPLIFY=F)


saveRDS(anno_out,"out.anno.RDS")

