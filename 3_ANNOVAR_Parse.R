#########
#
#Read in annotated ANNOVAR
#
#######
annos=Sys.glob("out/*")
annos=annos[order(sapply(strsplit(annos,"[/_]"),'[',2))]

anno_out=mapply(function(x,y){
  anno_in=read_tsv(x)  %>%
    dplyr::select(Chr,Start,End,Ref,Alt,Func.refGene,Gene.refGene,AAChange.refGene,
                  SIFT_score,Polyphen2_HDIV_score,
                  Polyphen2_HVAR_score,CADD_raw,
                  CLNSIG)
  
  
  y %>% #mutate(chr=as.character(chr)) %>%
    left_join(.,anno_in,by=c('chr'="Chr",'Start','End','REF'="Ref",'ALT'="Alt"))# %>%
    #filter(Func.refGene=="exonic") %>% 
    #dplyr::select(-Func.refGene)
  
  
  
},annos[1],vt,SIMPLIFY=F)


saveRDS(anno_out,"out.anno.RDS")
saveRDS(out,"original_vt.RDS")
saveRDS(vt,"vt.rds")



##Prep for Shinytable
afch=lapply(out,function(x){
  s=lapply(x,function(y)y[1,])
  bind_rows(s) %>% dplyr::select(Aligned,Reference,n_deleted,n_inserted,n_mutated,Reads_n,Reads_prop) %>%
    mutate(af.id=1:nrow(.)) %>%
    dplyr::select(af.id,everything())
  })


vt_full=mapply(function(x,y){
  anno_in=read_tsv(y)  %>%
    dplyr::select(Chr,Start,End,Ref,Alt,Func.refGene,Gene.refGene,AAChange.refGene,
                  SIFT_score,Polyphen2_HDIV_score,
                  Polyphen2_HVAR_score,CADD_raw,
                  CLNSIG)
  
  fin_dt=lapply(names(x),function(x1){
    x[[x1]] %>% #mutate(chr=as.character(chr)) %>%
      left_join(.,anno_in,by=c('chr'="Chr",'Start','End','REF'="Ref",'ALT'="Alt")) %>%
      dplyr::select(-Aligned,-Reference,-n_deleted,-n_inserted,-n_mutated,-Reads_n,-Reads_prop) %>%
      mutate(af.id=as.numeric(x1)) 
  })
  bind_rows(fin_dt)

  },out,annos[1],SIMPLIFY = F)


saveRDS(afch,"afch.rds")
saveRDS(vt_full,"vt_full.rds")
