########All Support Functions

#Process Allele frequency_table--Output as a dataframe with only the Aligned Sequences
allele_freq_tab=function(file,file_ext){
  
  ##Xlsx/txt will return list length 1
  if(file_ext=="xlsx"){
    sheets=readxl::excel_sheets(file)
    
    out_tab=lapply(sheets,function(x){
      mytable=readxl::read_xlsx(file,sheet = x)
      this <- mytable %>% 
        dplyr::select(Aligned_Sequence,
                      Reference_Sequence,
                      n_deleted,n_inserted,n_mutated,Reads_n=`#Reads`,Reads_prop=`%Reads`) %>%
        filter(n_mutated > 0) 
     # split(this,1:nrow(this))
      
    })
    names(out_tab)=sheets
    return(out_tab)
    
  }else if(file_ext=="txt"){
    mytable = read_tsv(x)
    this <- mytable %>% 
      dplyr::select(Aligned_Sequence,
                    Reference_Sequence,
                    n_deleted,n_inserted,n_mutated,Reads_n=`#Reads`,Reads_prop=`%Reads`) %>%
      filter(n_mutated > 0) 
   # out_tab=split(this,1:nrow(this))
    return(mytable)  
    
  }else {
    ##Filter on only those that have any number of mutations
    crispresso_file_name="Alleles_frequency_table.zip"
    
    myfiles = list.files(path=file, pattern = crispresso_file_name, recursive = T, full.names = T)
    myfiles.txt=sub(".zip",".txt",myfiles)
    names(myfiles) = gsub("^CRISPResso_on_","",basename(dirname(myfiles)))
    
    if(file.info(myfiles)$size>100000000){
      system(sprintf("unzip -o %s",myfiles))
      myfiles=myfiles.txt    
      }
    
    out_tab=lapply(myfiles,function(x){
      mytable=vroom::vroom(myfiles)
      this <- mytable %>% 
        dplyr::select(Aligned_Sequence,
                      Reference_Sequence,
                      n_deleted,n_inserted,n_mutated,Reads_n=`#Reads`,Reads_prop=`%Reads`) %>%

      filter(n_mutated > 0) 
    })
   names(out_tab)=names(myfiles)
    
  return(out_tab)  
    
  }
}


align_crispresso=function(Aligned_Sequence,Reference_Sequence,n_deleted,n_inserted,n_mutated,Reads_n,Reads_prop){
  
  ##Run Alignment
  alignment=pairwiseAlignment(DNAString(Aligned_Sequence),DNAString(gene_sequence),type="global-local")
  
  ##View Alignment
  
  #seq <- c(alignedPattern(alignment), alignedSubject(alignment))
  #DECIPHER::BrowseSeqs(seq)
  mmt=mismatchTable(alignment)
  
  
  
  ##INDELS FIRST TO ADD TO TABLE
  indT=indel(alignment)
  
  ins=data.frame(indT@insertion@unlistData)
  del=data.frame(indT@deletion@unlistData)
  
  
  
  vt_snp=tibble(chr=gene_coords$chr,Start=gene_coords$start+mmt$SubjectStart-1,
                End=gene_coords$start+mmt$SubjectEnd-1,
                REF=mmt$SubjectSubstring,ALT=mmt$PatternSubstring,
                Aligned=Aligned_Sequence,Reference=Reference_Sequence,
                n_deleted=n_deleted,n_inserted=n_inserted,n_mutated=n_mutated,
                Reads_n=Reads_n,Reads_prop=Reads_prop)
  
  
  if(nrow(ins)>0){
    
    ##INSERTION REF
    shift <- c(0L, head(cumsum(width(indT@insertion)[[1]]), n=-1L))
    ins_ranges <- shift(indT@insertion[[1]], shift)
   # end(ins_ranges)=ifelse(end(ins_ranges)>end(alignment@pattern@range),end(alignment@pattern@range),end(ins_ranges))
    #ins_seqs <- extractAt(alignment@pattern@unaligned[[1]], ins_ranges)
    
    
    aligned_pattern <- as(as.character(pattern(alignment)),
                          paste0(seqtype(pattern(alignment)), "String"))
    
    ins_seqs <- extractAt(aligned_pattern, ins_ranges)
    
    ins_seqs_dt=data.frame(ins_seqs@ranges)
    
    ins_alt=data.frame(ins_seqs)
    
    
    
    vt_ins=tibble(chr=gene_coords$chr,
                  Start=(gene_coords$start+alignment@subject@range@start+ins$start)-3,
                  End=gene_coords$start+ins$start+alignment@subject@range@start+ins$width-3,
                  REF="-",ALT=ins_alt$ins_seqs,
                  Aligned=Aligned_Sequence,Reference=Reference_Sequence,
                  n_deleted=n_deleted,n_inserted=n_inserted,n_mutated=n_mutated,
                  Reads_n=Reads_n,Reads_prop=Reads_prop)
    
    
  }
  
  if(nrow(del)>0){
    
    
    shift <- c(0L, head(cumsum(width(indT@deletion[[1]])), n=-1L))
    del_ranges <- shift(indT@deletion[[1]], shift)
    
    aligned_subject <- as(as.character(alignment@subject),
                          paste0(seqtype(subject(alignment)), "String"))
    
    del_seqs <- extractAt(aligned_subject, del_ranges)
    #end(del_ranges)=ifelse(end(del_ranges)>end(d@subject@range),end(d@pattern@range),end(del_ranges))
    
    
    
    #del_ranges=data.frame(start=del_ranges@start-sum(ins$width),end=del_ranges@start-sum(ins$width)+del_ranges@width)
    del_ranges=data.frame(del_seqs@ranges)
    #del_alt=data.frame(del_seqs)
    
    vt_del=tibble(chr=gene_coords$chr,
                  Start=(gene_coords$start+alignment@subject@range@start+del_ranges$start)-2,
                  End=(gene_coords$start+alignment@subject@range@start+del_ranges$end)-2,
                  REF=data.frame(del_seqs)$del_seqs,ALT="-",
                  Aligned=Aligned_Sequence,Reference=Reference_Sequence,
                  n_deleted=n_deleted,n_inserted=n_inserted,n_mutated=n_mutated,
                  Reads_n=Reads_n,Reads_prop=Reads_prop)
    
  }
  
  
  if(exists("vt_del") & exists("vt_ins")){
    bind_rows(vt_del,vt_ins,vt_snp) %>% arrange(Start) 
  }else if(exists("vt_ins")){
    bind_rows(vt_ins,vt_snp) %>% arrange(Start) 
  }else if(exists("vt_del")){
    bind_rows(vt_del,vt_snp) %>% arrange(Start) 
  }else{
    vt_snp %>% arrange(Start)
  }
  
}


  


aachange=function(POS,dt){
  REFnt=sapply(POS-gene_coords$start+1,function(x)as.character(subseq(DNAString(gene_sequence),x,x)))
  REFdt=data.frame(REF=REFnt,POS=POS) %>%
    left_join(.,dt,by=c("REF",c("POS"="Start")))
  newalt=ifelse(is.na(REFdt$ALT),REFdt$REF,REFdt$ALT)
  newref=paste(REFnt,collapse="")
  data.frame(REFnt=newref,ALTnt=paste(newalt,collapse=""),REF=as.character(translate(DNAString(newref))),ALT=as.character(translate(DNAString(paste(newalt,collapse="")))))
}

getGeneInfo=function(genename){
##Load the Gene Name into the name space right away!
if(file.exists(sprintf("%s_gene_coords.rds",genename)) & file.exists(sprintf("%s_gene_sequence.rds",genename))){

  assign('gene_coords',readRDS(sprintf("%s_gene_coords.rds",genename)),envir =globalenv())
  assign('gene_sequence',readRDS(sprintf("%s_gene_sequence.rds",genename)),envir = globalenv())
  
  }else{
  mart = useMart('ensembl', dataset="hsapiens_gene_ensembl")
  #mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl",host = "asia.ensembl.org")
  
  
  my_attrs = c(gene="external_gene_name",chr="chromosome_name",start="start_position",
               end="end_position",strand="strand",sequence="gene_exon_intron")
  
  gene_sequence_info = getBM(attributes = my_attrs,
                             filters = "external_gene_name", 
                             values = genename, mart = mart, verbose=T)
  
  
  gene_sequence_info = gene_sequence_info[,match(my_attrs, colnames(gene_sequence_info))]
  colnames(gene_sequence_info) = names(my_attrs)
  
  gene_coords = gene_sequence_info[gene_sequence_info$gene==genename, 
                                   c("chr","start","end","strand")]
  gene_sequence = gene_sequence_info$sequence[gene_sequence_info$gene==genename]
  
  saveRDS(gene_coords,sprintf("%s_gene_coords.rds",genename))
  saveRDS(gene_sequence,sprintf("%s_gene_sequence.rds",genename))
  
  assign('gene_coords',readRDS(sprintf("%s_gene_coords.rds",genename)),envir =globalenv())
  assign('gene_sequence',readRDS(sprintf("%s_gene_sequence.rds",genename)),envir = globalenv())
  
}
}

##Return AA annotation
return_AA_anno=function(split_dataframe,positions,annotin){
  lapply(split_dataframe,function(x){
    POS=positions
    AA1_change=aachange(POS,dt=x)
    reg_ex=paste0(gene_coords$chr[1],":",POS[1],"-",POS[3])
    tempa=annotin %>% filter(X2 %in% reg_ex) %>% dplyr::select(X3,REF_AA,ALT_AA,SIFT,Polyphen,CADD)
    
    data.frame(af.id=x$af.id[1],AA1_change) %>% left_join(.,tempa,by=c('REF'='REF_AA','ALT'="ALT_AA",'ALTnt'="X3"))%>%
      dplyr::rename(REF_AA=REF,ALT_AA=ALT,SIFT_AA=SIFT,Polyphen_AA=Polyphen,CADD_AA=CADD)
  })
}

