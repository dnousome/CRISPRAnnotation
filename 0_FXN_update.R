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
                      n_deleted,n_inserted,n_mutated,`#Reads`,`%Reads`) %>%
        filter(n_mutated > 0) 
      split(this,1:nrow(this))
      
    })
    names(out_tab)=sheets
    return(out_tab)
    
  }else if(file_ext=="txt"){
    mytable = read_tsv(x)
    this <- mytable %>% 
      dplyr::select(Aligned_Sequence,
                    Reference_Sequence,
                    n_deleted,n_inserted,n_mutated,`#Reads`,`%Reads`) %>%
      filter(n_mutated > 0) 
   # out_tab=split(this,1:nrow(this))
    return(out_tab)  
    
  }else {
    ##Filter on only those that have any number of mutations
    crispresso_file_name="Alleles_frequency_table.zip"
    myfiles = list.files(path=file, pattern = crispresso_file_name, recursive = T, full.names = T)
    names(myfiles) = gsub("^CRISPResso_on_","",basename(dirname(myfiles)))
    
    out_tab=lapply(myfiles,function(x){
      mytable=read_tsv(unz(x,"Alleles_frequency_table.txt"))
      this <- mytable %>% 
        dplyr::select(Aligned_Sequence,
                      Reference_Sequence,
                      n_deleted,n_inserted,n_mutated,`#Reads`,`%Reads`) %>%

      filter(n_mutated > 0) 
      #split(this,1:nrow(this))
    })
   names(out_tab)=names(myfiles)
    
  return(out_tab)  
    
  }
}

align_crispresso_p1=function(x){
  d=pairwiseAlignment(DNAString(x$Aligned_Sequence),DNAString(gene_sequence),type="global-local")
}

align_crispresso_p2=function(d){
  ##Check the alignment
  #seq <- c(alignedSubject(d),alignedPattern(d))
  #DECIPHER::BrowseSeqs(seq)
  
  mmt=mismatchTable(d)
  
  
  
  ##INDELS FIRST TO ADD TO TABLE
  indT=indel(d)

  ins=data.frame(indT@insertion@unlistData)
  del=data.frame(indT@deletion@unlistData)
  
 
  
  vt_snp=tibble(chr=gene_coords$chr,Start=gene_coords$start+mmt$SubjectStart-1,
                End=gene_coords$start+mmt$SubjectEnd-1,
                REF=mmt$SubjectSubstring,ALT=mmt$PatternSubstring,
                Aligned=x$Aligned_Sequence,Reference=x$Reference_Sequence,
                n_deleted=x$n_deleted,n_inserted=x$n_inserted,n_mutated=x$n_mutated,
                Reads_n=x$`#Reads`,Reads_prop=x$`%Reads`)
  
  
  if(nrow(ins)>0){
    
    ##INSERTION REF
    shift <- c(0L, head(cumsum(width(indT@insertion)[[1]]), n=-1L))
    ins_ranges <- shift(indT@insertion[[1]], shift)
    end(ins_ranges)=ifelse(end(ins_ranges)>end(d@pattern@range),end(d@pattern@range),end(ins_ranges))
    ins_seqs <- extractAt(d@pattern@unaligned[[1]], ins_ranges)
    
    ins_seqs_dt=data.frame(ins_seqs@ranges)

    ins_alt=data.frame(ins_seqs)

    
    
    vt_ins=tibble(chr=gene_coords$chr,
                  Start=(gene_coords$start+d@subject@range@start+ins$start)-3,
                  End=gene_coords$start+ins$start+d@subject@range@start+ins$width-3,
                  REF="-",ALT=ins_alt$ins_seqs,
                  Aligned=x$Aligned_Sequence,Reference=x$Reference_Sequence,
                  n_deleted=x$n_deleted,n_inserted=x$n_inserted,n_mutated=x$n_mutated,
                  Reads_n=x$`#Reads`,Reads_prop=x$`%Reads`)
    
    
  }
  
  if(nrow(del)>0){
      
      
      shift <- c(0L, head(cumsum(width(indT@deletion[[1]])), n=-1L))
      del_ranges <- shift(indT@deletion[[1]], shift)
      del_seqs <- extractAt(d@subject@unaligned[[1]], del_ranges)
      #end(del_ranges)=ifelse(end(del_ranges)>end(d@subject@range),end(d@pattern@range),end(del_ranges))
      
      del_ranges=data.frame(start=del_ranges@start-sum(ins$width),end=del_ranges@start-sum(ins$width)+del_ranges@width)
    
      
      vt_del=tibble(chr=gene_coords$chr,
                    Start=(gene_coords$start+d@subject@range@start+del_ranges$start)-2,
                    End=(gene_coords$start+d@subject@range@start+del_ranges$end)-2,
                    REF=data.frame(del_seqs)$del_seqs,ALT="-",
                    Aligned=x$Aligned_Sequence,Reference=x$Reference_Sequence,
                    n_deleted=x$n_deleted,n_inserted=x$n_inserted,n_mutated=x$n_mutated,
                    Reads_n=x$`#Reads`,Reads_prop=x$`%Reads`)
      
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
