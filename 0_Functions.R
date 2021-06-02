align_crispresso=function(x){
  d=pairwiseAlignment(DNAString(x$Aligned_Sequence),DNAString(gene_sequence),type="global-local")

  ##Check the alignment
  #seq <- c(alignedPattern(d), alignedSubject(d))
  #s1=DECIPHER::ConsensusSequence(seq)
  #DECIPHER::BrowseSeqs(seq)
  
  mmt=mismatchTable(d)
  
  
  
  ##INDELS FIRST TO ADD TO TABLE
  indT=indel(d)
  #print(c(paste0("INS=",length(indT@insertion@unlistData@start)),paste0("DEL=",length(indT@deletion@unlistData@start))))
  ins=data.frame(indT@insertion@unlistData)
  
  if(nrow(ins)>0){

      ins_ref=sapply(split(ins,1:nrow(ins)),function(y){

      as.character(subseq(
        DNAString(gene_sequence),(d@subject@range@start+y$start-2),
        d@subject@range@start+y$start-2))
      })
    
    ins_alt=lapply(split(ins,1:nrow(ins)),function(y){
      as.character(
        subseq(alignedPattern(d),(y$start),(y$end)))
      })
    
    ins_alt=mapply(function(z1,z2)paste0(z1,z2),ins_ref,ins_alt)
    
    
    vt_snp=tibble(chr=gene_coords$chr,Start=gene_coords$start+mmt$SubjectStart-1,
                  End=gene_coords$start+mmt$SubjectEnd-1,
                  REF=mmt$SubjectSubstring,ALT=mmt$PatternSubstring,
                  Aligned=x$Aligned_Sequence,Reference=x$Reference_Sequence,
                  n_deleted=x$n_deleted,n_inserted=x$n_inserted,n_mutated=x$n_mutated,
                  Reads_n=x$`#Reads`,Reads_prop=x$`%Reads`)
    
    vt_ins=tibble(chr=gene_coords$chr,
                  Start=(gene_coords$start+d@subject@range@start+ins$start)-3,
                  End=gene_coords$start+ins$start+d@subject@range@start+ins$width-3,
                  REF=ins_ref,ALT=ins_alt,
                  Aligned=x$Aligned_Sequence,Reference=x$Reference_Sequence,
                  n_deleted=x$n_deleted,n_inserted=x$n_inserted,n_mutated=x$n_mutated,
                  Reads_n=x$`#Reads`,Reads_prop=x$`%Reads`)
    
    #print(nrow(vt_snp)==x$n_mutated)
    bind_rows(vt_snp,vt_ins)
  }else{
    vt_snp=tibble(chr=gene_coords$chr,Start=gene_coords$start+mmt$SubjectStart-1,
                  End=gene_coords$start+mmt$SubjectEnd-1,
                  REF=mmt$SubjectSubstring,ALT=mmt$PatternSubstring,
                  Aligned=x$Aligned_Sequence,Reference=x$Reference_Sequence,
                  n_deleted=x$n_deleted,n_inserted=x$n_inserted,n_mutated=x$n_mutated,
                  Reads_n=x$`#Reads`,Reads_prop=x$`%Reads`)
    
    print(nrow(vt_snp)==x$n_mutated)
    vt_snp
  }
  
}



#PRocess ALLELEs frequency_table
allele_freq_tab=function(x){
  
  mytable = read_tsv(unz(x,"Alleles_frequency_table.txt"))


##Filter on only those that have any number of mutations
this <- mytable %>% 
  filter(n_mutated > 0) %>%
  filter(`#Reads`>5)



that=split(this,1:nrow(this))
lapply(that,align_crispresso)
}
