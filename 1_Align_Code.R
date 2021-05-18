#################
#
#CRISPRESSO Alignment Project
#By Darryl Nousome
#
#################
library(tidyverse)
library(biomaRt)
library(Biostrings)
library(GenomicRanges)


rm(list=ls())
setwd=("~/Documents/projects/ccbr1046/combined_flowcells/")


batch_dir="~/Documents/projects/ccbr1046/combined_flowcells/CRISPRessoBatch_on_batch"

crispresso_file_name="Alleles_frequency_table.zip"
myfiles = list.files(path=batch_dir, pattern = crispresso_file_name, recursive = T, full.names = T)

names(myfiles) <- gsub("^CRISPResso_on_","",basename(dirname(myfiles)))


##Load the Gene Name
mygenes = c("BRCA2")
mart = useMart('ensembl', dataset="hsapiens_gene_ensembl")
#mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl",host = "asia.ensembl.org")


#all_attr = listAttributes(mart)
my_attrs = c(gene="external_gene_name",chr="chromosome_name",start="start_position",
              end="end_position",strand="strand",sequence="gene_exon_intron")

gene_sequence_info = getBM(attributes = my_attrs,
                           filters = "external_gene_name", 
                           values = mygenes, mart = mart, verbose=T)


gene_sequence_info = gene_sequence_info[,match(my_attrs, colnames(gene_sequence_info))]
colnames(gene_sequence_info) = names(my_attrs)


#colnames(gene_sequence_info)
#nchar(gene_sequence_info$sequence)



gene_coords = gene_sequence_info[gene_sequence_info$gene==mygenes[1], 
                                 c("chr","start","end","strand")]
gene_sequence = gene_sequence_info$sequence[gene_sequence_info$gene==mygenes[1]]



##Read in the data
##Source the functions
source("0_Functions.R")
out=lapply(myfiles,allele_freq_tab)

##For annovar Keep on the REF/ALT
vt=lapply(out,function(x){
  bind_rows(x) %>% arrange(Start) %>% distinct() 
})

vt_annovar=lapply(out,function(x){
  bind_rows(x) %>% arrange(Start) %>% distinct() %>%
    dplyr::select(chr,Start,End,REF,ALT)
})

##Load this VT_SNP object into ANNOVAR And see?
lapply(names(vt_annovar),function(x){
  write_tsv(vt[[x]],paste0(x,"_toanno.avinput"),col_names = F)
})



###Run ANNOVAR output
##Use 2_runanno.sh











############Test code
d1=pairwiseAlignment(DNAString(this$Aligned_Sequence[27]),DNAString(gene_sequence),type="overlap")
matchPattern()
seq <- c(alignedPattern(d1), alignedSubject(d1))
DECIPHER::BrowseSeqs(seq)
s1=DECIPHER::ConsensusSequence(seq)



Views(d)
che@ranges
str_length(gene_sequence)
str_sub(gene_sequence,che@ranges@start,che@ranges@start+che@ranges@width)
this$Aligned_Sequence[3]
##Create some sort of GRanges object?



#gene_aln <- pairwiseAlignment(gsub("-","",this$Reference_Sequence[1]), subject = gene_sequence, type="local")
gene_aln <- pairwiseAlignment(this$Reference_Sequence[1], subject = gene_sequence, type="local")

# gene_aln <- pairwiseAlignment(gsub("-","",this$Reference_Sequence), subject = ref_sequence, type="local")
# gene_aln <- my_aln
myviews <- Views(gene_aln)

library(GenomicRanges)
range1 <- with(gene_coords,GRanges(chromosome_name, IRanges(start=start_position,end=end_position), strand="+"))
getSeq(DNAString(gene_sequence,start = gene_coords$start_position), range1)

d=DNAString(gene_sequence)
d[48009:48310]

ob=DNAString(gsub("-","",this$Reference_Sequence[1]))
ob
che=matchPattern(DNAString(gsub("-","",this$Aligned_Sequence[3])),DNAString(gene_sequence),with.indels = T,algorithm = "boyer-moore")

che=pairwiseAlignment(DNAString(this$Aligned_Sequence[3]),DNAString(this$Reference_Sequence[3]))


DNA_ALPHABET
che=matchPattern(DNAString(this$Reference_Sequence[1]),DNAString(gene_sequence))


chr_coordinate_data <- data.frame( start_pos_gene=start(myviews),
                                   end_pos_gene=end(myviews),
                                   start_pos_chr=start(myviews) + gene_coords$start,
                                   end_pos_gene=end(myviews) + gene_coords$start,
                                   chr = gene_coords$chr,
                                   strand = gene_coords$strand,
                                   stringsAsFactors = F)







sum(nmismatch(gene_aln) > 0)
sum(nindel(gene_aln)@insertion[,"Length"] > 0)
compareStrings(gene_aln)
class(gene_aln)
aligned(gene_aln)
alignedPattern(gene_aln)
alignedSubject(gene_aln)

pattern(my_aln)
aligned(pattern(my_aln))


sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)

my_aln <- pairwiseAlignment(pattern = gsub("-","",this$Aligned_Sequence[1]), subject = gsub("-","",this$Reference_Sequence[1]))
my_aln <- pairwiseAlignment(pattern = gsub("-","",this$Aligned_Sequence[1]), 
                            subject = gsub("-","",this$Reference_Sequence[1]),substitutionMatrix=sigma)


my_aln <- pairwiseAlignment(pattern = gsub("-","",this$Aligned_Sequence[1:10]), subject = gsub("-","",this$Reference_Sequence[1:10]), type="local")

myviews <- Views(gene_aln)
mismatchTable(my_aln)

lapply(my_aln, mismatchTable)

nmismatch(my_aln)
nindel(my_aln)
as.matrix(nindel(my_aln)@insertion)
mismatchSummary(my_aln)
compareStrings(my_aln)

