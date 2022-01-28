#!/usr/bin/env python3
import argparse,os,subprocess

def parse_args():
  parser = argparse.ArgumentParser(description='Input files')
  parser.add_argument('--fq1',help='fastq1')
  parser.add_argument('--fq2',help='fastq2')
  args = parser.parse_args()
  return(args)


def main():
  args=parse_args()
  cmd=["singularity exec -B /data/CCBR/rawdata/nousome/1046/Exon21b,/data/CCBR_Pipeliner/db/PipeDB/dev/","/data/CCBR/rawdata/nousome/1046/crispresso2_v2.1.1.sif","CRISPResso","--fastq_r1", args.fq1, "--fastq_r2",args.fq2,
       "--amplicon_seq CTTTGGGTGTTTTATGCTTGGTTCTTTAGTTTTAGTTGCTTTTGAATTTACAGTTTAGTGAATTAATAATCCTTTTGTTTTCTTAGAAAACACAACAAAACCATATTTACCATCACGTGCACTAACAAGACAGCAAGTTCGTGCTTTGCAAGATGGTGCAGAGCTTTATGAAGCAGTGAAGAATGCAGCAGACCCAGCTTACCTTGAGGTGAGAGAGTAAGAGGACATATAATGAGGCTTGAT",
"--trim_sequences --trimmomatic_options_string ILLUMINACLIP:/data/CCBR_Pipeliner/db/PipeDB/dev/adapters2.fa:0:90:10:0:true",
"--plot_window_size 50 --quantification_window_size 0 --max_rows_alleles_around_cut_to_plot 100 --min_frequency_alleles_around_cut_to_plot 0.001 --place_report_in_output_folder --write_detailed_allele_table"]
  cmd1=' '.join(cmd)
  os.system(cmd1)
    
if __name__=="__main__":
  main()
