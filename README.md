# CCBR1046 Project 
CRISPResso2 Annotation

Download all files 0_Functions.R 1_AnnotateCRISPResso.R, 2_runannovar.sh, and 3_Annotation_render.Rmd and the script to run it all CRISPR_Annotate.sh
  
 
### `CRISPR_Annotate.sh`
Script to run the full pipeline. 
Change the .sh to executable before running
```
usage: CRISPR_Annotate.sh [-h] 

arguments:
  --help                show this help message and exit
  --path/-p PATH              [input_params] Path to directory containing all Alleles_Frequency_Table.zip files
  --gene/-g Gene              [input_params] Gene of Interest to grab coordinates (works only with 1 at a time)
```
Example
```
chmod +x CRISPR_Annotate.sh
./CRISPR_Annotate.sh -p /scratch/nousome/ccbr1046/Annotate_B1B2repeat/ -g BRCA2

```