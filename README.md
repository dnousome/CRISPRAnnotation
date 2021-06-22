# CCBR1046 Project 
CRISPResso2 Annotation

Download all files 0_Functions.R 1_AnnotateCRISPResso.R, 2_runannovar.sh, and 3_Annotation_render.Rmd and the script to run it all CRISPR_Annotate.sh
  
 
### `CRISPR_Annotate.sh`
This will run the pipeline
```
usage: CRISPR_Annotate.sh [-h] 

arguments:
  --help                show this help message and exit
  --p PATH              [input_params] Path to directory containing all Alleles_Frequency_Table.zip files
  --g Gene              [input_params] Gene of Interest to grab coordinates (works only with 1 at a time)
```
Example
```
chmod +x RunAnnotate.sh
./RunAnnotate.sh -p /scratch/nousome/ccbr1046/Annotate_B1B2repeat/ -g BRCA2

```