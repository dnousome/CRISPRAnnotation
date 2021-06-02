# CCBR1046 Project 
CRISPResso2 Annotation

Download all files 0_Functions.R AnnotatedCRISPResso.R, RunAnnotate.sh, runannovar.sh, and Annotation_render.Rmd  
  
 
### `RunAnnotate.sh`
This will Run the pipeline
```
usage: RunAnnotate.sh [-h] 

Run muh pipelinezz

optional arguments:
  --help                show this help message and exit
  --p PATH              [input_params] Path to directory containing all Alleles_Frequency_Table.zip
  --g Gene              [input_params] Gene of Intereset to 
```
Example
```
chmod +x RunAnnotate.sh
./RunAnnotate.sh -p /scratch/nousome/ccbr1046/Annotate_B1B2repeat/ -g BRCA2

```