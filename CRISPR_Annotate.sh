#!/bin/bash

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"
         
case $key in
  -i|--inputpath)
  CRISPPath="$2"
  shift # past argument
  shift # past value
  ;;
  -g|--gene)
  GENE="$2"
  shift # past argument
  shift # past value
  ;;
  -s|--start)
  START="$2"
  shift # past argument
  shift # past value
  ;;
  -e|--end)
  END="$2"
  shift # past argument
  shift # past value
  ;;
  -p|--pamsite)
  PAMSITE="$2"
  shift # past argument
  shift # past value
  ;;
  -a|--pamsiteallele)
  PAMSITEALLELE="$2"
  shift # past argument
  shift # past value
  ;;
  -o|--output)
  OUTPUT="$2"
  shift # past argument
  shift # past value
  ;;
  *)    
  POSITIONAL+=("$1") # save it in an array for later
  shift 
  ;;
esac
done
set -- "${POSITIONAL[@]}"

echo "CRISPResso Path = ${CRISPPath}"
echo "Gene = ${GENE}"
echo "Start = ${START}"
echo "End = ${END}"
echo "PAMSite = ${PAMSITE}"
echo "PAMSiteAlleles = ${PAMSITEALLELE}"
echo "Output Name = ${OUTPUT}"

if [[ -n $1 ]]; then
echo "Last line of file specified as non-opt/last argument:"
tail -1 $1
fi

SUBMIT_SCRIPT="CRISPR_annot_submit_$OUTPUT.slurm"
subdate=$(date +%F_%H%M%S)
echo -e "#!/usr/bin/bash\nmodule load R/4.1\nRscript 1_AnnotateCRISPResso.R -i $CRISPPath -g $GENE -s $START -e $END -p $PAMSITE -a $PAMSITEALLELE -o $OUTPUT"> $SUBMIT_SCRIPT
echo "Submitting pipeline to cluster... "

primaryID=$(sbatch --cpus-per-task=8 --mem=20g --time 12:00:00 --partition norm --output submit_"$subdate".log --error error_"$subdate".log $SUBMIT_SCRIPT)
#primaryID=$(sbatch --cpus-per-task=2 --mem=20g --time 5-00:00:00 --partition ccr,norm --output submit.log --error submit.log $SUBMIT_SCRIPT)
echo "Primary Job ID: $primaryID"
