#!/bin/bash

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"
         
case $key in
  -p|--path)
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
  -pam|--pamsite)
  PAMSITE="$2"
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
echo "END = ${END}"
echo "PAMSite = ${PAMSITE}"

if [[ -n $1 ]]; then
echo "Last line of file specified as non-opt/last argument:"
tail -1 $1
fi

SUBMIT_SCRIPT="CRISPR_annot_submit.slurm"
subdate=$(date +%F_%H%M%S)
echo -e "#!/usr/bin/bash\nmodule load R/4.0\nRscript 1_AnnotateCRISPResso.R -p $CRISPPath -g $GENE -s $START -e $END -pam $PAMSITE"> $SUBMIT_SCRIPT
echo "Submitting pipeline to cluster... "

primaryID=$(sbatch --cpus-per-task=16 --mem=64g --time 48:00:00 --partition norm --output submit_"$subdate".log --error error_"$subdate".log $SUBMIT_SCRIPT)
#primaryID=$(sbatch --cpus-per-task=2 --mem=20g --time 5-00:00:00 --partition ccr,norm --output submit.log --error submit.log $SUBMIT_SCRIPT)
echo "Primary Job ID: $primaryID"
