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
  *)    # unknown option
  POSITIONAL+=("$1") # save it in an array for later
  shift 
  ;;
esac
done
set -- "${POSITIONAL[@]}"

echo "CRISPResso Path = ${CRISPPath}"
echo "Gene = ${GENE}"

if [[ -n $1 ]]; then
echo "Last line of file specified as non-opt/last argument:"
tail -1 $1
fi

SUBMIT_SCRIPT="CRISPanno.sh"

echo -e "#!/usr/bin/bash\nmodule load R/4.0\nRscript AnnotateCRISPResso.R -p $CRISPPath -g $GENE"> $SUBMIT_SCRIPT
echo "Submitting pipeline to cluster... "
primaryID=$(sbatch --cpus-per-task=2 --mem=64g --time 25:00:00 --partition norm --output submit.log --error submit.log $SUBMIT_SCRIPT)
#primaryID=$(sbatch --cpus-per-task=2 --mem=20g --time 5-00:00:00 --partition ccr,norm --output submit.log --error submit.log $SUBMIT_SCRIPT)
echo "Primary Job ID: $primaryID"
