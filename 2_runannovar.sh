#!/bin/bash

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -a|--annovarin)
    AVIN="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}"

echo "ANNOVAR AVINPUT = ${AVIN}"
if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi

OUTNAME="${AVIN/_toanno.avinput/}"    


module load annovar/2019-10-24
mkdir -p annovar
runcmd="table_annovar.pl $AVIN $ANNOVAR_DATA/hg38 -buildver hg38 -protocol gene,dbnsfp35a,dbnsfp41a,gnomad211_exome,gnomad30_genome,clinvar_20200419,cosmic92_coding,cosmic92_noncoding -operation g,f,f,f,f,f,f,f --remove --thread $SLURM_JOB_CPUS_PER_NODE --outfile annovar/$OUTNAME"
eval $runcmd