#!/bin/bash

DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/jobs"
DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification"

cd ${DIRSCRIPTS}

YEAR="2020"
for MONTH in `seq 1 12`; do
    export YEAR=${YEAR}
    export MONTH=${MONTH}
    sbatch --output=${DIRJOUT}"/jextract_GEFS_"${YEAR}${MONTH}".out" ${DIRSCRIPTS}"/jextract_GEFS.sh"
    echo " job jextract_GEFS_"${YEAR}${MONTH}" submitted OK at "$(date +"%T")
    sleep 1
  done
done


