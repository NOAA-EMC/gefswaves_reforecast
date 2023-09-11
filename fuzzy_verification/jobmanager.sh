#!/bin/bash

DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/jobs"
DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification"

cd ${DIRSCRIPTS}

YEAR="2021"
for MONTH in `seq 7 8`; do
  for DAY in `seq 1 31`; do
    export YEAR=${YEAR}
    export MONTH=${MONTH}
    export DAY=${DAY}
    export DATE=${YEAR}$(printf "%02.f" $MONTH)$(printf "%02.f" $DAY)
    sbatch --output=${DIRJOUT}"/jbuildfuzzydataset_"${YEAR}${MONTH}${DAY}".out" ${DIRSCRIPTS}"/jbuildfuzzydataset.sh"
    echo " job jbuildfuzzydataset_"${YEAR}${MONTH}${DAY}" submitted OK at "$(date +"%T")
    sleep 2
  done
done

