#!/bin/bash

DIRJOUT="/work/noaa/marine/ricardo.campos/data/archiveOPruns/GDASwave/jobs"
DIRSCRIPTS="/work/noaa/marine/ricardo.campos/data/archiveOPruns/GDASwave"

cd ${DIRSCRIPTS}

YEAR="2024"
for MONTH in `seq 1 12`; do
    export YEAR=${YEAR}
    export MONTH=${MONTH}
    sbatch --output=${DIRJOUT}"/postproc_AWSarchive_GDAS_Field_"${YEAR}${MONTH}".out" ${DIRSCRIPTS}"/postproc_AWSarchive_GDAS_Field.sh"
    echo " job postproc_AWSarchive_GDAS_Field_"${YEAR}${MONTH}" submitted OK at "$(date +"%T")
    sleep 2
done

