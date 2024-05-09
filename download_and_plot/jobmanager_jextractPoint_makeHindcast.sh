#!/bin/bash

DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/1preproc/extract_point_from_grib/jobs"
DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/1preproc/extract_point_from_grib"

cd ${DIRSCRIPTS}

for YEAR in `seq 2000 2019`; do
  for MONTH in `seq 1 12`; do

    export YEAR=${YEAR}
    export MONTH=${MONTH}
    sbatch --output=${DIRJOUT}"/jextractPoint_makeHindcast_"${YEAR}_${MONTH}".out" ${DIRSCRIPTS}"/jextractPoint_makeHindcast.sh"
    echo " job jextractPoint_makeHindcast_"${YEAR}_${MONTH}" submitted OK at "$(date +"%T")

  done
done

