#!/bin/bash

# ---------------
# Retrieve GDAS waves from AWS data archives, running download_AWSarchive_GDAS_Field.sh
# It generates a single output file for one month of continuous 6-hourly global data.
# Variables U10 and Hs.
# ---------------

# bash run_AWSarchive_GDAS_field.sh 2022 1 31 /work/noaa/marine/ricardo.campos/data/archiveOPruns/GDASwave

# Orion module load
module load nco

YEAR="$1"
MONTH="$2"
ndays="$3" # number of days for this month
DIRS="$4" # output path

cd ${DIRS}

YMTIME=${YEAR}`printf %2.2d $MONTH`
DIRW=${DIRS}/${YMTIME}
mkdir ${DIRW}
cd ${DIRW}

for DAY in $(seq $ndays); do
 WTIME=${YMTIME}`printf %2.2d $DAY`
 bash ${DIRS}/download_AWSarchive_GDAS_Field.sh ${WTIME} 00 ${DIRW}
 wait $!
 bash ${DIRS}/download_AWSarchive_GDAS_Field.sh ${WTIME} 12 ${DIRW}
 wait $!
done

# Merge all netcdf files for this month
ncrcat ${DIRW}/gdaswave.*.f*.nc -O ${DIRS}/gdaswave.${YMTIME}.global.0p25.nc
wait $!
sleep 1
rm -rf ${DIRW}

# permissions and groups
chmod -R 775 ${DIRS}/gdaswave.${YMTIME}.global.0p25.nc

