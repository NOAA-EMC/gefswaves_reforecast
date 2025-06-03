#!/bin/bash

########################################################################
# download_AWSarchive_GEFS_Field_st.sh
# download short term time steps not downloaded before
#
# VERSION AND LAST UPDATE:
#   v1.0  02/15/2023
#   v1.1  04/30/2025
#
# PURPOSE:
#  Script to download NOAA Global Ensemble Forecast System (GEFS), Wave 
#   Forecast from WAVEWATCH III operational. Download from AWS archive
#   and save the grib2 files without any conversion or processing. It
#   includes the control and all perturbed members of the ensemble.
#   Global wind and wave fields.
#
# USAGE:
#  Two input arguments, date and path, must be entered.
#  Example:
#    bash download_AWSarchive_GEFS_Field.sh 20220823 /home/ricardo/data/gefs
#
# OUTPUT:
#  Multiple grib2 files, for each time step and ensemble member.
#
# DEPENDENCIES:
#  wget
#
# AUTHOR and DATE:
#  02/15/2023: Ricardo M. Campos, first version 
#  04/30/2025: Ricardo M. Campos, flexible cycle time
#
# PERSON OF CONTACT:
#  Ricardo M. Campos: ricardo.campos@noaa.gov
#
#  If you are interested in operational forecasts from NOAA ftp, see:
#  https://github.com/NOAA-EMC/WW3-tools/tree/develop/opforecast
#
########################################################################

# Two input arguments
# date
CTIME="$1"
# cycle 00,06,12,18
HCYCLE="$2"
# destination path
DIRW="$3"

# server address
SERVER=https://noaa-gefs-pds.s3.amazonaws.com/
# ensemble members
ensblm="`seq -f "%02g" 0 1 30`"
# Forecast lead time (hours) to download
fleads="`seq -f "%03g" 3 6 168`"

cd ${DIRW}
for h in $fleads; do
  echo " ======== GEFS Forecast, AWS archive: ${CTIME} ${HCYCLE}Z  $h ========"
  for e in $ensblm; do
    echo "$e"
    
    FILE=$DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f$(printf "%03.f" $h).grib2

    # Skip if file exists and is large enough
    if [ -f "$FILE" ]; then
      TAM=$(du -sb "$FILE" | awk '{ print $1 }')
      if [ "$TAM" -ge 10000000 ]; then
        echo "File $FILE already exists and is large enough. Skipping download."
        continue
      fi
    fi

    # Initialize TAM and TRIES for download loop
    TAM=0
    TRIES=1

    while [ $TAM -lt 10000000 ] && [ $TRIES -le 130 ]; do
      if [ $TRIES -gt 5 ]; then
        sleep 30
      fi

      # Download logic
      if [ "$e" == "00" ]; then
        wget -l1 -H -t1 -nd -N -np -erobots=off --tries=3 \
          "${SERVER}gefs.${CTIME}/${HCYCLE}/wave/gridded/gefs.wave.t${HCYCLE}z.c${e}.global.0p25.f$(printf "%03.f" $h).grib2" \
          -O "$FILE" 2>&1
      else
        wget -l1 -H -t1 -nd -N -np -erobots=off --tries=3 \
          "${SERVER}gefs.${CTIME}/${HCYCLE}/wave/gridded/gefs.wave.t${HCYCLE}z.p${e}.global.0p25.f$(printf "%03.f" $h).grib2" \
          -O "$FILE" 2>&1
      fi
      wait $!

      # Check file existence and size
      if [ -f "$FILE" ]; then
        TAM=$(du -sb "$FILE" | awk '{ print $1 }')
      else
        TAM=0
      fi

      echo "$FILE (Size: $TAM)"
      TRIES=$((TRIES + 1))
    done
  done
done

# remove empty data
find $DIR -empty -type f -delete

echo " Done ${CTIME}."

