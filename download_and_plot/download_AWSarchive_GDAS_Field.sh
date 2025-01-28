#!/bin/bash

# ---------------
# Download GDAS wave files from AWS
# convert to netcdf format, select U10 and Hs, and compress.
# This code is run inside run_AWSarchive_GDAS_field.sh
# ---------------

# 3 input arguments
# cycle date and time
CTIME="$1" # YYYYMMDD
CHOUR="$2" # HH
# destination path
DIRW="$3"
# server address
SERVER=https://noaa-gfs-bdp-pds.s3.amazonaws.com/
# model
model="global.0p16"
# cutoff decimals to reduce file size
dp=3
# Forecast lead time (hours) to download
fleads="`seq -f "%03g" 0 1 5`"

cd ${DIRW}
for h in $fleads;do
    echo " ======== GDAS, AWS archive: ${CTIME} ${CHOUR}Z  $h ========"

    # size TAM and tries TRIES will control the process
    TAM=0
    TRIES=1

    while [ $TAM -lt 8000000 ] && [ $TRIES -le 130 ]; do
      # sleep 5 minutes between attemps
      if [ ${TRIES} -gt 5 ]; then
        sleep 30
      fi

      if [ ${TAM} -lt 8000000 ]; then
          # Main line, download
          wget -l1 -H -t1 -nd -N -np -erobots=off --tries=3 ${SERVER}gdas.${CTIME}/"$(printf "%02.f" $CHOUR)"/wave/gridded/gdaswave.t"$(printf "%02.f" $CHOUR)"z.${model}.f"$(printf "%03.f" $h)".grib2 -O $DIRW/gdaswave.${CTIME}.t"$(printf "%02.f" $CHOUR)"z.${model}.f"$(printf "%03.f" $h)".grib2 2>&1
          wait $!

          # test if the downloaded file exists
          test -f $DIRW/gdaswave.${CTIME}.t"$(printf "%02.f" $CHOUR)"z.${model}.f"$(printf "%03.f" $h)".grib2
          TE=$?
          if [ ${TE} -eq 1 ]; then
            TAM=0
          else
            # check size of each file
            TAM=`du -sb $DIRW/gdaswave.${CTIME}.t"$(printf "%02.f" $CHOUR)"z.${model}.f"$(printf "%03.f" $h)".grib2 | awk '{ print $1 }'` 2>&1
          fi
          echo $DIRW/gdaswave.${CTIME}.t"$(printf "%02.f" $CHOUR)"z.${model}.f"$(printf "%03.f" $h)".grib2
      fi

      TRIES=`expr $TRIES + 1`
    done
done

echo " Done ${CTIME}.${CHOUR}"

