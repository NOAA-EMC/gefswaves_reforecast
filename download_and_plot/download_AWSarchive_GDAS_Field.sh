#!/bin/bash

# ---------------
# Download GDAS wave files from AWS
# convert to netcdf format, select U10 and Hs, and compress.
# This code is run inside run_AWSarchive_GDAS_field.sh
# ---------------

# Orion module loads
module load wgrib2
module load nco
module load cdo

# 3 input arguments
# cycle date and time
CTIME="$1" # YYYYMMDD
CHOUR="$2" # HH
# destination path
DIRW="$3"
# server address
SERVER=https://noaa-gfs-bdp-pds.s3.amazonaws.com/
# cutoff decimals to reduce file size
dp=2
# Forecast lead time (hours) to download
fleads="`seq -f "%03g" 0 6 6`"

cd ${DIRW}
for h in $fleads;do
    echo " ======== GDAS, AWS archive: ${CTIME} ${CHOUR}Z  $h ========"

    # size TAM and tries TRIES will control the process
    TAM=0
    TRIES=1

    while [ $TAM -lt 10000000 ] && [ $TRIES -le 130 ]; do
      # sleep 5 minutes between attemps
      if [ ${TRIES} -gt 5 ]; then
        sleep 30
      fi

      if [ ${TAM} -lt 10000000 ]; then
          # Main line, download
          wget -l1 -H -t1 -nd -N -np -erobots=off --tries=3 ${SERVER}gdas.${CTIME}/"$(printf "%02.f" $CHOUR)"/wave/gridded/gdaswave.t"$(printf "%02.f" $CHOUR)"z.global.0p25.f"$(printf "%03.f" $h)".grib2 -O $DIRW/gdaswave.${CTIME}.t"$(printf "%02.f" $CHOUR)"z.global.0p25.f"$(printf "%03.f" $h)".grib2 2>&1
          wait $!

          # test if the downloaded file exists
          test -f $DIRW/gdaswave.${CTIME}.t"$(printf "%02.f" $CHOUR)"z.global.0p25.f"$(printf "%03.f" $h)".grib2
          TE=$?
          if [ ${TE} -eq 1 ]; then
            TAM=0
          else
            # check size of each file
            TAM=`du -sb $DIRW/gdaswave.${CTIME}.t"$(printf "%02.f" $CHOUR)"z.global.0p25.f"$(printf "%03.f" $h)".grib2 | awk '{ print $1 }'` 2>&1
          fi
          echo $DIRW/gdaswave.${CTIME}.t"$(printf "%02.f" $CHOUR)"z.global.0p25.f"$(printf "%03.f" $h)".grib2
      fi

      TRIES=`expr $TRIES + 1`
    done
done

# Post-processing: convert to netcdf, compress, and reduce decimals resolution, to save disk space. ------------------
for h in $fleads;do

  arqn=$DIRW/gdaswave.${CTIME}.t"$(printf "%02.f" $CHOUR)"z.global.0p25.f"$(printf "%03.f" $h)"
  test -f ${arqn}.grib2
  TE=$?
  if [ "$TE" -eq 1 ]; then
     echo " File ${arqn}.grib2 does not exist. Failed to download " 
  else
     wgrib2 ${arqn}.grib2 -netcdf ${arqn}.saux.nc  2>&1
     wait $!
     cdo select,name='WIND_surface','HTSGW_surface' ${arqn}.saux.nc ${arqn}.saux2.nc
     wait $!
     ncks --ppc default=.$dp ${arqn}.saux2.nc ${arqn}.nc  2>&1
     wait $!
     ncatted -a _FillValue,,o,f,NaN ${arqn}.nc  2>&1
     wait $!
     rm -f ${arqn}.grib2
     rm ${arqn}.saux*
     echo " File ${arqn} converted to netcdf and compressed with success. " 
  fi

done

echo " Done ${CTIME}.${CHOUR}"

