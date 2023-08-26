#!/bin/bash

########################################################################
# download_and_proc_AWSarchive_Field.sh
#
# VERSION AND LAST UPDATE:
#   v1.0  02/15/2023
#
# PURPOSE:
#  Script to download NOAA Global Ensemble Forecast System (GEFS), Wave 
#   Forecast from WAVEWATCH III operational. Download grib2 files from
#   AWS archive, convert to netcdf format and compress. Only the control
#   member and the ensemble mean are obtained.
#   Global wind and wave fields.
#
# USAGE:
#  Two input arguments, date and output path, must be entered.
#  Example:
#    bash download_and_proc_AWSarchive_Field.sh 20220823 /home/ricardo/data/gefs
#
# OUTPUT:
#  Two netcdf files, one with the control member (c00) and one with the
#    ensemble mean (mean). Two variables included in each file, wind speed
#    and significant wave height.
#
# DEPENDENCIES:
#  wget, pearl, NCO, wgrib2
#  For most linux systems, wget and pearl are already included. The 
#    program NCO can be downloaded via apt-get (debian/ubuntu):
#    sudo apt-get install nco
#    the program wgrib2 can be downloaded at
#    https://www.ftp.cpc.ncep.noaa.gov/wd51we/wgrib2/
#
#  Before starting the download, two auxiliar pearl scripts are downloaded,
#    get_grib.pl and get_inv.pl, which allows to fetch only specific
#    variables. If you already have those files, feel free to comment these
#    lines. They must be in the same directory you are running this script.
#
# AUTHOR and DATE:
#  02/15/2023: Ricardo M. Campos, first version 
#
# PERSON OF CONTACT:
#  Ricardo M. Campos: ricardo.campos@noaa.gov
#
#  If you are interested in operational forecasts from NOAA ftp, see:
#  https://github.com/NOAA-EMC/WW3-tools/tree/develop/opforecast
#
########################################################################

# date
CTIME="$1"
# path
DIRW="$2"
# server address
SERVER=https://noaa-gefs-pds.s3.amazonaws.com/

# Auxiliar pearl scripts
# https://www.cpc.ncep.noaa.gov/products/wesley/fast_downloading_grib.html
wget --no-check-certificate --no-proxy -l1 -H -t1 -nd -N -np -erobots=off --tries=3 ftp://ftp.cpc.ncep.noaa.gov/wd51we/fast_downloading_grib/get_inv.pl
chmod 775 get_inv.pl
wget --no-check-certificate --no-proxy -l1 -H -t1 -nd -N -np -erobots=off --tries=3 ftp://ftp.cpc.ncep.noaa.gov/wd51we/fast_downloading_grib/get_grib.pl
chmod 775 get_grib.pl

# convertion and compression
WGRIB2=$(which wgrib2)
NCKS=$(which ncks)
NCATTED=$(which ncatted)
NCRCAT=$(which ncrcat)

# cutoff decimals to reduce file size
dp=2
# variable names to be downloaded.
VARSGET=":WIND:surface:|:HTSGW:surface:"
NVARSGET="WIND_surface,HTSGW_surface"

cd $DIRW

# hours of forecast lead time to be dowloaded
h1="`seq -f "%03g" 0 3 240`"
h2="`seq -f "%03g" 246 6 384`"
fleads=${h1}' '${h2}
# ensemble members
ensblm="c00 mean"

for h in $fleads;do
  echo " ======== GEFS Forecast, AWS archive: ${CTIME} 00Z  $h ========"
  for e in $ensblm;do
    echo $e
    # size TAM and tries TRIES will control the process
    TAM=0
    TRIES=1

    while [ $TAM -lt 1000000 ] && [ $TRIES -le 130 ]; do
      # sleep 5 minutes between attemps
      if [ ${TRIES} -gt 5 ]; then
        sleep 30
      fi

      if [ ${TAM} -lt 1000000 ]; then
          # Main line, download
          if [ ${e} == "c00" ]; then
            $DIRW/get_inv.pl ${SERVER}gefs.${CTIME}/00/wave/gridded/gefs.wave.t00z.${e}.global.0p25.f"$(printf "%03.f" $h)".grib2.idx | egrep "($VARSGET)" | $DIRW/get_grib.pl ${SERVER}gefs.${CTIME}/00/wave/gridded/gefs.wave.t00z.${e}.global.0p25.f"$(printf "%03.f" $h)".grib2 $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f"$(printf "%03.f" $h)".grib2 2>&1              
            wait $!
          else
            wget -l1 -H -t1 -nd -N -np -erobots=off --tries=3 ${SERVER}gefs.${CTIME}/00/wave/gridded/gefs.wave.t00z.${e}.global.0p25.f"$(printf "%03.f" $h)".grib2 -O $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f"$(printf "%03.f" $h)".grib2 2>&1
            wait $!
          fi
          # test if the downloaded file exists
          test -f $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f"$(printf "%03.f" $h)".grib2
          TE=$?
          if [ ${TE} -eq 1 ]; then
            TAM=0
          else
            # check size of each file
            TAM=`du -sb $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f"$(printf "%03.f" $h)".grib2 | awk '{ print $1 }'` 2>&1
          fi
          echo $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f"$(printf "%03.f" $h)".grib2
      fi

      TRIES=`expr $TRIES + 1`
    done
  done
done
sleep 1

cd $DIRW

# Post-processing: compress, select variables, and reduce decimals resolution, to save disk space. ------------------
echo " "
echo " Post-Processing. Netcdf4 compression "${CTIME}

for h in $fleads;do
  for e in $ensblm;do

     arqn=$DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f"$(printf "%03.f" $h)"
     test -f ${arqn}.grib2
     TE=$?
     if [ ${TE} -eq 1 ]; then
       echo " File ${arqn}.grib2 does not exist. Download failed."
     else
       $WGRIB2 ${arqn}.grib2 -netcdf ${arqn}.saux.nc 2>&1
       wait $!
       $NCKS -A -v $NVARSGET ${arqn}.saux.nc ${arqn}.saux2.nc 2>&1
       wait $!
       $NCKS -4 -L 1 ${arqn}.saux2.nc ${arqn}.saux3.nc 2>&1
       wait $!
       $NCKS --ppc default=.$dp ${arqn}.saux3.nc ${arqn}.nc 2>&1
       wait $!
       $NCATTED -a _FillValue,,o,f,NaN ${arqn}.nc 2>&1
       wait $!
       rm -f ${arqn}.grib2
       rm -f ${arqn}.saux*
       echo " File ${arqn} converted to netcdf and compressed with success. "
       sleep 1
     fi

  done
done

# Merge all netcdf files
for e in $ensblm;do
  $NCRCAT $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f*.nc -O $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.nc 2>&1
  wait $!
  sleep 1
  rm -f $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f*.nc
  echo " All time steps of ensemble member ${e} have been merged. "
done

