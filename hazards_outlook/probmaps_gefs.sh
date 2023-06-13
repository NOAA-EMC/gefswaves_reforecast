#!/bin/bash

########################################################################
# probmaps_gefs.sh and probmaps_gefs.py
#
# VERSION AND LAST UPDATE:
#   v1.0  06/09/2023
#
# PURPOSE:
#  Script to download NOAA Global Ensemble Forecast System (GEFS), marine
#   forecast, and generate the Probability Maps (Global Hazard Outlooks)
#
# USAGE:
#  This shell script is divided in two parts: the first one to download
#    operational forecast files, and the second one to produce the 
#    probability maps using a python script, probmaps_gefs.py.
#  The shell script requires one argument for the path where the 
#    codes and files are saved (or symbolic link): probmaps_gefs.sh, 
#    probmaps_gefs.py, probmaps_gefs.yaml; and where a directory will be
#    created to download the operational files.
#  Before starting the download, two auxiliar pearl scripts are downloaded,
#    get_grib.pl and get_inv.pl, which allows to fetch only specific
#    variables. If you already have those files, feel free to comment these
#    lines. They must be in the same directory you are running this script.
#  By default it downloads the forecast from the current day cycle (00Z),
#    with 6-h of time step up to 384 hours (16 days)
#  The python script probmaps_gefs.py is called at the end of this script, 
#    and it requires the configuration file probmaps_gefs.yaml. It is
#    IMPORTANT to edit the outpath variable in probmaps_gefs.yaml, which
#    is the location where the final plots will be saved. Once the .yaml
#    is configured, there is no need to change in the daily basis, unless
#    you want to modify the destination path or any other configuration.
#
#  Example:
#    ./probmaps_gefs.sh /home/user/path
#  or
#    bash probmaps_gefs.sh /home/user/path
#
# OUTPUT:
#  There are two outputs:
#  (1) One cycle of NOAA's GEFS Waves operational forecast, including all 
#    the 30 members plus the control member, in grib2 format. A directory
#    gefsWave.$YEAR$MONTH$DAY$HOUR is generated containing the data, plus
#    a log file logGEFS_$YEAR$MONTH$DAY$HOUR.
#  (2) Figures containing the probability maps, saved in the directory
#    outpath informed in the probmaps_gefs.yaml file
#
# DEPENDENCIES:
#  pearl and python (see dependencies in probmaps_gefs.py). 
#    In most linux systems, pearl is already included.
#  For the python installation, it is recommended:
#    https://www.anaconda.com/download
#
# AUTHOR and DATE:
#  06/09/2023: Ricardo M. Campos, first version 
#
# PERSON OF CONTACT:
#  Ricardo M. Campos: ricardo.campos@noaa.gov
#
########################################################################

# Path where the following codes and files are saved (or symbolic link):
# probmaps_gefs.sh, probmaps_gefs.py, probmaps_gefs.yaml
# and where a directory will be created to download operational files
DIR="$1"
# The output path where figures/plots will be saved is written in the
# configuration file probmaps_gefs.yaml

# NOAA server address
SERVER=https://ftpprd.ncep.noaa.gov/
s1="global.0p25" # main grid

# variable names to be downloaded. Wind speed, significant wave height, and peak period.
VARSGET=":WIND:surface:|:HTSGW:surface:|:PERPW:surface:"
# Corresponding variables for the python processing
MVARS="U10 Hs Tp"

# Initial date cycle for the ftp
YEAR=`date +%Y`
MONTH=`date +%m`
DAY=`date +%d`
# pa=2 #  days into the past. pa=1 dowloades data from yesterday's cycle
# YEAR=`date --date=-$pa' day' '+%Y'`
# MONTH=`date --date=-$pa' day' '+%m'`
# DAY=`date --date=-$pa' day' '+%d'`
HOUR="00" # first cycle 00Z

cd $DIR

# Auxiliar pearl scripts
# https://www.cpc.ncep.noaa.gov/products/wesley/fast_downloading_grib.html
wget --no-check-certificate --no-proxy -l1 -H -t1 -nd -N -np -erobots=off --tries=3 ftp://ftp.cpc.ncep.noaa.gov/wd51we/fast_downloading_grib/get_inv.pl
chmod 775 get_inv.pl
wget --no-check-certificate --no-proxy -l1 -H -t1 -nd -N -np -erobots=off --tries=3 ftp://ftp.cpc.ncep.noaa.gov/wd51we/fast_downloading_grib/get_grib.pl
chmod 775 get_grib.pl

# create directory
mkdir -p $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR
# all information about fetching and processing the grib2 files will be saved in the log file 
echo "  " > $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR

# hours of forecast lead time to be dowloaded. For now, using 6-h or time resolution
fleads="`seq -f "%03g" 0 6 384`"
# number of ensemble members
ensblm="`seq -f "%02g" 0 1 30`"
#
for h in $fleads;do
  for e in $ensblm;do

    echo "  " >> $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR
    echo " ======== GEFS Forecast: $YEAR$MONTH$DAY$HOUR  $h ========" >> $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 

    # size TAM and tries TRIES will control the process
    TAM=0
    TRIES=1
    # while file has lower size than expected it does:
    while [ $TAM -lt 1000000 ] && [ $TRIES -le 130 ]; do
      # sleep 5 minutes between attemps
      if [ ${TRIES} -gt 5 ]; then
        sleep 300
      fi

      if [ ${TAM} -lt 1000000 ]; then
          echo "  " >> $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR
          echo " attempt number: $TRIES" >> $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR
          echo "  " >> $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR
          # main line where get_inv.pl and get_grib.pl are used to fech the grib2 file
          if [ ${e} == 00 ]; then
             $DIR/get_inv.pl $SERVER/data/nccf/com/gens/prod/gefs.$YEAR$MONTH$DAY/$HOUR/wave/gridded/gefs.wave.t${HOUR}z.c${e}.${s1}.f"$(printf "%03.f" $h)".grib2.idx | egrep "($VARSGET)" | $DIR/get_grib.pl $SERVER/data/nccf/com/gens/prod/gefs.$YEAR$MONTH$DAY/$HOUR/wave/gridded/gefs.wave.t${HOUR}z.c${e}.${s1}.f"$(printf "%03.f" $h)".grib2 $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR/gefswave${e}.t${HOUR}z.pgrib2f"$(printf "%03.f" $h)".grib2 >> $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1 
          else
             $DIR/get_inv.pl $SERVER/data/nccf/com/gens/prod/gefs.$YEAR$MONTH$DAY/$HOUR/wave/gridded/gefs.wave.t${HOUR}z.p${e}.${s1}.f"$(printf "%03.f" $h)".grib2.idx | egrep "($VARSGET)" | $DIR/get_grib.pl $SERVER/data/nccf/com/gens/prod/gefs.$YEAR$MONTH$DAY/$HOUR/wave/gridded/gefs.wave.t${HOUR}z.p${e}.${s1}.f"$(printf "%03.f" $h)".grib2 $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR/gefswave${e}.t${HOUR}z.pgrib2f"$(printf "%03.f" $h)".grib2 >> $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1 
          fi
          # test if the downloaded file exists
          test -f $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR/gefswave${e}.t${HOUR}z.pgrib2f"$(printf "%03.f" $h)".grib2 >> $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
          TE=$?
          if [ ${TE} -eq 1 ]; then
            TAM=0
          else
            # check size of each file
            TAM=`du -sb $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR/gefswave${e}.t${HOUR}z.pgrib2f"$(printf "%03.f" $h)".grib2 | awk '{ print $1 }'` >> $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
          fi
      fi

      TRIES=`expr $TRIES + 1`
    done
  done
done
sleep 2
# permissions
chmod -R 775 $DIR/gefsWave.$YEAR$MONTH$DAY$HOUR
# Download complete. Starting python processing

# activate python 
source /home/ricardo/python/anaconda3/setanaconda3.sh

echo "  " > $DIR/logPythonOP_$YEAR$MONTH$DAY$HOUR
echo " PYTHON ENSEMBLE PROCESSING: GLOBAL HAZARDS OUTLOOK - PROBABILITY MAPS, $YEAR$MONTH$DAY$HOUR " >> $DIR/logPythonOP_$YEAR$MONTH$DAY$HOUR
# loop through variables
for WW3VAR in ${MVARS[*]}; do
  # the 7 14 is the time intervall (days), so week 2 is from day 7 to day 14 (included)
  python3 probmaps_gefs.py $YEAR$MONTH$DAY$HOUR 7 14 ${WW3VAR} >> $DIR/logPythonOP_$YEAR$MONTH$DAY$HOUR 2>&1
  wait $!
  echo "  " >> $DIR/logPythonOP_$YEAR$MONTH$DAY$HOUR
  echo "     ${WW3VAR} Done for Week 2." >> $DIR/logPythonOP_$YEAR$MONTH$DAY$HOUR
  echo "  " >> $DIR/logPythonOP_$YEAR$MONTH$DAY$HOUR
done

echo "  PYTHON PROCESSING COMPLETE" >> $DIR/logPythonOP_$YEAR$MONTH$DAY$HOUR

