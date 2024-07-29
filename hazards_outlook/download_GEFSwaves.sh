#!/bin/bash

########################################################################
# download_GEFSwaves.sh
#
# VERSION AND LAST UPDATE:
#   v1.0  11/09/2022
#   v1.1  04/01/2024
#
# PURPOSE:
#  Script to download NOAA Global Ensemble Forecast System (GEFS), Wave 
#   Forecast from WAVEWATCH III operational. Download from ftp, not nomads.
#
# USAGE:
#  There is only one input argument needed, containing the path where 
#    this script is located (or symbolic linked) and where data will 
#    be downloaded. See DIR below.
#  The variable VARSGET has the GEFS variables to download. If you 
#    change it, re-evaluate the variable FSIZE, associated with the 
#    expected size of each grib2 file downloaded.
#  This script does not convert/modify format and does not compress files.
#  By default it downloads the operational forecast of the current day.
#    If you want files from previous days, uncomment the lines and 
#    set pa=1
#
# OUTPUT:
#  One cycle of NOAA's GEFS Waves operational forecast, including 31 
#    members (control + 30 perturbed members), in grib2 format.
#  A directory gefs.$YEAR$MONTH$DAY/$HOUR/wave/gridded is generated 
#    containing the data, plus a log file logGEFS_$YEAR$MONTH$DAY$HOUR.
#
# Example:
#  This program can be run with: 
#    bash download_GEFSwaves.sh /home/name/data/nccf/com/gens/prod
#
# DEPENDENCIES:
#  wget, perl
#  For most linux systems, wget and perl are already included.
#  Two codes get_grib.pl and get_inv.pl must be included for the download.
#  If they are not in the directory, they will be automatically downloaded.
#
# AUTHOR and DATE:
#  11/09/2022: Ricardo M. Campos, first simplified version (previous 
#    get_gefsWaves.sh in WW3-Tools).
#  04/01/2024: Ricardo M. Campos, download of .grib2 field outputs 
#    without any post-processing.
#
# PERSON OF CONTACT:
#  Ricardo M. Campos: ricardo.campos@noaa.gov
#
########################################################################

# Path where this script is saved (or symbolic link), and where downloaded
#  files will be placed.
DIR="$1"
cd $DIR

# NOAA server address
SERVER=https://ftpprd.ncep.noaa.gov/
s1="global.0p25" # main grid

# variable names to be downloaded. Wind speed and significant wave height.
VARSGET=":WIND:surface:|:HTSGW:surface:"

# Initial date cycle for the ftp
YEAR=`date +%Y`
MONTH=`date +%m`
DAY=`date +%d`
# pa=2 #  days into the past. pa=1 downloads data from yesterday's cycle
# YEAR=`date --date=-$pa' day' '+%Y'`
# MONTH=`date --date=-$pa' day' '+%m'`
# DAY=`date --date=-$pa' day' '+%d'`
HOUR="00" # forecast cycle 00Z

DIRS="$DIR/gefs.$YEAR$MONTH$DAY/$HOUR/wave/gridded"
# create directory
mkdir -p $DIRS
# log file 
echo "  " > $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR
echo " Download GEFSv12 Ensemble Wave Forecast. " > $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR
echo " ====================================================== " > $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR
echo "  " > $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR

# Auxiliar pearl scripts, get_grib.pl get_inv.pl
# https://www.cpc.ncep.noaa.gov/products/wesley/fast_downloading_grib.html
if [ ! -f "get_inv.pl" ]; then
    # File doesn't exist, download it using wget
    wget --no-check-certificate --no-proxy -l1 -H -t1 -nd -N -np -erobots=off --tries=3 https://ftp.cpc.ncep.noaa.gov/wd51we/fast_downloading_grib/get_inv.pl
    chmod +x get_inv.pl
else
    echo "get_inv.pl already exists in the directory." > $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR
fi

if [ ! -f "get_grib.pl" ]; then
    # File doesn't exist, download it using wget
    wget --no-check-certificate --no-proxy -l1 -H -t1 -nd -N -np -erobots=off --tries=3 https://ftp.cpc.ncep.noaa.gov/wd51we/fast_downloading_grib/get_grib.pl
    chmod +x get_grib.pl
else
    echo "get_grib.pl already exists in the directory." > $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR
fi
# symbolic links
ln -fs $DIR/get_inv.pl $DIRS/
ln -fs $DIR/get_grib.pl $DIRS/

# hours of forecast lead time to be dowloaded. Using 6-h of temporal resolution
fleads="`seq -f "%03g" 0 6 384`"
# number of ensemble members
ensblm="`seq -f "%02g" 0 1 30`"

# Dowload loop
for h in $fleads;do
  for e in $ensblm;do

    echo "  " >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR
    echo " ======== GEFS Forecast: $YEAR$MONTH$DAY$HOUR  $h ========" >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR 

    if [ ${e} == 00 ]; then
      fname=$DIRS/gefs.wave.t00z.c${e}.${s1}.f"$(printf "%03.f" $h)".grib2
    else
      fname=$DIRS/gefs.wave.t00z.p${e}.${s1}.f"$(printf "%03.f" $h)".grib2
    fi

    RUN=0
    test -f ${fname}
    TE=$?
    if [ ${TE} -eq 1 ]; then
      RUN=1
    else
      FSIZE=$(du -sb "${fname}" | awk '{ print $1 }') >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
      if [ ${FSIZE} -lt 800000 ]; then
        RUN=1
      fi
    fi

    if [ ${RUN} -eq 1 ]; then
      FSIZE=0
      TRIES=1
      # while file has lower size than expected it does:
      while [ "${FSIZE}" -lt 800000 ] && [ "${TRIES}" -le 130 ]; do
        # sleep 5 minutes between attemps
        if [ ${TRIES} -gt 5 ]; then
          sleep 300
        fi

        if [ ${FSIZE} -lt 800000 ]; then
            echo "  " >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR
            echo " attempt number: $TRIES" >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR
            echo "  " >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR
            # main line where get_inv.pl and get_grib.pl are used to fech grib2 files
            if [ ${e} == 00 ]; then
               $DIRS/get_inv.pl $SERVER/data/nccf/com/gens/prod/gefs.$YEAR$MONTH$DAY/$HOUR/wave/gridded/gefs.wave.t${HOUR}z.c${e}.${s1}.f"$(printf "%03.f" $h)".grib2.idx | egrep "($VARSGET)" | $DIRS/get_grib.pl $SERVER/data/nccf/com/gens/prod/gefs.$YEAR$MONTH$DAY/$HOUR/wave/gridded/gefs.wave.t${HOUR}z.c${e}.${s1}.f"$(printf "%03.f" $h)".grib2 $DIRS/gefs.wave.t${HOUR}z.c${e}.${s1}.f"$(printf "%03.f" $h)".grib2 >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1 
               test -f $DIRS/gefs.wave.t${HOUR}z.c${e}.${s1}.f"$(printf "%03.f" $h)".grib2 >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1   
            else
               $DIRS/get_inv.pl $SERVER/data/nccf/com/gens/prod/gefs.$YEAR$MONTH$DAY/$HOUR/wave/gridded/gefs.wave.t${HOUR}z.p${e}.${s1}.f"$(printf "%03.f" $h)".grib2.idx | egrep "($VARSGET)" | $DIRS/get_grib.pl $SERVER/data/nccf/com/gens/prod/gefs.$YEAR$MONTH$DAY/$HOUR/wave/gridded/gefs.wave.t${HOUR}z.p${e}.${s1}.f"$(printf "%03.f" $h)".grib2 $DIRS/gefs.wave.t${HOUR}z.p${e}.${s1}.f"$(printf "%03.f" $h)".grib2 >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1 
               test -f $DIRS/gefs.wave.t${HOUR}z.p${e}.${s1}.f"$(printf "%03.f" $h)".grib2 >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
            fi
            # test if downloaded file exists
            TE=$?
            if [ ${TE} -eq 1 ]; then
              FSIZE=0
            else
              # check size of each file
              if [ ${e} == 00 ]; then
                FSIZE=`du -sb $DIRS/gefs.wave.t${HOUR}z.c${e}.${s1}.f"$(printf "%03.f" $h)".grib2 | awk '{ print $1 }'` >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
              else
                FSIZE=`du -sb $DIRS/gefs.wave.t${HOUR}z.p${e}.${s1}.f"$(printf "%03.f" $h)".grib2 | awk '{ print $1 }'` >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
              fi

            fi
        fi

        TRIES=`expr $TRIES + 1`
      done
    fi

  done
done

# permissions
# chmod -R 775 $DIRS/gefsWave.$YEAR$MONTH$DAY$HOUR
# 
echo "  " >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR
echo " Download complete. " >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR
echo "  " >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR
echo " ====================================================== " >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR
echo "  " >> $DIRS/logGEFS_$YEAR$MONTH$DAY$HOUR

