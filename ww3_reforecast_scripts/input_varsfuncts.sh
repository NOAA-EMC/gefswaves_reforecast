#!/bin/bash

# input_varsfuncts.sh
# 
# VERSION AND LAST UPDATE:
#  v1.0  11/04/2022
#
# PURPOSE:
#  Description of the main paths, global input variables, WW3 variables, and
#  Group of shell functions that will be used for the GEFSv12 reforecasts
# 
# AUTHOR and DATE:
#  11/04/2022: Ricardo M. Campos, first version.
#
# PERSON OF CONTACT:
#  Ricardo M Campos: ricardo.campos@noaa.gov
#

export DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/ww3runs/scripts"
export WW3MAIN="/work/noaa/marine/ricardo.campos/models/WW3/wave_ens_reforecast/sorc/WW3.fd/model"
export DIRF="/work/noaa/marine/ricardo.campos/work/ww3runs/forcing"
export DIRWI="/work/noaa/marine/ricardo.campos/data/GEFSv12"
export DIRMD="/work/noaa/marine/ricardo.campos/work/ww3runs/grids"
export DIRRS="/work/noaa/marine/ricardo.campos/work/ww3runs/restart"
export DIRWORK="/work/noaa/marine/ricardo.campos/work/ww3runs/work"
export DIROUT="/work/noaa/marine/ricardo.campos/work/ww3runs/output"
export DIRJOUT="/work/noaa/marine/ricardo.campos/work/ww3runs/jobs"
export DIRESULT="/work/noaa/marine/ricardo.campos/work/ww3runs/results"
export MDEFS="glo_15mxt glo_15m so_20m ao_20m points"
export WIGRID="glo_15mxt"
export OPOINTS="/work/noaa/marine/ricardo.campos/work/ww3runs/fix/bstations_GEFSv12WW3grid.txt"
export OUTPARS_WAV="WND HS FP T01 T02 DIR DP SPR PHS PTP PDIR"
export OUTPARS_SWAVG="WND HS FP T01 T02 DIR DP SPR PHS PTP PDIR"
# forecast restart date, how many days ahead:
export FRESTT="1"
# spinup run: no(0), yes(1)
export SPINUP="0"

# number of prep input processing per job
NIDAYS=${NIDAYS:=365}
# number of forecast cycles per job, main job ww3_multi run
NDAYS=${NDAYS:=7}
# number of output processing per job
NODAYS=${NODAYS:=20}

# Orion modules
module load contrib noaatools
module use /apps/contrib/NCEP/libs/hpc-stack/modulefiles/stack
module load hpc/1.1.0
module load hpc-intel/2018.4
module load hpc-impi/2018.4
module load jasper/2.0.25
module load zlib/1.2.11
module load png/1.6.35
module load bacio/2.4.1
module load g2/3.4.1
module load netcdf/4.7.2

# WAVEWATCH variables and paths
export WWATCH3_NETCDF=NC4
export NETCDF_CONFIG=/apps/intel-2018/netcdf-4.7.2/bin/nc-config
# 
PATH=$PATH:${WW3MAIN}/bin
export PATH
PATH=$PATH:${WW3MAIN}/exe
export PATH

# DEFINE FUNCTIONS ------------------------

# check inputs and preps, including wind and ice netcdfs, ww3 input binary files, and restart files. Check file and size.
function checkwprep()
{
  NDATE=$1
  ENSM=$2

  local cprep=0
  LMDEFS=($MDEFS)

  wday=$(date -d "$NDATE" +%u)
  if [ "$wday" -eq 3 ]; then
    frg=(10-35)
    lSW2=93000000
    lSWND=1200000000
  else
    frg=(10-16)
    lSW2=19000000
    lSWND=820000000
  fi

  # check netcdf wind and ice are ok
  test -f ${DIRWI}/${NDATE:0:4}/wnd10m_${NDATE}00_${ENSM}.D1-10.nc
  TE=$?
  if [ "$TE" -eq "1" ]; then
    SW1=0
  else
    SW1=`du -sbL ${DIRWI}/${NDATE:0:4}/wnd10m_${NDATE}00_${ENSM}.D1-10.nc | awk '{ print $1}'`
  fi

  test -f ${DIRWI}/${NDATE:0:4}/wnd10m_${NDATE}00_${ENSM}.D${frg}.nc
  TE=$?
  if [ "$TE" -eq "1" ]; then
    SW2=0
  else
    SW2=`du -sbL ${DIRWI}/${NDATE:0:4}/wnd10m_${NDATE}00_${ENSM}.D${frg}.nc | awk '{ print $1}'`
  fi

  if [ $SW1 -gt 240000000 ] && [ $SW2 -gt $lSW2 ]; then
    cprep=`expr $cprep + 1`
  fi

  # check ww3 binary wind and ice are ok
  test -f ${DIRF}/ice/ice.${NDATE}00.glo_15mxt
  TE=$?
  if [ "$TE" -eq "1" ]; then
    SIC=0
  else
    SIC=`du -sbL ${DIRF}/ice/ice.${NDATE}00.glo_15mxt | awk '{ print $1}'`
  fi

  test -f ${DIRF}/${ENSM}/wind.${NDATE}00.glo_15mxt
  TE=$?
  if [ "$TE" -eq "1" ]; then
    SWND=0
  else
    SWND=`du -sbL ${DIRF}/${ENSM}/wind.${NDATE}00.glo_15mxt | awk '{ print $1}'`
  fi 

  if [ $SIC -gt 3000000 ] && [ $SWND -gt $lSWND ]; then
    cprep=`expr $cprep + 1`
  fi

  # if cprep=2 everything is fine. Less than 2, it means preproc has problems for this date.
  if [ "$cprep" -eq "2" ]; then  
    echo "0"
  else
    echo "1"
  fi
  # same result number as test -f
}


# check inputs and preps, including wind and ice netcdfs, ww3 input binary files, and restart files. Check file and size.
function checkwic()
{
  NDATE=$1
  ENSM=$2

  local cric=0
  LMDEFS=($MDEFS)

  declare -a lSRES=("2300000000" "160000000" "160000000")

  # check restart for each grid
  i=0
  for grid in ${LMDEFS[@]:1:3}; do
    test -f ${DIRRS}"/stream"${STRM}"/"${ENSM}"/"${NDATE}".030000.restart."${grid}
    TE=$?
    if [ "$TE" -eq "1" ]; then
      SRES=0
    else
      SRES=`du -sbL ${DIRRS}/stream${STRM}/${ENSM}/${NDATE}.030000.restart.${grid} | awk '{ print $1}'`
    fi

    if [ $SRES -gt ${lSRES[$i]} ]; then
      cric=`expr $cric + 1`
    fi
    i=`expr $i + 1`
  done

  # if cric=3 everything is fine.
  if [ "$cric" -eq 3 ]; then  
    echo "0"
  else
    echo "1"
  fi
  # same result number as test -f
}


# check outputs, including out_grd for each grid, out_pnt, and restart files. Check file and size.
function checkwout()
{
  NDATE=$1
  STRM=$2
  ENSM=$3

  local cout=0
  LMDEFS=($MDEFS)

  wday=$(date -d "$NDATE" +%u)
  if [ "$wday" -eq 3 ]; then
    declare -a lSOUT=("11000000000" "910000000" "910000000" "3800000000")
  else
    declare -a lSOUT=("5700000000" "410000000" "400000000" "1600000000")
  fi
  declare -a lSRES=("2300000000" "160000000" "160000000")
  # restart time
  RTIME="$(date -u --date="${NDATE} +${FRESTT} day" '+%Y%m%d')"

  # check out_grd for each grid
  i=0
  for grid in ${LMDEFS[@]:1:3}; do
    test -f ${DIROUT}"/stream"${STRM}"/"${ENSM}"/out_grd."${NDATE}"."${grid}
    TE=$?
    if [ "$TE" -eq "1" ]; then
      SOUT=0
    else
      SOUT=`du -sbL ${DIROUT}/stream${STRM}/${ENSM}/out_grd.${NDATE}.${grid} | awk '{ print $1}'`
    fi

    if [ $SOUT -gt ${lSOUT[$i]} ]; then
      cout=`expr $cout + 1`
    fi
    i=`expr $i + 1`
  done

  # check out_pnt
  test -f ${DIROUT}"/stream"${STRM}"/"${ENSM}"/out_pnt."${NDATE}"."${LMDEFS[4]}
  TE=$?
  if [ "$TE" -eq "1" ]; then
    SOUT=0
  else
    SOUT=`du -sbL ${DIROUT}/stream${STRM}/${ENSM}/out_pnt.${NDATE}.${LMDEFS[4]} | awk '{ print $1}'`
  fi

  if [ $SOUT -gt ${lSOUT[$i]} ]; then
    cout=`expr $cout + 1`
  fi

  # check restart for each grid, next day
  i=0
  for grid in ${LMDEFS[@]:1:3}; do
    test -f ${DIRRS}"/stream"${STRM}"/"${ENSM}"/"${RTIME}".030000.restart."${grid}
    TE=$?
    if [ "$TE" -eq "1" ]; then
      SRES=0
    else
      SRES=`du -sbL ${DIRRS}/stream${STRM}/${ENSM}/${RTIME}.030000.restart.${grid} | awk '{ print $1}'`
    fi

    if [ $SRES -gt ${lSRES[$i]} ]; then
      cout=`expr $cout + 1`
    fi
    i=`expr $i + 1`
  done

  # if cout=7 everything is fine. Less than 7, it means preproc has problems for this date.
  if [ "$cout" -eq "7" ]; then  
    echo "0"
  else
    echo "1"
  fi
  # same result number as test -f
}


# for the extra ensemble members running on Wednesdays, pick up a restard id generated from the initial members c00 to p04.
function pickrestartid()
{
  NDATE=$1
  STRM=$2
  ENSM=$3
  #
  local restartid='999'

  file=${DIRSCRIPTS}"/ensMembIDlist.txt"
  while IFS= read -r line
  do

    YYYY=$(awk -F' ' '{print $1}' <<< "$line")
    MM=$(awk -F' ' '{print $2}' <<< "$line")
    DD=$(awk -F' ' '{print $3}' <<< "$line")
    if [ "$NDATE" -eq "$YYYY$MM$DD" ]; then
      aux="`expr $((10#${ENSM:1:3})) - 1`"
      restartid=$(awk -F' ' '{print $'$aux'}' <<< "$line")     
    fi

  done <"$file"

  echo "$restartid"
}


# check results. Final grib2 and netcdf files. Run a python script to open the files and verify one-by-one, creating a log. Plot figures.
function checkwresult()
{
  NDATE=$1
  STRM=$2
  ENSM=$3
  #
  DIRW=${DIRESULT}"/stream"${STRM}"/"${ENSM}

  local crsl=0

  export OUTFILET="_field.grib2 _spec.nc _tab.nc"
  LOUTFILET=($OUTFILET)

  wday=$(date -d "$NDATE" +%u)
  if [ "$wday" -eq 3 ]; then
    declare -a lSROUT=("3500000000" "800000000" "12000000")
  else
    declare -a lSROUT=("1400000000" "480000000" "8200000")
  fi

  # check restart for each grid, next day
  i=0
  for grid in ${LOUTFILET[@][1:3]}; do
    test -f ${DIRW}"/ww3gefs."${NDATE}${grid}
    TE=$?
    if [ "$TE" -eq "1" ]; then
      SRSL=0
    else
      SRSL=`du -sbL ${DIRW}/ww3gefs.${NDATE}${grid} | awk '{ print $1}'`
    fi

    if [ $SRSL -gt ${lSROUT[$i]} ]; then
      crsl=`expr $crsl + 1`
    fi
    i=`expr $i + 1`
  done

  # if crsl=3 everything is fine.
  if [ "$crsl" -eq 3 ]; then  
    echo "0"
  else
    echo "1"
  fi
  # same result number as test -f

}

# ----------

