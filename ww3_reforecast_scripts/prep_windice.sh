#!/bin/bash

#
# prep_windice.sh
# 
# VERSION AND LAST UPDATE:
#  v1.0  09/02/2022
#
# PURPOSE:
#  Prepare ice and wind binary files for each given cycle and member(wind) 
#
# USAGE:
#  Inputs: 
#  WW3MAIN="$1" # WW3/model dir
#  ENSM="$2" # ensemble member
#  DIRF: dir forcing /work/noaa/marine/ricardo.campos/work/ww3runs/forcing
#  DIRWI: dir where wind and ice inputs are saved
#  DIRMD: dir mod_def
#  WTIME: cycle time YYYYMMDD, always at 00Z so hour does not need to be informed
#  ------------------------------------------------------------
#  Example: 
#  WW3MAIN="/work/noaa/marine/ricardo.campos/models/WW3/model"
#  DIRF="/work/noaa/marine/ricardo.campos/work/ww3runs/forcing"
#  DIRWI="/work/noaa/marine/ricardo.campos/data/GEFSv12"
#  DIRMD="/work/noaa/marine/ricardo.campos/work/ww3runs/grids"
#  or export variables: source input_vars.sh
#  ./prep_windice.sh 20000101 c00 >> prep_windice.out 2>&1 &
#  -----------------------------------------------------------
#
# AUTHOR and DATE:
#  09/02/2022: Ricardo M. Campos, first version.
#
# PERSON OF CONTACT:
#  Ricardo M Campos: ricardo.campos@noaa.gov
#

# date
WTIME="$1"
# ensemble member
ENSM="$2"

source /work/noaa/marine/ricardo.campos/work/ww3runs/scripts/input_varsfuncts.sh 

WW3EXE="${WW3MAIN}/exe"
DIRFO=${DIRF}/${ENSM}

YEAR="$(date -d "$WTIME" '+%Y')"
ITIME="${WTIME}00"
# forecast range
wday=$(date -d "$WTIME" +%u)
if [ "$wday" -eq 3 ]; then
  frg=(10-35)
else
  frg=(10-16)
fi

# ========== ICE =============
if [ "${ENSM}" = "c00" ]; then
  cd $DIRF/ice
  # check mode_def files
  test -f mod_def.ww3; TE=$?
  if [ "$TE" -eq "1" ]; then
    ln -fs ${DIRMD}/mod_def.${WIGRID} mod_def.ww3
  fi
  # runs ww3_prnc only if file does not exist
  test -f ${DIRF}/ice/ice.${ITIME}.${WIGRID}
  TE=$?
  if [ "$TE" -eq "1" ]; then

cat > ww3_prnc.inp << EOF
$
 'ICE' 'LL' T T
 lon lat time
 icecsfc
 'ice.${ITIME}.nc'
$
EOF

   ln -fs ${DIRWI}/${YEAR}/ice.${ITIME}.nc
   srun -n 1 ${WW3EXE}/ww3_prnc
   wait $!
   sleep 2
   rm ww3_prnc.inp
   rm times.ICE
   mv ice.ww3 ice.${ITIME}.${WIGRID}
   rm ice.${ITIME}.nc
   wait $!
   sleep 2
  fi
fi

# =========== WIND ===============
cd $DIRFO
# check mode_def files
test -f mod_def.ww3; TE=$?
if [ "$TE" -eq "1" ]; then
  ln -fs ${DIRMD}/mod_def.${WIGRID} mod_def.ww3
fi
# runs ww3_prnc only if file does not exist
test -f ${DIRFO}/wind.${ITIME}.${WIGRID}
TE=$?
if [ "$TE" -eq "1" ]; then

# LEG1
cat > ww3_prnc.inp << EOF
$
 'WND' 'LL' T T
 longitude latitude time
 uwnd vwnd 
 'wnd10m_${ITIME}_${ENSM}.D1-10.nc'
$
EOF

 ln -fs ${DIRWI}/${YEAR}/wnd10m_${ITIME}_${ENSM}.D1-10.nc
 wait $!
 sleep 2
 srun -n 1 ${WW3EXE}/ww3_prnc
 wait $!
 sleep 2
 rm ww3_prnc.inp
 rm times.WND
 mv wind.ww3 wind.${ITIME}.leg1.ww3
 rm wnd10m_${ITIME}_${ENSM}.D1-10.nc
 wait $!
 sleep 2

# LEG2
cat > ww3_prnc.inp << EOF
$
 'WND' 'LL' T F
 longitude latitude time
 uwnd vwnd 
 'wnd10m_${ITIME}_${ENSM}.D${frg}.nc'
$
EOF

 ln -fs ${DIRWI}/${YEAR}/wnd10m_${ITIME}_${ENSM}.D${frg}.nc
 wait $!
 sleep 2
 srun -n 1 ${WW3EXE}/ww3_prnc
 wait $!
 sleep 2
 rm ww3_prnc.inp
 rm times.WND
 mv wind.ww3 wind.${ITIME}.leg2.ww3
 rm wnd10m_${ITIME}_${ENSM}.D${frg}.nc
 wait $!
 sleep 2
 cat wind.${ITIME}.leg1.ww3 wind.${ITIME}.leg2.ww3 > wind.${ITIME}.${WIGRID}
 wait $!
 sleep 2
 rm wind.${ITIME}.leg1.ww3 wind.${ITIME}.leg2.ww3

fi


