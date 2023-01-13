#!/bin/bash

#
# prep_multi.sh
# 
# VERSION AND LAST UPDATE:
#  v1.0  08/28/2022
#
# PURPOSE:
#  Organize mod_defs, wind and ice binary inputs, and restart files, 
#    for each given cycle and member(wind) and stream.
#
# USAGE:
#  Inputs:
#  WW3MAIN: WW3/model dir
#  WTIME: forecast cycle date YYYYMMDD
#  STRM: stream
#  ENSM:  ensemble member
#  ------------------------------------------------------------
#  Example: 
#  export variables: source input_vars.sh
#  ./prep_multi.sh 20000101 1 c00 >> prep_multi.out 2>&1 &
#  -----------------------------------------------------------
#
# AUTHOR and DATE:
#  08/28/2022: Ricardo M. Campos, first version.
#
# PERSON OF CONTACT:
#  Ricardo M Campos: ricardo.campos@noaa.gov
#

# forecast cycle date YYYYMMDD
WTIME="$1"
# stream, 1 to 4
STRM="$2"
# ensemble member, c00, p01, p02 ...
ENSM="$3"
# input variables
source /work/noaa/marine/ricardo.campos/work/ww3runs/scripts/input_varsfuncts.sh 
DIRW="${DIRWORK}/stream${STRM}/${ENSM}"
DIRR="${DIRRS}/stream${STRM}/${ENSM}"
LMDEFS=($MDEFS)
cd "${DIRW}"
# clean directory and prepare for the next run
rm -f ${DIRW}/ice.*
rm -f ${DIRW}/wind.*
rm -f ${DIRW}/ww3_multi.inp
rm -f ${DIRW}/log.*
rm -f ${DIRW}/rmp_src_to_dst_conserv_*
rm -f ${DIRW}/out_*
rm -f ${DIRW}/*restart*

# Time Organizing
ITIME="${WTIME}00"
# forecast range
if [ "$SPINUP" -eq "1" ]; then
  fr=1
  FHOUR="03"
else
  FHOUR="00"
  wday=$(date -d "$WTIME" +%u)
  if [ "$wday" -eq "3" ]; then
    fr=35
  else
    fr=16
  fi
fi

# final time
FTIME="$(date -u --date="${WTIME} +${fr} day" '+%Y%m%d')"
# restart time
RTIME="$(date -u --date="${WTIME} +${FRESTT} day" '+%Y%m%d')"

# MOD_DEFS
for mdef in ${MDEFS[*]}; do
  test -f ${DIRW}/mod_def.${mdef}; TE=$?
  if [ "$TE" -eq "1" ]; then
    test -f ${DIRMD}/mod_def.${mdef}; TE=$? 
    if [ "$TE" -eq "0" ]; then
      ln -fs ${DIRMD}/mod_def.${mdef} ${DIRW}/mod_def.${mdef}
    else
      echo "${DIRMD}/mod_def.${mdef} file not found."; exit 10
    fi
  fi
done

# INPUTS: wind and ice
test -f ${DIRF}/${ENSM}/wind.${ITIME}.${LMDEFS[0]}
TE=$?
if [ "$TE" -eq "0" ]; then
  ln -fs ${DIRF}/${ENSM}/wind.${ITIME}.${LMDEFS[0]} ${DIRW}/wind.${LMDEFS[0]}
else
  echo "${DIRF}/${ENSM}/wind.${ITIME}.${LMDEFS[0]} file not found."; exit 10
fi

test -f ${DIRF}/ice/ice.${ITIME}.${LMDEFS[0]}
TE=$?
if [ "$TE" -eq "0" ]; then
  ln -fs ${DIRF}/ice/ice.${ITIME}.${LMDEFS[0]} ${DIRW}/ice.${LMDEFS[0]}
else
  echo "${DIRF}/ice/ice.${ITIME}.${LMDEFS[0]} file not found."; exit 10
fi

# RESTART FILE
for grid in ${MDEFS[*]}; do
  test -f ${DIRW}/restart.${grid}; TE=$?
  if [ "$TE" -eq "0" ]; then
    rm -f ${DIRW}/restart.${grid}
  fi
  test -f ${DIRR}/${WTIME}.030000.restart.${grid}; TE=$?
  if [ "$TE" -eq "0" ]; then
    ln -fs ${DIRR}/${WTIME}.030000.restart.${grid} restart.${grid}
  fi
done

# ww3_multi
cat > ww3_multi.inp << EOF
$
 3 1 T 3 T F
$
$
 '${LMDEFS[0]}' F F T T F F F F F
$
 '${LMDEFS[4]}'
$
$ the presence of 1) water levels 2) currents 3) winds 4) ice and
$ 5-7) assimilation data as in the file ww3_shel.inp.
 '${LMDEFS[1]}'  'no' 'no' '${LMDEFS[0]}' '${LMDEFS[0]}' 'no' 'no' 'no'  'no' 'no'  1  10  0.00 1.00  F
 '${LMDEFS[2]}'  'no' 'no' '${LMDEFS[0]}' '${LMDEFS[0]}' 'no' 'no' 'no'  'no' 'no'  2  20  0.00 1.00  F
 '${LMDEFS[3]}'  'no' 'no' '${LMDEFS[0]}' '${LMDEFS[0]}' 'no' 'no' 'no'  'no' 'no'  3  30  0.00 1.00  F
$
 ${WTIME} 030000   ${FTIME} ${FHOUR}0000
$
 T T
$
 ${WTIME} 030000  10800  ${FTIME} ${FHOUR}0000
 N
 ${OUTPARS_WAV}
$
 ${WTIME} 030000  3600  ${FTIME} ${FHOUR}0000
$
$    LON       LAT       NAME     
$ ----------------------------
EOF

cat ${OPOINTS} >> ww3_multi.inp

cat >> ww3_multi.inp << EOF
$
 0.00    0.00  'STOPSTRING'  999.   XXX  NCEP       0
$
 ${WTIME} 030000  0  ${FTIME} ${FHOUR}0000
$
 ${WTIME} 030000  0  ${FTIME} ${FHOUR}0000 T
 ${RTIME} 030000  1  ${RTIME} 030000
$
 ${WTIME} 030000  0  ${FTIME} ${FHOUR}0000
$
 ${WTIME} 030000  0  ${FTIME} ${FHOUR}0000
$
 'the_end'  0
$
 'STP'
$
EOF


