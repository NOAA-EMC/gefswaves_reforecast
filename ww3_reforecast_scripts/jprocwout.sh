#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -q batch
#SBATCH -t 07:50:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# Process output files, fields and points, for each given cycle and member(wind) and stream
# using ww3_ounp for spectra and tables, and ww3_grib for grid fields
#
# DIRJOUT="/work/noaa/marine/ricardo.campos/work/ww3runs/jobs"
# export STIME="20000108"
# export STRM="1"
# export ENSM="c00"
# sbatch --output=${DIRJOUT}/jprocwout_stream${STRM}.${ENSM}.out jprocwout.sh

ulimit -s unlimited
ulimit -c 0

# input variables and functions
source /work/noaa/marine/ricardo.campos/work/ww3runs/scripts/input_varsfuncts.sh

# paths
WW3EXE=${WW3MAIN}"/exe"
DIRW=${DIRESULT}"/stream"${STRM}"/"${ENSM}
FDIROUT=${DIROUT}"/stream"${STRM}"/"${ENSM}
LMDEFS=($MDEFS)

# if [ "$(ls $DIRW/*.grib2 | wc -l)" -gt 0 ]; then
#  # date time of the last result
#  LASTFTIME=$( ls ${DIRW}/ | grep ww3gefs.*_field.grib2 | tail -n 1 | awk -F'[_]' '{print $1}' | awk -F'[.]' '{print $2}')
#  # next cycle
#  STIME=$(date -u --date="${LASTFTIME} +1 day" '+%Y%m%d')
# fi

if [ -z "$STIME" ]; then
  echo " ERROR - STIME not assigned "
  exit 0
fi
export STIME="$STIME"

cd ${DIRW}
for iday in $(seq $NODAYS); do

  NDATE=$(date -u --date="${STIME} +$((${iday}-1)) day" '+%Y%m%d')
  YEAR="$(date -d "$NDATE" '+%Y')"
  wday=$(date -d "$NDATE" +%u)

  if [  $((10#${ENSM:1:3})) -gt 4 ] && [ "$wday" -ne "3" ]; then
    worun=0
    echo " No ww3 output binary files expected for this date and member, "${NDATE}" "${STRM}" "${ENSM}
  else
    worun=1
  fi
  if [ "$worun" -eq "1" ]; then
    cout=$(checkwout ${NDATE} ${STRM} ${ENSM})
    echo " "
    if [ "$cout" -eq "0" ]; then
      echo " OK - Expected binary ww3 outputs "${NDATE}" "${STRM}" "${ENSM}
    else
      echo " ERROR - missing ww3 output binary files, "${NDATE}" "${STRM}" "${ENSM}
      echo " ERROR - missing ww3 output binary files, "${NDATE}" "${STRM}" "${ENSM} >> ${DIRW}/ERROR_${STIME}_stream${STRM}_${ENSM}.txt
      exit 0
    fi

    # forecast range and high/low resolution time steps (high res for point output, low res for fields)
    if [ "$wday" -eq "3" ]; then
      fr=35
      hrts=840
      lrts=280
    else
      fr=16
      hrts=384
      lrts=128
    fi

    # MOD_DEFS
    for grid in ${LMDEFS[@]}; do
      test -f ${DIRW}/mod_def.${grid}; TE=$?
      if [ "$TE" -eq "1" ]; then
        test -f ${DIRMD}/mod_def.${grid}; TE=$?
        if [ "$TE" -eq "0" ]; then
          ln -fs ${DIRMD}/mod_def.${grid} ${DIRW}/mod_def.${grid}
        else
          echo "${DIRMD}/mod_def.${grid} file not found."; exit 10
        fi
      fi
    done
    # Output files
    for grid in ${LMDEFS[@]}; do
      test -f ${DIRW}/out_grd.${grid}
      TE=$?
      if [ "$TE" -eq "0" ]; then
        rm -f ${DIRW}/out_grd.${grid}
      fi
      test -f ${FDIROUT}/out_grd.${NDATE}.${grid}
      TE=$?
      if [ "$TE" -eq "0" ]; then
        ln -fs ${FDIROUT}/out_grd.${NDATE}.${grid} ${DIRW}/out_grd.${grid}
      fi
    done

    # POINT OUTPUTS ----------------------
    # TABLES
    ln -fs ${DIRMD}/mod_def.${LMDEFS[4]} ${DIRW}/mod_def.ww3
    ln -fs ${FDIROUT}/out_pnt.${NDATE}.${LMDEFS[4]} ${DIRW}/out_pnt.ww3

  # ww3_ounp
cat > ww3_ounp.inp << EOF
$ 
 ${NDATE} 030000 3600. ${hrts}
$
 -1
$
 ww3gefs.
 0
 4
 T 1 
 2
 0
 T
$
 2
$
EOF

    srun -n 1 ${WW3EXE}/ww3_ounp
    wait $!   
    mv ${DIRW}/ww3gefs.tab.nc ${DIRW}/ww3gefs.${NDATE}_tab.nc
    rm -f ${DIRW}/ww3_ounp.inp

  # SPECTRA
  # ww3_ounp
cat > ww3_ounp.inp << EOF
$ 
 ${NDATE} 030000 10800. ${lrts}
$
 -1
$
ww3gefs.
 0
 4
 T 1000
 1
 0
 T
 3 1. 0. 4.
$
EOF

    srun -n 1 ${WW3EXE}/ww3_ounp
    wait $!
    mv ${DIRW}/ww3gefs.spec.nc ${DIRW}/ww3gefs.${NDATE}_spec.nc
    rm -f ${DIRW}/ww3_ounp.inp
    rm -f ${DIRW}/mod_def.ww3
    rm -f ${DIRW}/out_pnt.ww3
    # --------------

    # FIELD OUTPUTS
    # gint
cat > ww3_gint.inp << EOF
$
 ${NDATE} 030000 10800. ${lrts}
$
 4
$
 '${LMDEFS[1]}'
 '${LMDEFS[2]}'
 '${LMDEFS[3]}'
 '${LMDEFS[0]}'
$
 0
$
EOF

    srun -n 1 ${WW3EXE}/ww3_gint
    wait $!
    rm -f ${DIRW}/ww3_gint.inp

    # ww3_grib
    mv ${DIRW}/out_grd.glo_15mxt ${DIRW}/out_grd.ww3
    ln -fs ${DIRMD}/mod_def.${LMDEFS[0]} ${DIRW}/mod_def.ww3

cat > ww3_grib.inp << EOF
$  
${NDATE} 030000 10800 ${lrts}
N
${OUTPARS_SWAVG}
$
${NDATE} 030000 7 11 255 0 0
$
EOF

    srun -n 1 ${WW3EXE}/ww3_grib

    mv ${DIRW}/gribfile ${DIRW}/ww3gefs.${NDATE}_field.grib2
    rm -f ${DIRW}/ww3_grib.inp
    rm -f ${DIRW}/mod_def.*
    rm -f ${DIRW}/out_grd.*
    rm -f ${DIRW}/WHTGRIDINT.bin
    wait $!
    sleep 5

    crsl=$(checkwresult ${NDATE} ${STRM} ${ENSM})
    if [ "$crsl" -eq "0" ]; then
      echo " OK - grib2 and netcdf files have been successfully generated, "${NDATE}" "${STRM}" "${ENSM}
    else
      echo " ERROR - grib2 and netcdf files, "${NDATE}" "${STRM}" "${ENSM}
      echo " ERROR - grib2 and netcdf files, "${NDATE}" "${STRM}" "${ENSM} >> ${DIRW}/ERROR_${STIME}_stream${STRM}_${ENSM}.txt
      exit 0
    fi
    echo " "

  fi
done

