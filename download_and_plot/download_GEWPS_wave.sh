#!/bin/bash

########################################################################
# download_GEWPS_wave.sh
#
# VERSION AND LAST UPDATE:
#   v1.0  06/18/2025
#
# PURPOSE:
#  Script to download Environment Canada's ensemble forecast
#   https://dd.meteo.gc.ca/
#   https://weather.gc.ca/ensemble/index_e.html
#   https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-2-10-0.shtml
#
# USAGE:
#  Two input arguments must be entered: cycle time (0,6,12,18) and the output path.
#  Example:
#    bash download_GEWPS_wave.sh 00 /home/ricardo/data/ECMWF_ENS
#  it will download the current day (00Z cycle)
#
# OUTPUT:
#  Multiple netcdf files containing the forecast steps, one file per ensemble member.
#
# DEPENDENCIES:
#  wget, CDO, NCO
#
# AUTHOR and DATE:
#  06/18/2025: Ricardo M. Campos, first version 
#
# PERSON OF CONTACT:
#  Ricardo M. Campos: ricardo.campos@noaa.gov
#
########################################################################

set -euo pipefail

# Function to check if file has been previously downloaded
checkfile() {
  local FILE="$1"

  if [[ -f "$FILE" ]]; then
    local TAM
    TAM=$(stat -c%s "$FILE" 2>/dev/null)

    if [[ "$TAM" -ge 10000000 ]]; then
      echo "0"  # File exists and is big enough
      return
    fi
  fi

  echo "1"  # File missing or too small
}

# Wave variables: Hs, Dp, Tp, Tm
WVARS=("HTSGW" "PWAVEDIR" "PWPER" "MZWPER")
# corresponding variable names after netcdf conversion
wvarname=("swh" "param46.0.10" "pp1d" "mp2")

# Ensemble members
ensbl=($(seq 1 21))
ensblm=($(seq -f "%02g" 0 1 20))

# Directory and cycle time
CHOUR="$1"
CHOUR=$(printf "%02.f" $CHOUR)
WDIR="$2"

cd $WDIR
if [ ! -d "$WDIR/work" ]; then
  mkdir "$WDIR/work"
fi
cd "$WDIR/work"

# EnvCanada server address
SERVER=https://dd.meteo.gc.ca
# Forecast lead time
fleads=$(echo $(seq -f "%03g" 0 3 168) $(seq -f "%03g" 174 6 384))

# Forecast cycle date
# DATE=$(date '+%Y%m%d')
pa=1
DATE=$(date --date="-${pa} day" '+%Y%m%d')

latmin=-82.
latmax=89.
dp=2

# Dowload loop
for h in $fleads;do

  for i in "${!WVARS[@]}"; do

    echo " ======== GEWPS Forecast: ${DATE} ${CHOUR}Z $h ${wv} ========"

    wv="${WVARS[$i]}"
    wvarn="${wvarname[$i]}"

    arqn="gewps.wave.${DATE}T${CHOUR}Z_"$(printf "%03.f" $h)"H_${wv}"

    FILE="$WDIR/GEWPS_wave_${DATE}${CHOUR}.00.nc"
    cfile=$(checkfile "$FILE")
    if [ "$cfile" -eq 0 ]; then
      echo "File $FILE already exists and is large enough. Skipping download 1."
      continue
    fi

    # size TAM and tries TRIES will control the process
    TAM=0
    TRIES=1
    while [ $TAM -lt 12000000 ] && [ $TRIES -le 130 ]; do
      # sleep 5 minutes between attemps
      if [ ${TRIES} -gt 5 ]; then
        sleep 30
      fi

      if [ ${TAM} -lt 12000000 ]; then

        wfile="${DATE}T${CHOUR}Z_MSC_GEWPS_${wv}_Sfc_LatLon0.25_PT"$(printf "%03.f" $h)"H.grib2"  
        wget -l1 -H -t1 -nd -N -np -erobots=off --tries=3 "${SERVER}/${DATE}/WXO-DD/model_gewps/25km/00/${wfile}" -O "$WDIR/work/${arqn}.grib2" 2>&1

        test -f "$WDIR/work/${arqn}.grib2"
        TE=$?
        if [ ${TE} -eq 1 ]; then
          TAM=0
        else
          # check size of each file
          TAM=`du -sb "$WDIR/work/${arqn}.grib2" | awk '{ print $1 }'`
        fi
        echo "wget $WDIR/work/${arqn}.grib2"

      fi
      TRIES=`expr $TRIES + 1`
    done

    test -f "$WDIR/work/${arqn}.grib2"
    TE=$?
    if [ ${TE} -eq 0 ]; then

      # wgrib2 "${arqn}.grib2" -netcdf ${arqn}.saux.nc 2>&1
      cdo -f nc4 copy "${arqn}.grib2" "${arqn}.saux.nc" 2>&1

      for j in "${!ensbl[@]}"; do
        e="${ensbl[$j]}"
        ne="${ensblm[$j]}"

        if [ "$j" -eq 0 ]; then
          ncks -v ${wvarn} ${arqn}.saux.nc -o ${arqn}_${ne}.saux.nc 2>&1
          cdo sellonlatbox,-180,180,-90,90 ${arqn}_${ne}.saux.nc ${arqn}_${ne}.saux1.nc 2>&1
          ncrename -v ${wvarn},${wv} ${arqn}_${ne}.saux1.nc ${arqn}_${ne}.nc 2>&1
        else
          ncks -v "${wvarn}_$e" ${arqn}.saux.nc -o ${arqn}_${ne}.saux.nc 2>&1
          cdo sellonlatbox,-180,180,-90,90 ${arqn}_${ne}.saux.nc ${arqn}_${ne}.saux1.nc 2>&1
          ncrename -v ${wvarn}_${e},${wv} ${arqn}_${ne}.saux1.nc ${arqn}_${ne}.nc  2>&1
        fi

        rm -f ${arqn}_${ne}.saux*.nc

      done

      rm -f $arqn.grib2
      rm -f $arqn.saux.nc

    fi

  done

  for i in "${!ensbl[@]}"; do
    ne="${ensblm[$i]}"

    FILE="$WDIR/GEWPS_wave_${DATE}${CHOUR}.${ne}.nc"
    cfile=$(checkfile "$FILE")
    if [ "$cfile" -eq 0 ]; then
      echo "File $FILE already exists and is large enough. Skipping download 2."
      continue
    fi

    for j in "${!WVARS[@]}"; do
      wv="${WVARS[$j]}"
      if [ "$j" -eq 0 ]; then
        iname="gewps.wave.${DATE}T${CHOUR}Z_$(printf "%03.f" "$h")H_${ne}.nc"
        mv "gewps.wave.${DATE}T${CHOUR}Z_$(printf "%03.f" "$h")H_${wv}_${ne}.nc" "$iname"
      else
        src="gewps.wave.${DATE}T${CHOUR}Z_$(printf "%03.f" "$h")H_${wv}_${ne}.nc"
        ncks -A -v $wv "$src" "$iname"  2>&1
        rm -f $src
      fi
    done

  done

done


for i in "${!ensbl[@]}"; do

  ne="${ensblm[$i]}"

  FILE="$WDIR/GEWPS_wave_${DATE}${CHOUR}.${ne}.nc"
  cfile=$(checkfile "$FILE")
  if [ "$cfile" -eq 0 ]; then
    echo "File $FILE already exists and is large enough. Skipping download 3."
    continue
  fi

  arqn="GEWPS_wave_"${DATE}${CHOUR}.${ne}
  ncrcat gewps.wave.${DATE}T${CHOUR}Z_*H_${ne}.nc ${arqn}.saux.nc  2>&1
  rm -f gewps.wave.${DATE}T${CHOUR}Z_*H_${ne}.nc

  ncks -4 -L 1 -d lat,${latmin},${latmax} ${arqn}".saux.nc" ${arqn}".saux1.nc" 2>&1
  ncks --ppc default=.$dp ${arqn}".saux1.nc" "${arqn}.nc" 2>&1
  ncatted -a _FillValue,,o,f,NaN "${arqn}.nc" 2>&1
  chmod 775 "${arqn}.nc"
  mv "${arqn}.nc" $WDIR

  rm -f ${arqn}.saux*
  rm -f ${arqn}.*idx*
  rm -f *ncks* *tmp
  echo " File ${arqn} converted to netcdf and compressed with success. "
  sleep 1

done


