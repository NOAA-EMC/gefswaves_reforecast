#!/bin/bash

########################################################################
# download_ECMWF_ENS_wind.sh
#
# VERSION AND LAST UPDATE:
#   v1.0  06/13/2025
#
# PURPOSE:
#  Script to download ECMWF Ensemble Forecast Data from OpenData
#  https://www.ecmwf.int/en/forecasts/datasets/open-data
#  using python https://pypi.org/project/ecmwf-opendata/
#  https://www.ecmwf.int/sites/default/files/elibrary/2021/81274-ifs-documentation-cy47r3-part-vii-ecmwf-wave-model_1.pdf
#
# USAGE:
#  Two input arguments must be entered: cycle time (0,6,12,18) and the output path.
#  Example:
#    bash download_ECMWF_ENS_wind.sh 00 /home/ricardo/data/ECMWF_ENS
#  it will download the current day (00Z cycle)
#
# OUTPUT:
#  Multiple netcdf files containing the forecast steps, one file per ensemble member.
#
# DEPENDENCIES:
#  wget, python (ecmwf.opendata), CDO, NCO
#  user must activate python/environment (see source below)
#
# AUTHOR and DATE:
#  06/13/2025: Ricardo M. Campos, first version 
#
# PERSON OF CONTACT:
#  Ricardo M. Campos: ricardo.campos@noaa.gov
#
########################################################################

set -euo pipefail

# Directory, cycle time, and python environment setup
CHOUR="$1"
CHOUR=$(printf "%02.f" $CHOUR)
DIR="$2"
cd "$DIR"
source /home/ricardo/work/python/anaconda3/setanaconda3.sh

# Function to convert and compress GRIB2 to NetCDF
function ccompress() {
  local arqn=$1
  local DIR=$2
  cd "$DIR"

  cdo -f nc4 copy "${arqn}.grib2" "${arqn}.saux1.nc"
  ncks -4 -L 1 -d lat,${latmin},${latmax} "${arqn}.saux1.nc" "${arqn}.saux2.nc"
  ncks -C -O -x -v height "${arqn}.saux2.nc" "${arqn}.saux3.nc"
  ncks --ppc default=.${dp} "${arqn}.saux3.nc" "${arqn}.nc"

  rm -f "${arqn}".saux* "${arqn}".*idx* "${arqn}".*ncks* "${arqn}".*tmp "${arqn}.grib2" "${sname}"

  echo "File ${arqn} converted to NetCDF and compressed."
}

# Forecast date
# DATE=$(date '+%Y%m%d')
pa=1
DATE=$(date --date="-${pa} day" '+%Y%m%d')

# Parameters
wparam_wind="10u/10v/msl"
latmin=-82.
latmax=89.
# decimals/resolution for compression
dp=2

# Lead times: 3h to 144h, then 6h to 360h
fleads=$(echo $(seq -f "%g" 0 3 144) $(seq -f "%g" 150 6 360) | tr ' ' '/')

# Ensemble members
ensblm=$(seq -f "%02g" 0 1 50)

# Loop over members
for e in $ensblm; do
  member_num=$((10#$e))
  arqn="ECMWF_ENS_wave_${DATE}${CHOUR}.${e}"
  FILE="${DIR}/${arqn}.nc"

  # Skip if existing file is large enough
  if [[ -f "$FILE" ]]; then
    size=$(du -sb "$FILE" | awk '{print $1}')
    if (( size >= 260000000 )); then
      echo " $FILE exists and is large enough. Skipping."
      continue
    else
      rm -f "$FILE"
    fi
  fi

  sname="download_ECMWF_ENS_wave_m${e}.py"
  echo " Running ${sname} ..."

  if (( member_num == 0 )); then
    etype="cf"
    number_line=""
  else
    etype="pf"
    number_line="        \"number\": \"${member_num}\","
  fi

  # Create download script
  cat > "${sname}" << EOF
#!/usr/bin/env python3
from ecmwf.opendata import Client
client = Client()

client.retrieve({
    "class": "od",
    "date": "${DATE}",
    "time": "${CHOUR}",
    "stream": "enfo",
    "type": "${etype}",
${number_line}
    "step": "${fleads}",
    "levtype": "sfc",
    "param": "${wparam_wind}",
    "target": "${arqn}.grib2"
})
EOF

  rm -f "${arqn}".saux* "${arqn}".*idx* "${arqn}".*ncks* "${arqn}".*tmp "${arqn}.grib2"
  chmod +x "${sname}"
  python3 "${sname}"

  ccompress "${arqn}" "${DIR}"
  chmod 775 "$FILE"
done

