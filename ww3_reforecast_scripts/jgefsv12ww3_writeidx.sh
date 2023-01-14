#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# This job code uses wgrib2 to write .idx files
#
# DIRJOUT="/work/noaa/marine/ricardo.campos/work/ww3runs/jobs"
# export STRM="1"
# export YEAR="2000"
# export ENSM="c00"
# sbatch --output=${DIRJOUT}/jgefsv12ww3_writeidx_${STRM}_${YEAR}_${ENSM}.out jgefsv12ww3_writeidx.sh

ulimit -s unlimited
ulimit -c 0

export STRM="${STRM}"
export YEAR="${YEAR}"
export ENSM="${ENSM}"

# input variables and functions
source /work/noaa/marine/ricardo.campos/work/ww3runs/scripts/input_varsfuncts.sh

module load wgrib2
# files
cd ${DIRESULT}"/stream"${STRM}"/"${ENSM}"/"

echo " Starting at "$(date +"%T")
echo " "

for month in `seq 1 12`; do
  for day in `seq 1 31`; do

    echo " "
    WTIME=$YEAR`printf %2.2d $month``printf %2.2d $day`

    test -f ${DIRESULT}"/stream"${STRM}"/"${ENSM}"/gefs.wave."${WTIME}"."${ENSM}".global.0p25.grib2" 
    TE=$?
    if [ "$TE" -eq "0" ]; then
   
          wgrib2 -s ${DIRESULT}"/stream"${STRM}"/"${ENSM}"/gefs.wave."${WTIME}"."${ENSM}".global.0p25.grib2" > ${DIRESULT}"/stream"${STRM}"/"${ENSM}"/gefs.wave."${WTIME}"."${ENSM}".global.0p25.idx"
          wait $!
          sleep 2
          echo " OK "${DIRESULT}"/stream"${STRM}"/"${ENSM}"/gefs.wave."${WTIME}"."${ENSM}

    else
      echo "NO FILE "${DIRESULT}"/stream"${STRM}"/"${ENSM}"/gefs.wave."${WTIME}"."${ENSM}".global.0p25.grib2"
    fi
    echo " "
  done
done

echo " "
echo " Complete at "$(date +"%T")

