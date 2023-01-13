#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# This job code runs checkww3result.py to verify the simulation is correct.
#
# DIRJOUT="/work/noaa/marine/ricardo.campos/work/ww3runs/jobs"
# export STRM="1"
# export YEAR="2000"
# export ENSM="c00"
# sbatch --output=${DIRJOUT}/runwcheck_${STRM}_${YEAR}_${ENSM}.out jrunwcheck.sh

ulimit -s unlimited
ulimit -c 0

# input variables and functions
source /work/noaa/marine/ricardo.campos/work/ww3runs/scripts/input_varsfuncts.sh
# python env
source /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
sh /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
# figures
opath="/work/noaa/marine/ricardo.campos/work/ww3runs/results/check/figs/"
# work directory
cd /work/noaa/marine/ricardo.campos/work/ww3runs/results/check

echo " Starting at "$(date +"%T")
echo " "

for month in `seq 1 12`; do
  for day in `seq 1 31`; do

      WTIME=$YEAR`printf %2.2d $month``printf %2.2d $day`
      wday=$(date -d "$WTIME" +%u)
      wday="${wday//[$'\t\r\n ']}"

      test -f ${DIRESULT}"/stream"${STRM}"/"${ENSM}"/ww3gefs."${WTIME}"_field.grib2"
      TE=$?
      if [ "$TE" -eq "0" ]; then
        test -f ${opath}"/"${YEAR}"/CheckWW3GEFS_"${WTIME}"_"${ENSM}".png"
        TE=$?
        if [ "$TE" -eq "1" ]; then
          /work/noaa/marine/ricardo.campos/progs/python/anaconda3/bin/python3 checkww3result.py ${WTIME} "stream"${STRM}"_"${ENSM}
          sleep 2
          wait $!
          echo "   ok "${DIRESULT}"/stream"${STRM}"/"${ENSM}"/ww3gefs."${WTIME}"_field.grib2"
        fi
      else
        if [  $((10#${ENSM:1:3})) -lt 5 ]; then
         echo "NO FILE "${DIRESULT}"/stream"${STRM}"/"${ENSM}"/ww3gefs."${WTIME}"_field.grib2"
        else
          if [ "$wday" -eq 3 ] ; then 
            echo "NO FILE "${DIRESULT}"/stream"${STRM}"/"${ENSM}"/ww3gefs."${WTIME}"_field.grib2"
          fi
        fi
      fi
  done
done

echo " "
echo " Complete at "$(date +"%T")

