#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# This job run the code probmaps_gefs_refcst_monthly.sh on Orion
# Each job submission runs one month of consecutive forecast cycles
#
# DIRJOUT="/work/noaa/marine/ricardo.campos/work/products/probmaps/reruns_experiments/version1p2/jobs"
# export YEAR=2023
# export MONTH=1
# sbatch --output=${DIRJOUT}/probmaps_gefs_refcst_${YEAR}${MONTH}${DAY}.out probmaps_gefs_refcst.sh

ulimit -s unlimited
ulimit -c 0

export YEAR=${YEAR}
export MONTH=${MONTH}
HOUR=0

echo " Starting at "$(date +"%T")

# INPUT ARGUMENT
# .yaml configuration file containing paths and information for this
#   shell script as well as for the python code.
PYCYAML="/work/noaa/marine/ricardo.campos/work/products/probmaps/reruns_experiments/probmaps_gefs_refcst.yaml"
# Read the YAML as a text file:
#  GEFS data path
gefspath_line=$(grep 'gefspath' "${PYCYAML}")
GEFSMDIR=$(echo "$gefspath_line" | awk -F': ' '{print $2}')
#  Python script (probability maps)
pyscript_line=$(grep 'pyscript' "${PYCYAML}")
PYSCRIPT=$(echo "$pyscript_line" | awk -F': ' '{print $2}')
#  Variable names, for the python processing (probability maps)
mvars_line=$(grep 'mvars' "${PYCYAML}")
MVARS=$(echo "$mvars_line" | awk -F': ' '{gsub(/"/, "", $2); print $2}')
#  Output path
outpath_line=$(grep 'outpath' "${PYCYAML}")
OUTPATH=$(echo "$outpath_line" | awk -F': ' '{print $2}')

# python env
source /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
sh /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh

echo "  "
echo " PYTHON PROCESSING: GLOBAL HAZARDS OUTLOOK - PROBABILITY MAPS"
echo "  "

for DAY in `seq 1 31`; do

  DATEF=${YEAR}$(printf "%02.f" $MONTH)$(printf "%02.f" $DAY)
  DATE=${DATEF}$(printf "%02.f" $HOUR)

  # Check data
  test -f $GEFSMDIR"/gefs.wave."${DATEF}".00.global.0p25.nc"
  TE=$?
  if [ ${TE} -eq 1 ]; then
    echo " skip ${WW3VAR} ${DATE} no data"
  else
    # loop through variables
    for WW3VAR in ${MVARS[*]}; do
      python3 ${PYSCRIPT} ${PYCYAML} $DATE 7 14 ${WW3VAR}
      echo " Probability maps for ${WW3VAR} ${DATE} Ok."
      wait $!
    done
    wait $!
    sleep 1
  fi

done

echo "  "
echo " COMPLETE at "$(date +"%T")

