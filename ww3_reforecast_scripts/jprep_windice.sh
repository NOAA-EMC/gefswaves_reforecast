#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion
#SBATCH --exclusive

# Prepare ice and wind binary files for each given cycle and member(wind)
#
# DIRJOUT="/work/noaa/marine/ricardo.campos/work/ww3runs/jobs"
# export STIME="20000201"
# export ENSM="c00"
# sbatch --output=${DIRJOUT}/jprep_windice_${STIME}_${ENSM}.out jprep_windice.sh

ulimit -s unlimited
ulimit -c 0

# input variables and functions
source /work/noaa/marine/ricardo.campos/work/ww3runs/scripts/input_varsfuncts.sh

for iday in $(seq $NIDAYS); do

  NDATE=$(date -u --date="${STIME} +$((${iday}-1)) day" '+%Y%m%d')

  if [ $((10#${ENSM:1:3})) -gt 4 -a $(date -d "$NDATE" +%u) -ne 3 ]; then
    echo " "
    echo " Not Wednesday ${NDATE}. Extra ensemble member ${ENSM} only available on Wednesdays. Moving to the next one ..."
    echo " "
  else

    sh ${DIRSCRIPTS}/prep_windice.sh ${NDATE} ${ENSM} &&
    wait $!
    sleep 2

  fi

done

