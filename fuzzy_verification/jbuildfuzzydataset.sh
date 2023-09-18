#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -q batch
#SBATCH -t 05:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# This job run the code buildfuzzydataset.py on Orion
#
# DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/jobs"
# export YEAR=2023
# export MONTH=1
# export DAY=15
# sbatch --output=${DIRJOUT}/jbuildfuzzydataset_${YEAR}${MONTH}${DAY}.out jbuildfuzzydataset.sh

ulimit -s unlimited
ulimit -c 0

export YEAR=${YEAR}
export MONTH=${MONTH}
export DAY=${DAY}

DATE=${YEAR}$(printf "%02.f" $MONTH)$(printf "%02.f" $DAY)
export DATE=${YEAR}$(printf "%02.f" $MONTH)$(printf "%02.f" $DAY)

DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification"

echo " Starting at "$(date +"%T")

# python env
source /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
sh /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
# work dir
cd ${DIRSCRIPTS}
# run
/work/noaa/marine/ricardo.campos/progs/python/anaconda3/bin/python3 ${DIRSCRIPTS}/buildfuzzydataset.py ${DATE}

wait $!
echo " Complete at "$(date +"%T")

