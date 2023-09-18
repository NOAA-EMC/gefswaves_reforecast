#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128G
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# This job code runs extractPoint_makeHindcast.py on Orion
#
# DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/1preproc/extract_point_from_grib/jobs"
# export YEAR=2000
# export MONTH=1
# sbatch --output=${DIRJOUT}/jextractPoint_makeHindcast_${YEAR}_${MONTH}.out jextractPoint_makeHindcast.sh

ulimit -s unlimited
ulimit -c 0

export YEAR=${YEAR}
export MONTH=${MONTH}
DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/1preproc/extract_point/"

echo " Starting at "$(date +"%T")

# python env
source /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
sh /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
# work dir
cd ${DIRSCRIPTS}
# run
/work/noaa/marine/ricardo.campos/progs/python/anaconda3/bin/python3 ${DIRSCRIPTS}/extractPoint_makeHindcast.py ${YEAR} ${MONTH}

wait $!
echo " Complete at "$(date +"%T")

