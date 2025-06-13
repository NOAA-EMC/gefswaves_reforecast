#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128G
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# This job script run the code extract_GEFS.py on Orion for one month
#
# DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/jobs"
# export YEAR=2023
# export MONTH=1
# sbatch --output=${DIRJOUT}/jextract_GEFS_${YEAR}${MONTH}.out jextract_GEFS.sh

echo " Starting at "$(date +"%T")

ulimit -s unlimited
ulimit -c 0

DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification"
DIRO="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/data" # final netcdf4 output path

export YEAR=${YEAR}
export MONTH=${MONTH}

# python env
source /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
sh /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh

echo "  "
echo " EXTRACT GEFS Points"
echo "  "

# work dir
cd ${DIRSCRIPTS}

for DAY in `seq 1 31`; do
  for HOUR in 00 12; do
    export MONTH=${MONTH}
    export DAY=${DAY}
    export HOUR=${HOUR}

    DATE=${YEAR}$(printf "%02d" $MONTH)$(printf "%02d" $DAY)$(printf "%02d" $HOUR)
    export DATE=${DATE}

    FILE=$DIRO/GEFS.PointExtract.${DATE}.nc
    # Skip if file exists and is large enough
    if [ -f "$FILE" ]; then
      TAM=$(du -sb "$FILE" | awk '{ print $1 }')
      if [ "$TAM" -ge 400000000 ]; then
        echo "File $FILE already exists and is large enough. Skipping processing."
        continue
      else
        rm -f ${FILE}
      fi
    fi

    echo " "
    echo " Running for ${DATE}"
    echo " "

    python3 ${DIRSCRIPTS}/extract_GEFS.py ${DATE}
    wait $!
    echo " extract_GEFS.py for ${DATE} OK at $(date +"%T")"
    echo " "
    sleep 1
  done
done

echo " Complete at "$(date +"%T")

