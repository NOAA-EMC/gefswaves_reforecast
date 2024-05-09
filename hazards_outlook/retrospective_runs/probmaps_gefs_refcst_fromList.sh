#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# This job run the code probmaps_gefs_refcst_fromList.sh on Orion
# Each job submission runs consecutive forecast cycles, with dates
#  specified in buoys_sel.txt
# All the information is specified in probmaps_gefs_refcst.yaml
#
# DIRJOUT="/work/noaa/marine/ricardo.campos/work/products/probmaps/reruns_experiments/version1p2/jobs"
# sbatch --output=${DIRJOUT}/probmaps_gefs_refcst_fromList.out probmaps_gefs_refcst_fromList.sh

ulimit -s unlimited
ulimit -c 0

echo " Starting at "$(date +"%T")

# INPUT ARGUMENTS
# .yaml configuration file containing paths and information for this
#   shell script as well as for the python code.
PYCYAML="/work/noaa/marine/ricardo.campos/work/products/probmaps/reruns_experiments/probmaps_gefs_refcst.yaml"
# Text file with list of events to reforecast, with date (YYYYMMDD) and string
LISTEVENTS="/work/noaa/marine/ricardo.campos/work/products/probmaps/reruns_experiments/list_events.txt"
# Cycle hour is fixed (only one per day)
HOUR=0

# Read the YAML as a text file:
#  GEFS data path
gefspath_line=$(grep 'gefspath' "${PYCYAML}")
GEFSMDIR=$(echo "$gefspath_line" | awk -F': ' '{print $2}')
#  Python script (probability maps)
pyscript_line=$(grep 'pyscript' "${PYCYAML}")
PYSCRIPT=$(echo "$pyscript_line" | awk -F': ' '{print $2}')
#  Python script (hindcast reference)
pyhindcst_line=$(grep 'pyhindcst' "${PYCYAML}")
PYHINDCST=$(echo "$pyhindcst_line" | awk -F': ' '{print $2}')
#  Variable names, for the python processing (probability maps)
mvars_line=$(grep 'mvars' "${PYCYAML}")
MVARS=$(echo "$mvars_line" | awk -F': ' '{gsub(/"/, "", $2); print $2}')
#  Output path
outpath_line=$(grep 'outpath' "${PYCYAML}")
OUTPATH=$(echo "$outpath_line" | awk -F': ' '{print $2}')

# READ LIST
fdates=()
fstrings=()
# Read the file line by line
while IFS= read -r line; do
  # Skip lines starting with #
  if [[ $line =~ ^# ]]; then
    continue
  fi
  # Extract date and string from each line
  fdate=$(echo "$line" | awk '{print $1}')
  fstring=$(echo "$line" | awk '{print $2}')
  # Add date and fstring to respective arrays
  fdates+=("$fdate")
  fstrings+=("$fstring")
done < ${LISTEVENTS}

# python env
source /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
sh /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh

echo "  "
echo " PYTHON PROCESSING: GLOBAL HAZARDS OUTLOOK - PROBABILITY MAPS"
echo "  "

# Iterate over dates array
for (( i=0; i<${#fdates[@]}; i++ )); do
  # Access each date and string using array indices
  echo "Date: ${fdates[i]}, String: ${fstrings[i]}"

  DATEF=${fdates[i]}
  DATE=${DATEF}$(printf "%02.f" $HOUR)

  # Check data
  test -f $GEFSMDIR"/gefs.wave."${DATEF}".00.global.0p25.nc"
  TE=$?
  if [ ${TE} -eq 1 ]; then
    echo " skip ${WW3VAR} ${DATE} no data"
  else
    # loop through variables
    for WW3VAR in ${MVARS[*]}; do
      python3 ${PYSCRIPT} ${PYCYAML} $DATE 7 14 ${WW3VAR} ${fstrings[i]}
      wait $!
      python3 ${PYHINDCST} ${PYCYAML} $DATE 7 14 ${WW3VAR} ${fstrings[i]}
      wait $!
      echo " Probability maps for ${WW3VAR} ${DATE} Ok."
    done
    wait $!
    sleep 1
    mkdir ${OUTPATH}/${fstrings[i]}
    mkdir ${OUTPATH}/${fstrings[i]}/Hs
    mkdir ${OUTPATH}/${fstrings[i]}/WS10
    mv ${OUTPATH}/ProbMap_Hs_*${fstrings[i]}*.png ${OUTPATH}/${fstrings[i]}/Hs
    mv ${OUTPATH}/HindcastReference_Hs_*${fstrings[i]}*.png ${OUTPATH}/${fstrings[i]}/Hs
    mv ${OUTPATH}/ProbMap_WS10_*${fstrings[i]}*.png ${OUTPATH}/${fstrings[i]}/WS10
    mv ${OUTPATH}/HindcastReference_WS10_*${fstrings[i]}*.png ${OUTPATH}/${fstrings[i]}/WS10
    wait $!
    sleep 1
  fi

done

echo "  "
echo " COMPLETE at "$(date +"%T")

