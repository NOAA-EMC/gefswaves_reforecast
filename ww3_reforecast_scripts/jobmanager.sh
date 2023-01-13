#!/bin/bash

#
# jobmanager.sh
# 
# VERSION AND LAST UPDATE:
#  v1.0  10/12/2022
#
# PURPOSE:
#  Main script to run WW3 GEFSv12 reforecast (automatically submit the jobs),
#   calling other functions and scripts. 
#  It works as a main manager of the reforecast run.
#
# USAGE:
#  Two input arguments: stream and ensemble member.
#  nohup ./jobmanager.sh 1 c00 >> logmanager/jobmanager_1_c00.out 2>&1 &
#
# AUTHOR and DATE:
#  10/12/2022: Ricardo M. Campos, first version.
#
# PERSON OF CONTACT:
#  Ricardo M Campos: ricardo.campos@noaa.gov
# 

# input stream
STRM="$1"
# input ensemble member
ENSM="$2"
# consecutive jobs to be submitted.
cjobs=3

export STRM="$STRM"
export ENSM="$ENSM"

# input variables and functions
source /work/noaa/marine/ricardo.campos/work/ww3runs/scripts/input_varsfuncts.sh 
# work directory
DIRW="${DIRWORK}/stream${STRM}/${ENSM}"
# output directory
DIRO="${DIROUT}/stream${STRM}/${ENSM}"
# grids
LMDEFS=($MDEFS)
echo " Starting at "$(date +"%T")
echo " "

for ijw in $(seq $cjobs); do

  echo "  "
  # date time of the last output file
  LASTFTIME=$( ls ${DIRO}/ | grep out_grd.*${LMDEFS[1]} | tail -n 1 | awk -F'[.]' '{print $2}' )
  
  # next cycle
  if [  $((10#${ENSM:1:3})) -gt 4 ]; then
    STIME=$(date -u --date="${LASTFTIME} +7 day" '+%Y%m%d')
  else
    STIME=$(date -u --date="${LASTFTIME} +1 day" '+%Y%m%d') 
    crest=$(checkwic ${STIME} ${ENSM})
    if [ "$crest" -eq 0 ]; then
      echo " OK - Restart files. "
    else
      echo " error - Restart files missing or wrong: "${STIME}" . Back to cycle "${LASTFTIME}" ."
      # re-running previous cycle, something happened and the restart was not generated
      for grid in ${LMDEFS[@]:1:3}; do

        test -f ${DIRRS}"/stream"${STRM}"/"${ENSM}"/"${STIME}".030000.restart."${grid}
        TE=$?
        if [ "$TE" -eq "0" ]; then
          rm -f ${DIRRS}/stream${STRM}/${ENSM}/${STIME}.030000.restart.${grid}
        fi

        test -f ${DIROUT}"/stream"${STRM}"/"${ENSM}"/out_grd."${LASTFTIME}"."${grid}
        TE=$?
        if [ "$TE" -eq "0" ]; then
          rm -f ${DIROUT}/stream${STRM}/${ENSM}/out_grd.${LASTFTIME}.${grid}
        fi

      done

      test -f ${DIROUT}"/stream"${STRM}"/"${ENSM}"/out_pnt."${LASTFTIME}"."${LMDEFS[4]}
      TE=$?
      if [ "$TE" -eq "0" ]; then
        rm -f ${DIROUT}/stream${STRM}/${ENSM}/out_pnt.${LASTFTIME}.${LMDEFS[4]}
      fi

      rm -f ${DIRJOUT}"/jww3gefs_"${STIME}"_stream"${STRM}"."${ENSM}".out"
      STIME="${LASTFTIME}"
      rm -f ${DIRJOUT}"/jww3gefs_"${STIME}"_stream"${STRM}"."${ENSM}".out"

    fi 
  fi
  
  # next cycle
  export STIME="$STIME"
  # submit new job
  rm -f ${DIRJOUT}"/jww3gefs_"${STIME}"_stream"${STRM}"."${ENSM}".out"
  sbatch --output=${DIRJOUT}"/jww3gefs_"${STIME}"_stream"${STRM}"."${ENSM}".out" ${DIRSCRIPTS}"/stream"${STRM}"/jwgefs0"${STRM}".sh"
  echo " job "${STIME}"_stream"${STRM}"."${ENSM}" submitted OK at "$(date +"%T")
  sleep 600
  wait $!

  # Wait while job is waiting in queue
  checkwaiting=0
  while [ "${checkwaiting}" -eq "0" ]; do
    test -f ${DIRJOUT}"/jww3gefs_"${STIME}"_stream"${STRM}"."${ENSM}".out"
    TE=$?
    if [ "$TE" -eq "0" ]; then
      test -f ${DIRW}/runningnow.txt
      TE=$?
      if [ "$TE" -eq "0" ]; then  
        # job has started and .out exist
        checkwaiting=1
        echo " job "${STIME}"_stream"${STRM}"."${ENSM}" started at "$(date +"%T")
        sleep 15000
        wait $!
      else
        # job waiting in queue, wait 5 min more
        sleep 300
        wait $!
      fi
    else
      # job waiting in queue, wait 5 min more
      sleep 300
      wait $!
    fi
  done

  # on hold while ww3 is still running
  checkrunning=0
  while [ "${checkrunning}" -eq "0" ]; do
    test -f ${DIRJOUT}"/jww3gefs_"${STIME}"_stream"${STRM}"."${ENSM}".out"
    TE=$?
    if [ "$TE" -eq "0" ]; then
      test -f ${DIRW}/runningnow.txt
      TE=$?
      if [ "$TE" -eq "1" ]; then  
        # ww3 not running
        checkrunning=1
        echo " job "${STIME}"_stream"${STRM}"."${ENSM}" finished at "$(date +"%T")
      else
        # job waiting in queue, wait 5 min more
        sleep 300
        wait $!
      fi
    else
      # job waiting in queue, wait 5 min more
      sleep 300
      wait $!
    fi
  done

  echo "  "

done

echo " Complete at "$(date +"%T")

