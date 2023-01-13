#!/bin/bash --login
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=40
#SBATCH -q batch
#SBATCH -t 07:50:00
#SBATCH -A marine-cpu
#SBATCH -p orion
#SBATCH --exclusive

# Main job to run WW3 GEFSv12 reforecast, which is run inside (executed by)
#  the manager script jobmanager.sh.
# 
# DIRJOUT="/work/noaa/marine/ricardo.campos/work/ww3runs/jobs"
# export STIME="20000101"
# export STRM="1"
# export ENSM="c00"
# sbatch --output=${DIRJOUT}/jww3gefs_${STIME}_stream${STRM}.${ENSM}.out jwgefs.sh

ulimit -s unlimited
ulimit -c 0

# input variables and functions
source /work/noaa/marine/ricardo.campos/work/ww3runs/scripts/input_varsfuncts.sh 
# work directory
DIRW="${DIRWORK}/stream${STRM}/${ENSM}"
# flag file to indicate model is running and job is active
echo " "${STIME}" "${STRM}" "${ENSM}" " > ${DIRW}/runningnow.txt
# output log file
echo " Starting at "$(date +"%T") > ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
echo " " >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt

cd "${DIRW}"
for iday in $(seq $NDAYS); do

  NDATE=$(date -u --date="${STIME} +$((${iday}-1)) day" '+%Y%m%d')

  # flag file to indicate model is running and job is active
  echo " "${NDATE}" "${STRM}" "${ENSM}" " >> ${DIRW}/runningnow.txt

  echo "${NDATE} ${STRM} ${ENSM} " >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
  echo " " >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
  echo " "
  echo " ------------------------------------------------------------- "
  echo " CYCLE TIME ${NDATE}"
  echo " ------------------------------------------------------------- "
  echo " "

  cout=$(checkwout ${NDATE} ${STRM} ${ENSM})
  if [ "$cout" -eq 0 ] || [ $((10#${ENSM:1:3})) -gt 4 -a $(date -d "$NDATE" +%u) -ne 3 ]; then

    if [ "$cout" -eq 0 ]; then
      echo " cycle ${NDATE} ${STRM} ${ENSM} already run. Ok. Moving to the next one ..." >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
      echo " " >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt 
    fi

    if [ $((10#${ENSM:1:3})) -gt 4 -a $(date -d "$NDATE" +%u) -ne 3 ]; then
      echo " Not Wednesday ${NDATE}. Extra ensemble member ${ENSM} only available on Wednesdays. Moving to the next one ..." >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
      echo " " >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt 
    fi

  else

    if [  $((10#${ENSM:1:3})) -gt 4 ]; then

      restartid=$(pickrestartid ${NDATE} ${STRM} ${ENSM})
      wait $!
      sleep 10
      if [ "$restartid" -eq "00" ]; then
        RENSM="c"${restartid}
      else
        RENSM="p"${restartid}
      fi

      LMDEFS=($MDEFS)
      # symbolic links restart files
      for grid in ${LMDEFS[@]:1:3}; do
        test -f ${DIRRS}"/stream"${STRM}"/"${RENSM}"/"${NDATE}".030000.restart."${grid}
        TE=$?
        if [ "$TE" -eq "0" ]; then
          ln -fs ${DIRRS}"/stream"${STRM}"/"${RENSM}"/"${NDATE}".030000.restart."${grid} ${DIRRS}"/stream"${STRM}"/"${ENSM}"/"${NDATE}".030000.restart."${grid}
        else
          echo " ERROR - Missing restart "${DIRRS}"/stream"${STRM}"/"${RENSM}"/"${NDATE}".030000.restart."${grid} >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
          rm -f ${DIRW}/runningnow.txt
          exit 0
        fi
      done

    fi

    # Prep module ----
    cprep=1; TRIES=1
    while [ ${cprep} -eq 1 ] && [ ${TRIES} -le 3 ]; do
      cprep=$(checkwprep ${NDATE} ${ENSM})
      if [ "$cprep" -eq 1 ]; then
        echo " error - Input files missing: wind and/or ice. Running prep_windice.sh ..." >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
        sh ${DIRSCRIPTS}/stream${STRM}/prep_windice.sh ${NDATE} ${ENSM} &&
        wait $!
        sleep 5
      fi
      TRIES=`expr $TRIES + 1`
    done

    crest=$(checkwic ${NDATE} ${ENSM})
    if [ "$crest" -eq 0 ]; then
      echo " OK - Restart files. " >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
    else
      echo " ERROR - Restart files missing or wrong: "${NDATE}" ." >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
      rm -f ${DIRW}/runningnow.txt
      for grid in ${LMDEFS[@]:1:3}; do
        test -f ${DIRRS}"/stream"${STRM}"/"${ENSM}"/"${NDATE}".030000.restart."${grid}
        TE=$?
        if [ "$TE" -eq "0" ]; then
          rm -f ${DIRRS}/stream${STRM}/${ENSM}/${NDATE}.030000.restart.${grid}
        fi
      done
      exit 0
    fi
    # ----

    cprep=$(checkwprep ${NDATE} ${ENSM})
    if [ "$cprep" -eq 0 ]; then
      echo " OK - Input and IC files: restart, wind, and ice. Running prep_multi ..." >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
      sh ${DIRSCRIPTS}/stream${STRM}/prep_multi.sh ${NDATE} ${STRM} ${ENSM} &&
      wait $!
      sleep 10
    else
      echo " ERROR - Input files missing, "${NDATE}" .">> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
      rm -f ${DIRW}/runningnow.txt
      exit 0
    fi
    # ----

    echo "   running ww3_multi ... "${TRIES} >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
    srun -n 400 ${WW3MAIN}/exe/ww3_multi
    wait $!
    sleep 40
    # organize outputs
    echo "   organizing outputs ... "${TRIES} >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
    sh ${DIRSCRIPTS}/stream${STRM}/org_out.sh ${NDATE} ${STRM} ${ENSM} &&
    wait $!
    sleep 20

    # Final check the entire cycle went well.
    cout=$(checkwout ${NDATE} ${STRM} ${ENSM})
    if [ "$cout" -eq 0 ]; then
      echo " OK - Expected outputs. Cycle complete at "$(date +"%T") >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
    else
      echo " ERROR - Missing or wrong expected outputs." >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt       
      rm -f ${DIRW}/runningnow.txt
      for grid in ${LMDEFS[@]:1:3}; do
        test -f ${DIROUT}"/stream"${STRM}"/"${ENSM}"/out_grd."${NDATE}"."${grid}
        TE=$?
        if [ "$TE" -eq "0" ]; then
          rm -f ${DIROUT}/stream${STRM}/${ENSM}/out_grd.${NDATE}.${grid}
        fi
      done
      test -f ${DIROUT}"/stream"${STRM}"/"${ENSM}"/out_pnt."${NDATE}"."${LMDEFS[4]}
      TE=$?
      if [ "$TE" -eq "0" ]; then
        rm -f ${DIROUT}/stream${STRM}/${ENSM}/out_pnt.${NDATE}.${LMDEFS[4]}
      fi
      exit 0
    fi

    echo "  " >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt  
    echo " "
    echo " ------------------------------------------------------------- "
    echo " DONE"
    echo " ------------------------------------------------------------- "
    echo " "

  fi

done

echo " Complete at "$(date +"%T") >> ${DIRJOUT}/log_${STIME}_stream${STRM}_${ENSM}.txt
sleep 10
rm -f ${DIRW}/runningnow.txt
wait $!
echo " End job."

