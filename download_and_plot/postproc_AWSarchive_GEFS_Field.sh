#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# DIRJOUT="/work/noaa/marine/ricardo.campos/data/archiveOPruns/GEFSv12Waves_AWS/jobs"
# export YEAR=2022
# export MONTH=1
# export DAY=1
# sbatch --output=${DIRJOUT}/postproc_AWSarchive_GEFS_Field_${YEAR}${MONTH}${DAY}.out postproc_AWSarchive_GEFS_Field.sh

echo " Starting at "$(date +"%T")

ulimit -s unlimited
ulimit -c 0

export YEAR=${YEAR}
export MONTH=${MONTH}
export DAY=${DAY}

module load nco
module load cdo
module load wgrib2

# date
CTIME=${YEAR}`printf %2.2d $MONTH``printf %2.2d $DAY`
DIRS="/work/noaa/marine/ricardo.campos/data/archiveOPruns/GEFSv12Waves_AWS" # archive path
DIRO="/work/noaa/marine/ricardo.campos/data/archiveOPruns/GEFSv12Waves_AWS/netcdf" # final netcdf4 output path

# cutoff decimals to reduce file size
dp=2

DIRW="${DIRS}/GEFSv12Waves_${CTIME}"
cd $DIRW

# hours of forecast lead time to be dowloaded
fleads="`seq -f "%03g" 0 6 384`"
# ensemble members
ensblm="`seq -f "%02g" 0 1 30`"

# Post-processing: compress, select variables, and reduce decimals resolution, to save disk space. ------------------
echo " "
echo " Post-Processing. Netcdf4 compression "${CTIME}

for h in $fleads;do
  for e in $ensblm;do
     arqn=$DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f"$(printf "%03.f" $h)"
     test -f ${arqn}.grib2
     TE=$?
     if [ ${TE} -eq 1 ]; then
       echo " File ${arqn}.grib2 does not exist. Download failed."
       echo " File ${arqn}.grib2 does not exist. Download failed." >> "${DIRS}/file_doesnt_extist.txt"
     else
       wgrib2 ${arqn}.grib2 -netcdf ${arqn}.saux.nc 2>&1
       wait $!
       cdo select,name='WIND_surface','HTSGW_surface','PERPW_surface' ${arqn}.saux.nc ${arqn}.saux2.nc
       wait $!
       ncks -4 -L 1 ${arqn}.saux2.nc ${arqn}.saux3.nc 2>&1
       wait $!
       ncks --ppc default=.$dp ${arqn}.saux3.nc ${arqn}.nc 2>&1
       wait $!
       ncatted -a _FillValue,,o,f,NaN ${arqn}.nc 2>&1
       wait $!
       rm -f ${arqn}.saux*
       rm -f ${arqn}.*idx*
       echo " File ${arqn} converted to netcdf and compressed with success. "
       sleep 1
     fi
  done
done

# Merge all netcdf files lead time
for e in $ensblm;do
  ncrcat $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f*.nc -O $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.nc 2>&1
  wait $!
  sleep 1
  mv $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.nc ${DIRO}
  rm -f $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f*.nc
  echo " All time steps of ensemble member ${e} have been merged. "
done

echo "Completed postproc_AWSarchive_GEFS_Field.sh "${CTIME}"  at "$(date +"%T")

