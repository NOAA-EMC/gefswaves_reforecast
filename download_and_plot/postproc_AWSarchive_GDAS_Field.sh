#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# DIRJOUT="/work/noaa/marine/ricardo.campos/data/archiveOPruns/GDASwave/jobs"
# export YEAR=2022
# export MONTH=1
# sbatch --output=${DIRJOUT}/postproc_AWSarchive_GDAS_Field_${YEAR}${MONTH}.out postproc_AWSarchive_GDAS_Field.sh

echo " Starting at "$(date +"%T")

ulimit -s unlimited
ulimit -c 0

# 2 input arguments
export YEAR=${YEAR}
export MONTH=${MONTH}

# Orion module loads
module load wgrib2
module load nco
module load cdo

# limits for domain selection
lonmin=-180.
lonmax=0.

# Paths
OPATH="/work/noaa/marine/ricardo.campos/data/archiveOPruns/GDASwave"
# Data path
DIRW=${OPATH}"/"${YEAR}`printf %2.2d $MONTH`
# destination path
DIRS=${OPATH}"/netcdf"

# cutoff decimals to reduce file size
dp=3
# Forecast lead time (hours)
fleads="`seq -f "%03g" 0 1 5`"
# Cycle (hours)
cycles="`seq -f "%02g" 0 6 18`"

cd ${DIRW}

# Post-processing: convert to netcdf, compress, and reduce decimals resolution, to save disk space. ------------------
for DAY in `seq 1 31`; do
  CTIME=${YEAR}`printf %2.2d $MONTH``printf %2.2d $DAY`
  for CHOUR in $cycles;do
    for h in $fleads;do

      arqn=$DIRW/gdaswave.${CTIME}.t"$(printf "%02.f" $CHOUR)"z.global.0p16.f"$(printf "%03.f" $h)"
      test -f ${arqn}.grib2
      TE=$?
      if [ "$TE" -eq 1 ]; then
         echo " File ${arqn}.grib2 does not exist. Failed to download " 
      else
         wgrib2 ${arqn}.grib2 -netcdf ${arqn}.saux.nc  2>&1
         wait $!
         cdo select,name='WIND_surface','HTSGW_surface','PERPW_surface' ${arqn}.saux.nc ${arqn}.saux1.nc
         wait $!
         cdo sellonlatbox,-180,180,-90,90 ${arqn}.saux1.nc ${arqn}.saux2.nc
         wait $!
         ncks -4 -L 1 -d longitude,${lonmin},${lonmax} ${arqn}.saux2.nc ${arqn}.saux3.nc
         wait $!
         ncks --ppc default=.$dp ${arqn}.saux3.nc ${arqn}.nc  2>&1
         wait $!
         ncatted -a _FillValue,,o,f,NaN ${arqn}.nc  2>&1
         wait $!
         # rm -f ${arqn}.grib2
         rm ${arqn}.saux*
         echo " File ${arqn} converted to netcdf and compressed with success. "
      fi

    done
  done

  ncrcat $DIRW/gdaswave.${CTIME}.t*z.global.0p16.f*.nc -O ${DIRS}/gdaswave.${CTIME}.global.0p16.nc
  echo " Done ${CTIME}"
done

echo "Completed postproc_AWSarchive_GDAS_Field.sh "${YEAR}`printf %2.2d $MONTH`"  at "$(date +"%T")

