#!/bin/bash

start_date="20201001"
end_date="20250331"
min_size=10000000  # in bytes
output_file="missing_or_small_files.txt" > "$output_file"

# Date loop
current=$(date -d "$start_date" +%s)
end=$(date -d "$end_date" +%s)

while [ "$current" -le "$end" ]; do
    ymd=$(date -u -d "@$current" +%Y%m%d)
    dir="GEFSv12Waves_${ymd}"

    # Loop over ensemble members
    for em in $(seq -w 0 30); do
        # Loop over forecast hours
        for fh in $(seq -w 0 6 384); do
            file="${dir}/gefs.wave.${ymd}.${em}.global.0p25.f${fh}.grib2"

            if [ ! -f "$file" ]; then
                echo "$file" >> "$output_file"
            else
                size=$(du -sb "$file" | cut -f1)
                if [ "$size" -lt "$min_size" ]; then
                    echo "$file" >> "$output_file"
                fi
            fi
        done
    done

    current=$(( current + 86400 ))  # next day
done

echo "Done. Missing or small files listed in $output_file"

