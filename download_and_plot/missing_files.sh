#!/bin/bash

start_date="20201001"
end_date="20250331"

min_size=480000000

output_file="missing_or_small_nc_files.txt" > "$output_file"

# Date loop
current=$(date -d "$start_date" +%s)
end=$(date -d "$end_date" +%s)

while [ "$current" -le "$end" ]; do
    ymd=$(date -u -d "@$current" +%Y%m%d)

    for em in $(seq -w 0 30); do
        file="gefs.wave.${ymd}.${em}.global.0p25.nc"

        if [ ! -f "$file" ]; then
            echo "$file" >> "$output_file"
        else
            size=$(du -sb "$file" | cut -f1)
            if [ "$size" -lt "$min_size" ]; then
                echo "$file" >> "$output_file"
            fi
        fi
    done

    current=$(( current + 86400 ))
done

echo "Done. Missing or small .nc files are listed in $output_file"

