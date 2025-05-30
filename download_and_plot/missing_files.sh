#!/bin/bash

# Output file
missing_file_list="missing_files.txt" > "$missing_file_list"  # Clear file if exists

# Date range
start_date="20201001"
end_date="20250331"

current=$(date -d "$start_date" +%s)
end=$(date -d "$end_date" +%s)

# Loop over each day
while [ "$current" -le "$end" ]; do
    ymd=$(date -u -d "@$current" +%Y%m%d)
    
    # Loop over each ensemble member (00 to 30)
    for em in $(seq -w 0 30); do
        filename="gefs.wave.${ymd}.${em}.global.0p25.nc"
        
        if [ ! -f "$filename" ]; then
            echo "$filename" >> "$missing_file_list"
        fi
    done

    # Advance to next day
    current=$(( current + 86400 ))
done

echo "Done. Missing files are listed in $missing_file_list"
