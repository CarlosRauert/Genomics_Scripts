#!/bin/sh
source /data/cephfs-1/home/users/rauertc_c/work/miniforge3/etc/profile.d/conda.sh
conda init
conda activate Gatsv
#set variables
DIR="/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs"
OUTPUT_DIR="/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/GetScaled/"
SCRIPT="/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/20241125_GatsV_Run.R"

# Loop over Sample IDs
#for item in $(ls "$DIR" | sed 's/\..*//' | sort -u); do
#    # Set variable with .purple.sv.vcf.gz extension
#    variable="${DIR}/${item}.purple.sv.vcf.gz"
#    # Run the R script with the required arguments
#    Rscript "$SCRIPT" "$variable" "$item" "$OUTPUT_DIR"
#done

# Loop over Sample IDs for a specified range (e.g., 1001:2000)
start=$1
end=$2

# Extract the range of items based on the specified start and end
count=0
for item in $(ls "$DIR" | sed 's/\..*//' | sort -u | sed -n "${start},${end}p"); do
    # Increment the counter
    count=$((count + 1))

    # Set variable with .purple.sv.vcf.gz extension
    variable="${DIR}/${item}.purple.sv.vcf.gz"
    echo "Processing Count: $count, Sample ID: $item"
    # Run the R script with the required arguments
    Rscript "$SCRIPT" "$variable" "$item" "$OUTPUT_DIR"

done