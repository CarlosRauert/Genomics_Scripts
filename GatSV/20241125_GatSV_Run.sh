sbatch --ntasks=32 --mem=30G /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/20241125_GatSV_RunStartEnd.sh 5001 5500
find /data/cephfs-1/work/groups/dubois/users/rauertc_c/GatSV/Out/GetScaled/ -type f -iname '*classified*' | wc -l


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
start=4001
end=4500

# Extract the range of items based on the specified start and end
count=0
for item in $(ls "$DIR" | sed 's/\..*//' | sort -u | sed -n "${start},${end}p"); do
    # Increment the counter
    count=$((count + 1))

    # Set variable with .purple.sv.vcf.gz extension
    variable="${DIR}/${item}.purple.sv.vcf.gz"

    # Run the R script with the required arguments
    Rscript "$SCRIPT" "$variable" "$item" "$OUTPUT_DIR" 

done