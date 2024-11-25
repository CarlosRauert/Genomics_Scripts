# get all functions
Rscript /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/20241125_Annotate_Gridss_Functions.R
Rscript /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/annotation_scripts.R

#set variables
DIR="/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs"
OUTPUT_DIR="/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/GetScaled"
SCRIPT="/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/20241125_GatsV_Run.R"

# Loop over Sample IDs
for item in $(ls "$DIR" | sed 's/\..*//' | sort -u); do
    # Set variable with .purple.sv.vcf.gz extension
    variable="${item}.purple.sv.vcf.gz"
    # Run the R script with the required arguments
    Rscript "$SCRIPT" "$variable" "$item" "$OUTPUT_DIR"
done