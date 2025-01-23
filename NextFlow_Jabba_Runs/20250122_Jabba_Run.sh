# Try test profile
cd NF_Jabba
nextflow run . -profile test 

# Does not work due to incompatible sample sheet. Ergo get some TCGA fastqs and make proper sample sheet.
cd /data/cephfs-1/home/users/rauertc_c/scratch/NF_Jabba/Bams/AB36
export GDC_TOKEN=/data/cephfs-1/home/users/rauertc_c/work/GDC_Datatransfer/gdc-user-token.2025-01-22T15_32_13.891Z.txt
gdc-client download -m /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/NextFlow_Jabba_Runs/20250122_Test_Manifest_AB36.txt -t /data/cephfs-1/home/users/rauertc_c/work/GDC_Datatransfer/gdc-user-token.2025-01-22T15_32_13.891Z.txt

find "$(pwd)" -type f
