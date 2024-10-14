conda activate NF_CS
cd /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq

ln -s scratch/ChIPseq/cache_Files/.apptainer ~/.apptainer 
ln -s scratch/ChIPseq/cache_Files/.nextflow ~/.nextflow
ln -s scratch/ChIPseq/cache_Files/.local ~/.local
-c /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/ChIPseq/20241014_edit_bih.config

export NXF_APPTAINER_CACHEDIR=/data/cephfs-1/work/projects/dubois-lpwgs/singularity
export NXF_SINGULARITY_CACHEDIR=/data/cephfs-1/work/projects/dubois-lpwgs/singularity
# test
nextflow run nf-core/chipseq -profile test --outdir /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq -c /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/ChIPseq/20241014_edit_bih.config
# run with sample sheet
nextflow run nf-core/chipseq --genome GRCh37 \
    --input /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/ChIPseq/20241014_ChIPseq_Sample_sheet.csv \
    --outdir /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq/20241014_NF_CS_Out \
    -c /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/ChIPseq/20241014_edit_bih.config \
    --macs_gsize 2700000000
