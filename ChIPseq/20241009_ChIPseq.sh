conda activate NF_CS
cd /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq

export NXF_APPTAINER_CACHEDIR=/data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq
export NXF_SINGULARITY_CACHEDIR=/data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq

nextflow run nf-core/chipseq -profile test --outdir /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq -c /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/ChIPseq/20241014_edit_bih.config

