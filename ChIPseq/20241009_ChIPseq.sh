conda activate NF_CS
cd /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq

export NXF_APPTAINER_CACHEDIR=/data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq
export NXF_SINGULARITY_CACHEDIR=/data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq

nextflow run nf-core/chipseq -profile test,bih --outdir /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq

