conda activate NF_CS
cd /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq

nextflow run nf-core/chipseq -profile test,bih --outdir /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq --scratch "/data/cephfs-1/scratch/groups/dubois/users/rauertc_c/"
