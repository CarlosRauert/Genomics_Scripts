for file in /data/cephfs-1/home/users/rauertc_c/scratch/ChIPseq/20241009_Acc_Lists*; do
    prefetch --option-file ${file} --output-directory /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq/SRA
done