for file in /data/cephfs-1/home/users/rauertc_c/scratch/ChIPseq/20241009_Acc_Lists/*; do
    if [ -f "$file" ]; then
        echo $file
        prefetch --option-file ${file} --output-directory /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq/SRA 
    fi
done

cd /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq/FastQ

for dir in /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq/SRA/*
do
    echo $dir
    fasterq-dump $dir
done

# Loop through all .fastq files in the current directory
for file in /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq/FastQ/*.fastq
do
    # Compress each file with gzip
    gzip "$file"
done