ssh -l duboisf_c hpc-login-2.cubi.bihealth.org
tmux new -s nfnWGS_BT474
tmux a -t nfnWGS_BT474

mkdir -p /data/cephfs-1/scratch/groups/dubois/users/duboisf_c/data/20241010nanoporeCL/epi2merunBT474
mkdir -p /data/cephfs-2/unmirrored/groups/dubois/me_proj/data/20241010nanoporeCL/epi2meBT474
srun --time 5-00 --mem=4G --ntasks=4 --pty bash -i
export NXF_APPTAINER_CACHEDIR=/data/cephfs-1/work/projects/dubois-lpwgs/singularity
export NXF_SINGULARITY_CACHEDIR=/data/cephfs-1/work/projects/dubois-lpwgs/singularity
export TOWER_ACCESS_TOKEN=eyJ0aWQiOiA5MzMxfS5kYjEzYmU3NDE3NjNhZjViMzBjODMwODE2NTUyZjczNjUxZGJkYTQy
cd /data/cephfs-1/scratch/groups/dubois/users/duboisf_c/data/20241010nanoporeCL/epi2merunBT474
conda activate env_nf
# check if ref is present
ls /data/cephfs-1/scratch/groups/dubois/me_proj/refs/
cd  /data/cephfs-1/scratch/groups/dubois/me_proj/refs/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
cd /data/cephfs-1/scratch/groups/dubois/users/duboisf_c/data/20241010nanoporeCL/epi2merunBT474

# v2.4.0

nextflow run /data/cephfs-1/work/groups/dubois/users/duboisf_c/nanopore/pipes/wf-human-variation --sample_name 'BT474' \
	--cnv --use_qdnaseq --sv --mod --phased --sample BT474 \
    --bam /data/cephfs-2/unmirrored/groups/dubois/me_proj/data/20241005bih_lrWGS/20240917_A4302_Lib252/20240917_A4302_Lib252_Pool_NB17-21/20240917_1201_2F_PAW94493_d3b9f6cb/bam_pass/barcode17  \
    --ref /data/cephfs-1/scratch/groups/dubois/me_proj/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz --bam_min_coverage 4 \
    --out_dir /data/cephfs-2/unmirrored/groups/dubois/me_proj/data/20241010nanoporeCL/epi2meBT474 --threads 64 --modkit_threads 32 \
    -c /data/cephfs-1/work/groups/dubois/users/duboisf_c/lpWGS/bih_hpc_min_20240906.config --override_basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0' -with-tower \
    -profile singularity -resume
