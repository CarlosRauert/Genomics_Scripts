cd /data/cephfs-1/home/users/rauertc_c/scratch/NF_Jabba/Fastqs/AB36

nextflow run nf-core/bamtofastq \
-profile bih \
--input /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/NextFlow_Jabba_Runs/20250130_fastq_samplesheet_test.csv \
--outdir './results'