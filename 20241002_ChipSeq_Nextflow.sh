ssh rauertc_c@hpc-login-1.cubi.bihealth.org -p 22
tmux new-session -A -s NF_CS
srun --mem=8G --ntasks=12 --pty bash -i
conda activate NF_CS
cd work 
conda install nf-core
nextflow run nf-core/chipseq -profile test,bih --outdir /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq/Nextflow_Outs --scratch /data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ --apptainer.pullTimeout="40m"

#ERROR ~ Error executing process > 'NFCORE_CHIPSEQ:CHIPSEQ:PICARD_MERGESAMFILES (6)'#
#
#Caused by:
#  Failed to pull singularity image
#    command: apptainer pull  --name quay.io-biocontainers-picard-2.27.4--hdfd78af_0.img.pulling.1727882740443 docker://quay.io/biocontainers/picard:2.27.4--hdfd78af_0 > /dev/null       
#    status : 255
#    hint   : Try and increase apptainer.pullTimeout in the config (current is "20m")
#  #  message:
# #     FATAL:   Failed to create an image cache handle: failed initializing caching directory: couldn't create cache directory /data/cephfs-1/home/users/rauertc_c/.apptainer/cache: mkdir /data/cephfs-1/home/users/rauertc_c/.apptainer: file exists
#
#
#
# -- Check '.nextflow.log' file for details