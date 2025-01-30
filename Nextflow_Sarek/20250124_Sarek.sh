nextflow run nf-core/sarek \
   -profile test,singularity \
   --outdir ./out 

nextflow run nf-core/sarek \
   -r 3.5.0 \
   -profile singularity \
   