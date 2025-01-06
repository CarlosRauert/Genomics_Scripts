ssh rauertc_c@hpc-login-1.cubi.bihealth.org -p 22
srun --mem=8G --ntasks=12 --pty bash -i
conda create --name Genomics
conda activate Genomics
conda install r-devtools jupyter networkx
R
devtools::install_github("mskilab/gGnome")
#* DONE (gGnome)
#Warning messages:
#1: In i.p(...) :
#  installation of package ‘/tmp/RtmpA3o5W9/file59a9e611ac12b/VariantAnnotation_1.51.1.tar.gz’ had non-zero exit status
#2: packages ‘S4Vectors’, ‘BiocGenerics’, ‘GenomeInfoDb’, ‘IRanges’, ‘rtracklayer’, ‘GenomicRanges’ are not available for this version of R
conda install bioconductor-s4vectors bioconductor-genomicranges bioconductor-iranges bioconductor-biocgenerics bioconductor-genomeinfodb bioconductor-rtracklayer -c conda-forge
tmux new-session -A -s Genomics2
BiocManager::install(c("S4Vectors","GenomicRanges","IRanges","BiocGenerics","GenomeInfoDb","rtracklayer"),force=TRUE)