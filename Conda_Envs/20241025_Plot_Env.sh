conda ycreate -n Plot_Env r-base=4.0.2
conda install r-devtools bioconductor-s4vectors bioconductor-genomicranges bioconductor-iranges bioconductor-biocgenerics bioconductor-genomeinfodb bioconductor-rtracklayer
devtools::install_github("mskilab/gGnome")
