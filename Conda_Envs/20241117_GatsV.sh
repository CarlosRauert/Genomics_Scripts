conda create -n Gatsv r-base
conda install libxml2
BiocManager::install(c("BiocGenerics","GenomeInfoDb", "GenomicRanges", "IRanges", "S4Vectors"))
install.packages(c("caTools", "e1071", "parallel", "data.table", "rlang", "ROCR", "stats4", "stringr", "here", "optparse"))