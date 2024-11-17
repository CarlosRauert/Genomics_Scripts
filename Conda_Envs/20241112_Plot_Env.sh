conda ycreate -n Plot_Env r-base=4.0.2
install.packages("devtools")
export PKG_CONFIG_PATH=/data/cephfs-1/home/users/rauertc_c/work/miniforge3/envs/Plot_Env/lib/pkgconfig:$PKG_CONFIG_PATH
devtools::install_github("mskilab/gGnome")
# Skipping 7 packages not available: [Biostrings, S4Vectors, BiocGenerics, GenomeInfoDb, IRanges, rtracklayer, GenomicRanges]
packageurl <- "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.2-18.tar.gz"
install.packages(packageurl, repos = NULL, type = "source")
devtools::install_github("mskilab/gGnome", dependencies = TRUE)
install.packages("https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-53.tar.gz", repos = NULL, type = "source")
install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install(c("Biostrings", "S4Vectors", "BiocGenerics", 
                       "GenomeInfoDb", "IRanges", "rtracklayer", 
                       "GenomicRanges"))
devtools::install_github("mskilab/gTrack", dependencies = TRUE)
BiocManager::install("BSgenome")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")


