library(gGnome)
library(rtracklayer)


file1 <- "/data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq/Adipocyte_Nuclei_GSM916066/Adipocyte.bw"
file2 <- "/data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq/Adipose_Tissue_Encode/ENCFF083HLK.bigWig"

import("/data/cephfs-1/scratch/groups/dubois/users/rauertc_c/ChIPseq/Adipocyte_Nuclei_GSM916066/GSM916066_BI.Adipose_Nuclei.H3K27ac.7.wig.gz", format=)

ChipSeq_Nucleus <- gTrack(file1, name= "Adipocyte Nuclei \n raw signal density", cex.label =0.3, height=20)
ChipSeq_Tissue <- gTrack(file2, name=" Adipose Tissue \n Fold enrichment", height=20)



plotGwalk_ChipSeq <- function(gWalk, Name){
  Gwalk_Obj <- readRDS(gWalk)
  walks_l = Gwalk_Obj[walk.id %in% names(sort(Gwalk_Obj$lengths, decreasing = T)[1:3])]
  gencode <- track.gencode(cached.dir="/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics", stack.gap = 1e5, cex.label = 0.8, height = 20)
  gencode$legend=FALSE
  pdf(file = paste0('/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/',Name,'/CycleWalks_ChipSeq.pdf'), height=15, width=9)                                                                                       
  plot(c(gencode, walks_l$gtrack(name = "circular walks"), ChipSeq_Tissue, ChipSeq_Nucleus), walks_l$footprint+2e5, cex.label=0.4, cex.axis=0.2)
  title(main=Name)
  par(cex.lab=0.1)
  dev.off()                                                                      
}

plotGwalk_ChipSeq("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/A1KU/A1KU_19_Cycles.rds","A1KU")
plotGwalk_ChipSeq("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/A1L0/A1L0_18_Cycles.rds","A1L0")
plotGwalk_ChipSeq("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/A1L3/walks_circ_CNmin_60.5.rds","A1L3")
plotGwalk_ChipSeq("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/A2J4/A2J4_17_Cycles.rds","A2J4")
plotGwalk_ChipSeq("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/A2QS/A2QS_48_Cycles.rds","A2QS")
plotGwalk_ChipSeq("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/A3LS/walks_circ_CNmin_37.25.rds","A3LS")
plotGwalk_ChipSeq("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/A3LT/A3LT_125_Cycles.rds","A3LT")

R <- readRDS("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1KU/walksCirc_CNmin19.rds")
Py <- readRDS("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1KU/A1KU_19_Cycles.rds")

pdf(file = paste0('/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/',"A1KU",'/CycleWalks_Compare.pdf'), height=15, width=9)                                                                                       
plot(c(gencode, R$gtrack(name="R"),Py$gtrack(name="Python")),R$footprint+2e5, cex.label=0.4, cex.axis=0.2)
par(cex.lab=0.1)
dev.off()

plot(ChipSeq_Nucleus)
