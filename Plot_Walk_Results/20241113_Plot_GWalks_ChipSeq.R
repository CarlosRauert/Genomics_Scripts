library(gGnome)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)  # Example for human hg19 genome
library(BSgenome.Hsapiens.1000genomes.hs37d5)

seqinfo_hg19 <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
wigToBigWig("/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/Adipocyte_Nuclei_GSM916066/GSM916066_BI.Adipose_Nuclei.H3K27ac.7.wig.gz", seqinfo_hg19)

file1 <- "/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/Adipocyte_Nuclei_GSM916066/GSM916066_BI.Adipose_Nuclei.H3K27ac.7.bw"
file2 <- "/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/Adipose_Tissue_Encode/ENCFF083HLK.bigWig"
ChipSeq_Nucleus <- gTrack(file1, name= "Adipocyte Nuclei \n raw signal density", cex.label =0.3, height=20)
ChipSeq_Tissue <- gTrack(file2, name=" Adipose Tissue \n Fold enrichment", height=20)
Lipo863B_REP1 <- gTrack("/data/cephfs-1/home/users/rauertc_c/work/ChipSeq/Lipo863B_mod_REP1.mLb.clN.bigWig",
                        name="Lipo836B", chr.sub=FALSE)
LP6_REP1 <- gTrack("/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/LP6_mod_REP1.mLb.clN.bigWig",
                        name="LP6", chr.sub=FALSE)
LPS141 <- gTrack("/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/LPS141_mod_REP1.mLb.clN.bigWig",
                        name="LPS141", chr.sub=FALSE)
SW872_REP1 <- gTrack("/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/SW872_mod_REP1.mLb.clN.bigWig",
                        name="SW872_REP1", chr.sub=FALSE)        
SW872_REP2 <- gTrack("/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/SW872_mod_REP2.mLb.clN.bigWig",
                        name="SW872_REP2", chr.sub=FALSE)
Tumor_1 <- gTrack("/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/Tumor_1_mod_REP1.mLb.clN.bigWig",
                        name="Tumor 1", chr.sub=FALSE)
Tumor_2 <- gTrack("/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/Tumor_2_mod_REP1.mLb.clN.bigWig",
                        name="Tumor 2", chr.sub=FALSE)

TissueEncode <- import()
AdiposeNuclei <- import(file1)
L863BR1 <- import("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241014_NF_CS_Out/bwa/merged_library/bigwig/Lipo863B_REP1.mLb.clN.bigWig")
seqlevels(L863BR1) <- paste0("chr", seqlevels(L863BR1))
export(L863BR1, "/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/Lipo863B_mod_REP1.mLb.clN.bigWig")
LP6_REP1 <- import("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241014_NF_CS_Out/bwa/merged_library/bigwig/LP6_REP1.mLb.clN.bigWig")
seqlevels(LP6_REP1) <- paste0("chr", seqlevels(LP6_REP1))
export(LP6_REP1, "/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/LP6_mod_REP1.mLb.clN.bigWig")
LPS141 <- import("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241014_NF_CS_Out/bwa/merged_library/bigwig/LPS141_REP1.mLb.clN.bigWig")
seqlevels(LPS141) <- paste0("chr", seqlevels(LPS141))
export(LPS141, "/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/LPS141_mod_REP1.mLb.clN.bigWig")
SW872R1 <- import("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241014_NF_CS_Out/bwa/merged_library/bigwig/SW872_REP1.mLb.clN.bigWig")
seqlevels(SW872R1) <- paste0("chr", seqlevels(SW872R1))
export(SW872R1, "/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/SW872_mod_REP1.mLb.clN.bigWig")
SW872R2 <- import("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241014_NF_CS_Out/bwa/merged_library/bigwig/SW872_REP2.mLb.clN.bigWig")
seqlevels(SW872R2) <- paste0("chr", seqlevels(SW872R2))
export(SW872R2, "/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/SW872_mod_REP2.mLb.clN.bigWig")
Tumor_1 <- import("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241014_NF_CS_Out/bwa/merged_library/bigwig/Tumor_1_REP1.mLb.clN.bigWig")
seqlevels(Tumor_1) <- paste0("chr", seqlevels(Tumor_1))
export(Tumor_1, "/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/Tumor_1_mod_REP1.mLb.clN.bigWig")
Tumor_2 <- import("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241014_NF_CS_Out/bwa/merged_library/bigwig/Tumor_2_REP1.mLb.clN.bigWig")
seqlevels(Tumor_2) <- paste0("chr", seqlevels(Tumor_2))
export(Tumor_2, "/data/cephfs-1/home/users/rauertc_c/work/genomics/ChIPseq/Tumor_2_mod_REP1.mLb.clN.bigWig")




plotGwalk_ChipSeq <- function(gWalk, Name){
  Gwalk_Obj <- readRDS(gWalk)
  walks_l = Gwalk_Obj[walk.id %in% names(sort(Gwalk_Obj$lengths, decreasing = T)[1:3])]
  gencode <- track.gencode(cached.dir="/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics", stack.gap = 1e5, cex.label = 0.8, height = 20)
  gencode$legend=FALSE
  pdf(file = paste0('/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/',Name,'/20241311_CycleWalks_ChipSeq_Test.pdf'), height=15, width=9)                                                                                       
  #plot(c(gencode, walks_l$gtrack(name = "circular walks"), ChipSeq_Tissue, ChipSeq_Nucleus), walks_l$footprint+2e5, cex.label=0.4, cex.axis=0.2)
  plot(c(gencode, walks_l$gt, LP6_REP1, LPS141, SW872_REP1, Tumor_1, Lipo863B_REP1), walks_l$footprint+2e5, cex.label=0.4, cex.axis=0.2)
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
