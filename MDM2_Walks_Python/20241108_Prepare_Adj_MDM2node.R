library(gGnome)

Args <- commandArgs(trailingOnly=TRUE)
CaseID <- Args[1]
CNxt <- as.integer(Args[2])

CaseID <- "A1KU"
CNxt <- 38

print(paste0("Read ",CaseID," and ",CNxt ," into R"))

prepare_adj <- function(CaseID, CNxt){
  Ggraph <-  gG(jabba=paste0("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/Jabba_Calls_TCGA/",CaseID,"/jabba.simple.rds"))
  highcopyX = Ggraph[cn>CNxt]        # subset to high copy only
  A = highcopyX$adj
  colnames(A) = rownames(A) = highcopyX$gr$snode.id
  A = as.matrix(A)
  write.csv(A, paste0("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/",CaseID,"/",CaseID,"_adj_Matrix_",CNxt,".csv"))
  nodes_list <- highcopyX$gr
  # Define MDM2 genomic coordinates
  mdm2_region <- GRanges(seqnames = "chr12", ranges = IRanges(start = 69202258, end = 69233629))
  seqlevelsStyle(mdm2_region) <- seqlevelsStyle(nodes_list)
  # Find overlaps between MDM2 and the GRangesList (nodes_list)
  mdm2_overlaps <- findOverlaps(nodes_list, mdm2_region, type = "any")
  mdm2_plusStrand <- mdm2_overlaps@from[1:(length(mdm2_overlaps@from)/2)]
  mdm2_minusStrand <- setdiff(mdm2_overlaps@from, mdm2_plusStrand)
  mdm2nodes <- rbind(mdm2_plusStrand, mdm2_minusStrand)
  write.csv(A, paste0("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/",CaseID,"/",CaseID,"_MDM2_nodes",CNxt,".csv"))
}

prepare_adj(CaseID, CNxt)