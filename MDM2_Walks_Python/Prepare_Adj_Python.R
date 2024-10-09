library(gGnome)

Args <- commandArgs(trailingOnly=TRUE)
CaseID <- Args[1]
CNxt <- as.integer(Args[2])

print(paste0("Read ",CaseID," and ",CNxt ," into R"))

prepare_adj <- function(CaseID, CNxt){
  Ggraph <-  gG(jabba=paste0("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/Jabba_Calls_TCGA/",CaseID,"/jabba.simple.rds"))
  highcopyX = Ggraph[cn>CNxt]        # subset to high copy only
  A = highcopyX$adj
  colnames(A) = rownames(A) = highcopyX$gr$snode.id
  A = as.matrix(A)
  write.csv(A, paste0("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/",CaseID,"/",CaseID,"_adj_Matrix_",CNxt,".csv"))
}

prepare_adj(CaseID, CNxt)