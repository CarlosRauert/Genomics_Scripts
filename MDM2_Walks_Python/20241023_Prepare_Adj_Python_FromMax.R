library(gGnome)
library(parallel)
CaseID_L <- list("CPCT02010386T",
                "CPCT02010680T",
                "CPCT02060104T",
                "CPCT02060191T",
                "CPCT02070051T",
                "CPCT02070366T",
                "CPCT02080206T",
                "CPCT02080227T",
                "CPCT02090057T",
                "CPCT02340046T",
                "GAYA01080001T",
                "WIDE01010374T",
                "WIDE01010421T",
                "WIDE01010452T",
                "WIDE01010517T",
                "WIDE01010527T",
                "WIDE01010736T",
                "WIDE01010743T",
                "WIDE01010781T",
                "WIDE01010822T",
                "WIDE01010851T",
                "WIDE01010861T",
                "WIDE01010917T",
                "WIDE01010941T",
                "WIDE01010958T",
                "WIDE01011015T",
                "A3M1")

lapply(CaseID_L, function(CaseID){
  Ggraph <-  gG(jabba=paste0("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/Jabba_Calls_TCGA/",CaseID,"/jabba.simple.rds"))
  CNxt <- max(na.omit(Ggraph$dt$cn))
  print(paste(CaseID,CNxt))
})

lapply(CaseID_L, function(CaseID){
  print(CaseID)
  dir.create(paste0("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/",CaseID))
  Ggraph <-  gG(jabba=paste0("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/Jabba_Calls_TCGA/",CaseID,"/jabba.simple.rds"))
  CNxt <- max(na.omit(Ggraph$dt$cn))
  highcopyX = Ggraph        # subset to high copy only
  A = highcopyX$adj
  colnames(A) = rownames(A) = highcopyX$gr$snode.id
  A = as.matrix(A)
  write.csv(A, paste0("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/",CaseID,"/",CaseID,"_adj_Matrix_Max",CNxt,".csv"))
})