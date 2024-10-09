library(gGnome)

Args <- commandArgs(trailingOnly=TRUE)
CaseID <- Args[1]
CNxt <- as.integer(Args[2])

print(paste0("Read ",CaseID," and ",CNxt ," into R"))

processFile = function(filepath) {
  LineList <- list()
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    LineList <- append(LineList, line)
  }
  
  close(con)
  return(LineList)
}

PythonLists2RLists <- function(StrList){
  NewList <- list()
  for (i in 1:length(StrList)){
    LineString <- StrList[i]
    LineString <- substring(LineString, 3, nchar(LineString)-2)
    LineList <- as.list(strsplit(LineString, "', '"))
    NewList <- append(NewList, LineList)
  }
  for (x in 1:length(NewList)){
    NewList[x] <- list(as.integer(unlist(NewList[x])))
  }
  return(NewList)
}

Python_Output <- processFile(paste0("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/",CaseID,"/",CaseID,"_",CNxt,"_Cycles.txt"))
Cycles_List <- PythonLists2RLists(Python_Output)
Increment <- length(Cycles_List)*1
Increment_Plot <- length(Cycles_List)*5

gencode <- track.gencode(cached.dir="/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics", stack.gap = 1e5, cex.label = 0.8, height = 20)
gencode$legend=FALSE
Ggraph <-  gG(jabba=paste0("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/Jabba_Calls_TCGA/",CaseID,"/jabba.simple.rds"))
highcopyX = Ggraph[cn>CNxt]        # subset to high copy only                                                                                                                                                   

Cycles_gW <- gW(snode.id=Cycles_List, graph=highcopyX, circular = TRUE)
saveRDS(Cycles_gW, paste0("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/",CaseID,"/",CaseID,"_",CNxt,"_Cycles.rds"))
pdf(file = paste0('/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/',CaseID,'/CycleWalks_CNmin_', CNxt,'.pdf'), width=9, height=4)
plot(c(gencode, highcopyX$gt, Cycles_gW$gtrack(name = "circular walks")), Cycles_gW$footprint+4e5, cex.label=0.3)
title(main=paste0("MDM2 circular walks with CN greater than ",CNxt))
par(cex.lab=0.1)
dev.off()  




