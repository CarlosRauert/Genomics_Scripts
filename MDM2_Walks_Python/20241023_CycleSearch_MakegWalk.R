library(gGnome)
library(gTrack)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(readxl)
library(parallel)

CaseID_l = list("CPCT02010386T", "CPCT02010680T", "CPCT02060104T", "CPCT02060191T", "CPCT02070051T", "CPCT02070366T", "CPCT02080206T", "CPCT02080227T", "CPCT02090057T", "CPCT02340046T")
CNxt_l = list("91", "11", "36", "11", "9", "9", "9", "10", "18", "13", "10", "91")

Lipo863B_REP1 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241014_NF_CS_Out/bwa/merged_library/bigwig/Lipo863B_REP1.mLb.clN.bigWig",
                        name="Lipo836B", chr.sub=FALSE)
LP6_REP1 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241014_NF_CS_Out/bwa/merged_library/bigwig/LP6_REP1.mLb.clN.bigWig",
                        name="LP6", chr.sub=FALSE)
LPS141 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241014_NF_CS_Out/bwa/merged_library/bigwig/LPS141_REP1.mLb.clN.bigWig",
                        name="LPS141", chr.sub=FALSE)
SW872_REP1 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241014_NF_CS_Out/bwa/merged_library/bigwig/SW872_REP1.mLb.clN.bigWig",
                        name="SW872_REP1", chr.sub=FALSE)        
SW872_REP2 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241014_NF_CS_Out/bwa/merged_library/bigwig/SW872_REP2.mLb.clN.bigWig",
                        name="SW872_REP2", chr.sub=FALSE)
Tumor_1 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241014_NF_CS_Out/bwa/merged_library/bigwig/Tumor_1_REP1.mLb.clN.bigWig",
                        name="Tumor 1", chr.sub=FALSE)
Tumor_2 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241014_NF_CS_Out/bwa/merged_library/bigwig/Tumor_2_REP1.mLb.clN.bigWig",
                        name="Tumor 2", chr.sub=FALSE)

xGr=4

lapply(1:length(CaseID_l), function(xGr){
  CaseID <- CaseID_l[xGr]
  CNxt <- CNxt_l[xGr]
  Cycles_File=paste0("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/", CaseID,"/optimized_cycles_",CNxt,".txt")

  print(paste0("Read ",CaseID," and ",CNxt ," into R"))
  lines <- readLines(Cycles_File, n = 5)  # Adjust "Cycles_File.txt" to your file path
  if(length(lines)>0){
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
 
    Python_Output <- processFile(Cycles_File)
    Cycles_List <- PythonLists2RLists(Python_Output)

    Ggraph <-  gG(jabba=paste0("/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/Jabba_Calls_TCGA/",CaseID,"/jabba.simple.rds"))
    highcopyX = Ggraph[cn>CNxt]    
    Cycles_gW <- gW(snode.id=Cycles_List, graph=highcopyX, circular = TRUE)
    saveRDS(Cycles_gW, paste0("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/",CaseID,"/",CaseID,"_",CNxt,"_OptimizedCycles.rds"))

    nodes_list <- Cycles_gW$grl
    # Define MDM2 genomic coordinates
    mdm2_region <- GRanges(seqnames = "chr12", ranges = IRanges(start = 69202258, end = 69233629))
    seqlevelsStyle(mdm2_region) <- seqlevelsStyle(nodes_list)
    # Find overlaps between MDM2 and the GRangesList (nodes_list)
    mdm2_overlaps <- findOverlaps(nodes_list, mdm2_region, type = "any")
    # Use logical indexing to filter walks that overlap with MDM2
    mdm2_walks <- Cycles_gW[unique(queryHits(mdm2_overlaps))]
    nodes_list_mdm2 <- mdm2_walks$grl
    walk_widths <- sapply(nodes_list_mdm2, function(gr) sum(width(gr)))
    min_size <- 50000    # 50 kb
    max_size <- 5000000  # 5 mb
    valid_walk_indices <- which(walk_widths >= min_size & walk_widths <= max_size)
    filtered_walks <- Cycles_gW[valid_walk_indices]
    saveRDS(filtered_walks, paste0("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/",CaseID,"/",CaseID,"_",CNxt,"_filteredCycles.rds"))
    if (length(filtered_walks)>0){
    gencode <- track.gencode(cached.dir="/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics", stack.gap = 2e5, cex.label = 0.5, height = 2000)
    gencode$legend=FALSE                                                                                                                                                
  
    pdf(file = paste0('/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/',CaseID,'/Optimized_Walks_CS_', CNxt,'.pdf'),width=20)
    plot(c(gencode, Lipo863B_REP1, filtered_walks$gtrack(height=4000)), filtered_walks$footprint+2e5, chr.sub=FALSE)
    title(main=paste0("MDM2 circular walks with CN greater than ",CNxt))
    dev.off()
    }
  } 
})