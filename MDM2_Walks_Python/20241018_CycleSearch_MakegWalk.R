library(gGnome)
library(gTrack)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(readxl)
library(parallel)
Cycles_L <- list("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1KU/optimized_cycles_12.txt",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1KW/optimized_cycles_38.txt",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1L0/optimized_cycles_36.txt",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1L2/optimized_cycles_20.txt",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1L3/optimized_cycles_48.txt",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A240/optimized_cycles_15.txt",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A2IZ/optimized_cycles_20.txt",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A2J0/optimized_cycles_14.txt",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A23R/optimized_cycles_55.txt",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A2J4/optimized_cycles_13.txt",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A2QS/optimized_cycles_38.txt",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A3LS/optimized_cycles_29.txt",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A3LT/optimized_cycles_125.txt",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A3LW/optimized_cycles_8.txt",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A3LY/optimized_cycles_14.txt")

CaseID_l = list("A1KU", "A1KW", "A1L0", "A1L2", "A1L3", "A240", "A2IZ", "A2J0", "A23R", "A2J4", "A2QS", "A3LS", "A3LT", "A3LW", "A3LY")
CNxt_l = list("12", "38", "36", "20", "48", "15", "20", "14", "55", "13", "38", "29", "125", "8", "14")

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

xGr=1

lapply(1:14, function(xGr){
  Cycles_File=Cycles_L[[xGr]][1]
  CaseID <- CaseID_l[xGr]
  CNxt <- CNxt_l[xGr]

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
})



