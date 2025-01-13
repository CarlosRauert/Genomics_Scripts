library(gUtils)
library(gGnome)
library(gTrack)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(parallel)
library(cluster)
library(dplyr)
library(kmed)

gencode <- track.gencode(cached.dir="/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics", stack.gap = 2e5, cex.label = 0.5, height = 3000)
gencode$legend=FALSE

Lipo863B_REP1 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/ChIPseq/Lipo863B_mod_REP1.mLb.clN.bigWig",
                        name="Lipo836B", chr.sub=FALSE)
LP6_REP1 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/ChIPseq/LP6_mod_REP1.mLb.clN.bigWig",
                        name="LP6", chr.sub=FALSE)
LPS141 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/ChIPseq/LPS141_mod_REP1.mLb.clN.bigWig",
                        name="LPS141", chr.sub=FALSE)
SW872_REP1 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/ChIPseq/SW872_mod_REP1.mLb.clN.bigWig",
                        name="SW872_REP1", chr.sub=FALSE)        
SW872_REP2 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/ChIPseq/SW872_mod_REP2.mLb.clN.bigWig",
                        name="SW872_REP2", chr.sub=FALSE)
Tumor_1 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/ChIPseq/Tumor_1_mod_REP1.mLb.clN.bigWig",
                        name="Tumor 1", chr.sub=FALSE)
Tumor_2 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/ChIPseq/Tumor_2_mod_REP1.mLb.clN.bigWig",
                        name="Tumor 2", chr.sub=FALSE)

CaseID_l = c(
  "A1KU", "A1KW", "A1L0", "A1L2", "A1L3", "A240", "A2IZ", "A2J0", "A2J4", "A2QS", "A3LS", "A3LT", "A3LW", "A3LY",
  "CPCT02010386T", "CPCT02010680T", "CPCT02060104T", "CPCT02060191T", "CPCT02070051T", "CPCT02070366T", 
  "CPCT02080206T", "CPCT02080227T", "CPCT02090057T", "CPCT02340046T")
CNxt_l = c(
  12, 38, 9, 20, 38, 12, 20, 11, 17, 38, 29, 125, 8, 14,
  91, 11, 36, 11, 9, 9, 9, 10, 18, 13)

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

i=15
i = 11
GetMedoidWalks <- function(i){
  CaseID=CaseID_l[i]
  CNxt=CNxt_l[i]
  if (!file.exists(paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID))) {
    dir.create(paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID))
  }
  print("reading file")
  Cycles_File=paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/MDM2_Walks_Out/",CaseID,"/",CaseID,"_",CNxt,"_Cycles.txt")
  print("file read")
  Python_Output <- processFile(Cycles_File)
  Cycles_List <- PythonLists2RLists(Python_Output)
  Ggraph <-  gG(jabba=paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/Jabba_Calls_TCGA/",CaseID,"/jabba.simple.rds"))
  highcopyX = Ggraph[cn>CNxt]
  nodes_list <- highcopyX$gr
  print(CaseID)
  print(CNxt)
  # Get only Cycles that contain MDM2
  print("filtering")
  # Define MDM2 genomic coordinates
  mdm2_region <- GRanges(seqnames = "chr12", ranges = IRanges(start = 69202258, end = 69233629))
  seqlevelsStyle(mdm2_region) <- seqlevelsStyle(nodes_list)
  # Find overlaps between MDM2 and the GRangesList (nodes_list)
  mdm2_overlaps <- findOverlaps(nodes_list, mdm2_region, type = "any")
  mdm2_plusStrand <- mdm2_overlaps@from[1:(length(mdm2_overlaps@from)/2)]
  mdm2_minusStrand <- lapply(mdm2_plusStrand, function(x) -x)
  mdm2nodes <- rbind(mdm2_plusStrand, mdm2_minusStrand)
  # Filter
  which_MDM2 <- sapply(Cycles_List, function(x) all(mdm2_minusStrand %in% x) | all(mdm2_plusStrand %in% x))
  Cycles_MDM2 <- Cycles_List[which_MDM2]

  # Get only Cycles that fulfill width criteria
  # Get list of node lengths in bp
  Widths=width(nodes_list)
  names(Widths)=nodes_list$snode.id
  # Find Cycle Indices that fulfill criteria
  CycleLengths_List <- lapply(1:length(Cycles_MDM2), function(i){
    WidSum <- sum(Widths[as.character(Cycles_MDM2[[i]])])
    return(WidSum)
  })
  # Filter
  Cycles_Wid <- Cycles_MDM2[which(CycleLengths_List>50000 & CycleLengths_List<5000000)]

  # Prepare binary feature matrix of Nodes in Ggraph
  feature_Mat <- lapply(1:length(Cycles_Wid), function(i){
    binary <- as.numeric(ifelse(nodes_list$snode.id %in% Cycles_Wid[[i]], 1, 0))
    return(binary)
  })
  feature_Mat <- do.call(rbind, feature_Mat)
  rownames(feature_Mat)<-as.character(1:length(Cycles_Wid))

  # Conduct K-medoid Clustering
  data <- feature_Mat
  logmat <- matrix(as.logical(data), nrow=nrow(data), ncol=ncol(data))
  print("computing distance")
  coocur_dist <- cooccur(logmat)
  if (!file.exists(paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/fastkmed_silplots/"))) {
    dir.create(paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/fastkmed_silplots/"))
  }
  asws <- numeric(19)
  for (i in 2:20) {
      print(i)
      asw <- tryCatch({
          km <- fastkmed(coocur_dist, ncluster = i)
          silkmi <- sil(coocur_dist, km$medoid, km$cluster)
          dir_path <- paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/", CaseID, "/fastkmed_silplots/")
          if (!file.exists(dir_path)) {
              dir.create(dir_path, recursive = TRUE)
          }
          pdf(paste0(dir_path, "silplot_k", i, ".pdf"), width = 7, height = 5)
          print(silkmi$plot)
          dev.off()
          mean(silkmi$result$silhouette)
      }, error = function(e) {
          message("Error at k = ", i, ": ", e$message)
          return(0)
      })
      asws[i - 1] <- asw
  }
  ks <- 2:20
  pdf(paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/fastkmed_optimalK.pdf"))
  plot(ks, as.numeric(unlist(asws)), type="b", xlab="k", ylab="avg sil width")
  dev.off()
  optik <- which.max(asws)+1
  optikm <- fastkmed(coocur_dist, ncluster=optik)
  print(paste0("optik: ",optik))
  # Plot and save Medoid walks
  medoid_walks <- Cycles_Wid[optikm$medoid]
  Cycles_gW <- gW(snode.id=medoid_walks, graph=highcopyX, circular = TRUE)
  saveRDS(Cycles_gW, paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/fastkmed_medoidwalks.rds"))
  pdf(file = paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/fastkmed_medoidwalks.pdf"))
  plot(c(gencode, Cycles_gW$gtrack(height=4000)), Cycles_gW$footprint+2e5, chr.sub=FALSE)
  title(main=paste0("Medoid Walks after K-medoid clustering, \n CaseID: ",CaseID,",\n min CN: ",CNxt))
  dev.off()

  # return number of clusters
  return(optik)
}

optiklist <- mclapply(1:24, GetMedoidWalks, mc.cores=32)