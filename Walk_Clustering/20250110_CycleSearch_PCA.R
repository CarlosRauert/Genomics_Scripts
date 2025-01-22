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
library(logisticPCA)
library(stats)

gencode <- track.gencode(cached.dir="/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics", stack.gap = 2e5, cex.label = 0.3, height = 4000)
gencode$legend=FALSE

Lipo863B_REP1 <- gTrack("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/ChIPseq/Lipo863B_mod_REP1.mLb.clN.bigWig",
                        name="Lipo836B", chr.sub=FALSE, height=3000)
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

i=1
i=15
i = 11
i=2
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
  fesature_Mat <- do.call(rbind, feature_Mat)
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

test <- list(c(1,1,1),c(2,2,2),c(1,1,1),c(2,2,2),c(3,3,3))



# conduct logistic pca

# remove mirror walks
names(Cycles_Wid) <- 1:length(Cycles_Wid)
abs_Cycles <- lapply(Cycles_Wid, abs)
abs_Cycles <- lapply(abs_Cycles, sort)
my_list <- abs_Cycles
# Identify duplicate values
duplicate_indices <- which(!duplicated(my_list))
filtered_Wid <- Cycles_Wid[duplicate_indices]

# Prepare binary feature matrix of Nodes in Ggraph
feature_Mat <- lapply(1:length(filtered_Wid), function(i){
  binary <- as.numeric(ifelse(nodes_list$snode.id %in% filtered_Wid[[i]], 1, 0))
  return(binary)
})
feature_Mat <- do.call(rbind, feature_Mat)
rownames(feature_Mat)<-as.character(1:length(filtered_Wid))

# Print the results
log_pca <- logisticPCA(feature_Mat, k=10)
pca_scores <- log_pca$PCs
set.seed(123)
kmeans_result <- kmeans(pca_scores, centers = 7, iter.max=999)
clusters <- kmeans_result$cluster
plot_data <- data.frame(PC1 = pca_scores[, 1], PC2 = pca_scores[, 2], Cluster = as.factor(clusters), Index = 1:nrow(pca_scores))
# Generate the plot
pdf(paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/logpca_kmeans.pdf"))  
print(ggplot(plot_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +
  geom_text(aes(label = Index), vjust = -1, size = 3) +
  labs(title = "K-means Clusters on First Two Principal Components",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal())
dev.off()

i=2
asws <- list()
# silhouette method
for (i in 2:800){
  print(i)
  ikmeans_result <- kmeans(pca_scores[,c(1,2)], centers=i)
  silhouette_values <- silhouette(ikmeans_result$cluster, dist(pca_scores))
  asw <- mean(silhouette_values[,"sil_width"])
  asws <- append(asws, asw)
}

set.seed(123)
kmeans_result <- kmeans(pca_scores[,c(1,2)], centers = which.max(asws))
clusters <- kmeans_result$cluster
plot_data <- data.frame(PC1 = pca_scores[, 1], PC2 = pca_scores[, 2], Cluster = as.factor(clusters), Index = 1:nrow(pca_scores))
# Generate the plot
pdf(paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/logpca_kmeans.pdf"))  
print(ggplot(plot_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +
  geom_text(aes(label = Index), vjust = -1, size = 3) +
  labs(title = "K-means Clusters on First Two Principal Components",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal())
dev.off()

pdf(paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/pcakmeans_optimalK.pdf"))
plot(2:316, as.numeric(unlist(asws)), type="b", xlab="k", ylab="avg sil width")
dev.off()

















# compute pca on absolute node ids while keeping original IDs
names(Cycles_Wid) <- 1:length(Cycles_Wid)
abs_Cycles <- lapply(Cycles_Wid, abs)
abs_Cycles <- lapply(abs_Cycles, sort)
nonduplicate_indices <- which(!duplicated(abs_Cycles))
filtered <- abs_Cycles[nonduplicate_indices]
# get feature mat
feature_Mat <- lapply(1:length(filtered), function(i){
  binary <- as.numeric(ifelse(unique(abs(nodes_list$snode.id)) %in% filtered[[i]], 1, 0))
  return(binary)
})
feature_Mat <- do.call(rbind, feature_Mat)
rownames(feature_Mat)<-as.character(1:length(filtered))
# do pca
log_pca <- logisticPCA(feature_Mat, k=10)
pca_scores <- log_pca$PCs
# silhouette method
set.seed(123)
asws <- list()
for (i in 2:30){
  print(i)
  ikmeans_result <- kmeans(pca_scores[,1:2], centers=i, iter.max=400000)
  silhouette_values <- silhouette(ikmeans_result$cluster, dist(pca_scores))
  asw <- mean(silhouette_values[,"sil_width"])
  asws <- append(asws, asw)
}
which.max(asws)
write.csv(pca_scores[,1:2], paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/logpca_PC1_2.csv"))

pdf(paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/logpca_kmeans__Silhouettes.pdf"))
plot(2:(length(asws)+1),asws, type="l", xlab="number of clusters", ylab="silhouette coefficient")
dev.off()

kmeans_result <- kmeans(pca_scores[,1:2], centers = 27)
clusters <- kmeans_result$cluster
plot_data <- data.frame(PC1 = pca_scores[, 1], PC2 = pca_scores[, 2], Cluster = as.factor(clusters), Index = 1:nrow(pca_scores))
pdf(paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/logpca_kmeans_.pdf"))  
print(ggplot(plot_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +
  geom_text(aes(label = Index), vjust = -1, size = 3) +
  labs(title = "K-means Clusters on First Two Principal Components",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal())
dev.off()

# Calculate distances to cluster centers
distances <- sapply(1:nrow(pca_scores), function(i) {
  cluster <- kmeans_result$cluster[i]
  centroid <- kmeans_result$centers[cluster, ]
  sum((pca_scores[i, 1:2] - centroid)^2)
})

# Find the elements with the least distance to their cluster centers
closest_elements <- sapply(1:max(kmeans_result$cluster), function(cluster) {
  cluster_indices <- which(kmeans_result$cluster == cluster)
  cluster_distances <- distances[cluster_indices]
  cluster_indices[which.min(cluster_distances)]
})

# closest_elements now contains the indices of the elements closest to their cluster centers
print(closest_elements)

# plot Gwalks of closest elements
center_names <- names(filtered)[closest_elements]
center_walks <- Cycles_Wid[as.integer(center_names)]
Centers_gW <- gW(snode.id=center_walks, graph=highcopyX, circular = TRUE)
saveRDS(Centers_gW, paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/kmeans_centerwalks.rds"))
pdf(file = paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/kmeans_centerwalks.pdf"), height=40)
plot(c(gencode, Centers_gW$gtrack(height=7000), Lipo863B_REP1), Centers_gW$footprint+2e5, chr.sub=FALSE)
title(main=paste0("Center Walks after K-Means clustering, \n CaseID: ",CaseID,",\n min CN: ",CNxt))
dev.off()

# when running with i=11, you get A3LS. when running with i=15 you get  the first CPCT


# visualize only Cluster 1

cluster1_names=names(filtered)[which(kmeans_result$cluster==1)]
cluster1_walks <- Cycles_Wid[as.integer(cluster1_names)]
cluster1_gW <- gW(snode.id=cluster1_walks, graph=highcopyX, circular = TRUE)
pdf(file = paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/kmeans_cluster1walks.pdf"))
plot(cluster1_gW$gtrack(height=7000), cluster1_gW$footprint+2e5, chr.sub=FALSE)
title(main=paste0("Cluster 1 Walks after K-Means clustering, \n CaseID: ",CaseID,",\n min CN: ",CNxt))
dev.off()

# visualize only Cluster 2

cluster2_names=names(filtered)[which(kmeans_result$cluster==2)]
cluster2_walks <- Cycles_Wid[as.integer(cluster2_names)]
cluster2_gW <- gW(snode.id=cluster2_walks, graph=highcopyX, circular = TRUE)
pdf(file = paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/kmeans_cluster2walks.pdf"))
plot(cluster2_gW$gtrack(height=7000), cluster2_gW$footprint+2e5, chr.sub=FALSE)
title(main=paste0("Cluster 2 Walks after K-Means clustering, \n CaseID: ",CaseID,",\n min CN: ",CNxt))
dev.off()

# visualize only Cluster 3

cluster3_names=names(filtered)[which(kmeans_result$cluster==3)]
cluster3_walks <- Cycles_Wid[as.integer(cluster3_names)]
cluster3_gW <- gW(snode.id=cluster3_walks, graph=highcopyX, circular = TRUE)
pdf(file = paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/kmeans_cluster3walks.pdf"))
plot(cluster3_gW$gtrack(height=7000), cluster3_gW$footprint+2e5, chr.sub=FALSE)
title(main=paste0("Cluster 3 Walks after K-Means clustering, \n CaseID: ",CaseID,",\n min CN: ",CNxt))
dev.off()

#Visualize all walks

all_walks<- gW(snode.id=Cycles_Wid, graph=highcopyX, circular = TRUE)
pdf(file = paste0("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/",CaseID,"/allwalks.pdf"), height=20)
plot(all_walks$gtrack(height=30000), all_walks$footprint+2e5, chr.sub=FALSE)
title(main=paste0("Circular walks returned by Johnson's algorithm, \n CaseID: ",CaseID,",\n min CN: ",CNxt))
dev.off()