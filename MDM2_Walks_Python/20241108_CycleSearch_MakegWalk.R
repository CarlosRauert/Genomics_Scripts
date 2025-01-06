library(gGnome)
library(gTrack)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(readxl)
library(parallel)
library(vegan)

toInstall <- c("cluster", "fpc", "mclust")
install.packages(toInstall, dependencies=TRUE)
library(cluster)
library(fpc)
library(mclust)

install.packages("proxy")
library(proxy)
library(dplyr)
CaseID="A1KU"
Cycles_File="/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1KU/A1KU_12_Cycles.txt"
CNxt=12

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

nodes_list <- highcopyX$gr

# Get only Cycles that contain MDM2
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

i=1
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

i=1
# Prepary binary feature matrix of Nodes in Ggraph
feature_Mat <- lapply(1:length(Cycles_Wid), function(i){
  binary <- as.numeric(ifelse(nodes_list$snode.id %in% Cycles_Wid[[i]], 1, 0))
  return(binary)
})
feature_Mat <- do.call(rbind, feature_Mat)
rownames(feature_Mat)<-as.character(1:)
# Cluster Cycles based on Jaccard similarity and perform hierarchical clustering
jaccard_Distances <- dist(feature_Mat, method="binary")
hc <- hclust(jaccard_Distances, method="average")
pdf(file="/data/cephfs-1/home/users/rauertc_c/work/genomics/20241206_Jaccard/20241206_Jaccard_Den.pdf")
plot(hc, main="Dendrogram of TCGA-A1KU cyclic DNA",ylab="Jaccard Dissimilarity",xlab="cycles detected", labels=FALSE)
dev.off()

dend <- as.dendrogram(hc)
# Function to recursively extract heights
extract_heights <- function(node) {
  if (is.leaf(node)) {
    return(NULL)  # Leaves don't have branches
  } else {
    return(c(attr(node, "height"), sapply(node, extract_heights)))
  }
}
# Get all branch heights
branch_heights <- unlist(extract_heights(dend))
print(branch_heights)
# Sort heights and count unique occurrences
height_counts <- table(branch_heights)
# Cumulative number of branches with respect to height
cumulative_branches <- cumsum(as.numeric(height_counts))
cumulative_branches <- rep(1592,414)-cumulative_branches
# Create a data frame for plotting
height_data <- data.frame(
  height = as.numeric(names(height_counts)),
  branches = cumulative_branches
)
# Plot cumulative branches vs. height
pdf(file="/data/cephfs-1/home/users/rauertc_c/work/genomics/20241206_Jaccard/20241206_Jaccard_elbow.pdf")
print(ggplot(height_data, aes(x = height, y = branches)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Cumulative Number of Branches vs. Height",
    x = "Height",
    y = "Number of Branches"
  ) +
  theme_minimal())
dev.off()




k <- 8  # Number of clusters
clusters <- cutree(hc, k)
pdf(file="/data/cephfs-1/home/users/rauertc_c/work/genomics/20241206_Jaccard/20241206_Jaccard_Den_Rect.pdf")
plot(hc, ylab="Jaccard Dissimilarity", xlab="cycles detected", labels=FALSE)
rect.hclust(hc, k = k, border = "red")  # Highlight the clusters
dev.off()

saveRDS(filtered_walks, paste0("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/",CaseID,"/",CaseID,"_",CNxt,"_filteredCycles.rds"))
if (length(filtered_walks)>0){
  gencode <- track.gencode(cached.dir="/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics", stack.gap = 2e5, cex.label = 0.5, height = 2000)
  gencode$legend=FALSE                                                                                                                                                
  
  pdf(file = paste0('/data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/MDM2_Walks_Out/',CaseID,'/Optimized_Walks_CS_', CNxt,'.pdf'),width=20)
  plot(c(gencode, Lipo863B_REP1, filtered_walks$gtrack(height=4000)), filtered_walks$footprint+2e5, chr.sub=FALSE)
  title(main=paste0("MDM2 circular walks with CN greater than ",CNxt))
  dev.off()
}



