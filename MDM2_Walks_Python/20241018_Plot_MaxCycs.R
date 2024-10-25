# Required Libraries
library(gTrack)
library(GenomicRanges)
library(parallel)
library(dplyr)
#https://genome.ucsc.edu/cgi-bin/hgGene?db=hg19&hgg_gene=MDM2

# Function to determine the largest Walk from each group

dir.create("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/Optimized_Cycles")

fileList=list("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1KU/A1KU_12_Cycles.rds",
            "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1KW/A1KW_38_Cycles.rds",
            "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1L0/A1L0_9_Cycles.rds",
            "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1L2/A1L2_20_Cycles.rds",
            "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1L3/A1L3_38_Cycles.rds",
            "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A240/A240_12_Cycles.rds",
            "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A2IZ/A2IZ_20_Cycles.rds",
            "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A2J0/A2J0_11_Cycles.rds",
            "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A2J4/A2J4_17_Cycles.rds",
            "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A2QS/A2QS_38_Cycles.rds",
            "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A3LS/A3LS_29_Cycles.rds",
            "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A3LT/A3LT_125_Cycles.rds",
            "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A3LW/A3LW_8_Cycles.rds",
            "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A3LY/A3LY_14_Cycles.rds")

i=1

Get_Gwalk_SimilarGroups <- function(i) {
  print(i)
  gw_file <- fileList[i]
  gw <- readRDS(gw_file[[1]])
  # Extract the segments/nodes of each walk as a GRangesList
  nodes_list <- gw$grl
  # Define MDM2 genomic coordinates
  mdm2_region <- GRanges(seqnames = "chr12", ranges = IRanges(start = 69202258, end = 69233629))
  seqlevelsStyle(mdm2_region) <- seqlevelsStyle(nodes_list)
  # Find overlaps between MDM2 and the GRangesList (nodes_list)
  mdm2_overlaps <- findOverlaps(nodes_list, mdm2_region, type = "any")
  # Use logical indexing to filter walks that overlap with MDM2
  mdm2_walks <- gw[unique(queryHits(mdm2_overlaps))]
  nodes_list_mdm2 <- mdm2_walks$grl
  # Calculate the widths of each walk
  walk_widths <- sapply(nodes_list_mdm2, function(gr) sum(width(gr)))
  # Define the size range in base pairs
  min_size <- 50000    # 50 kb
  max_size <- 5000000  # 5 mb
  # Identify walks that are within the specified size range
  valid_walk_indices <- which(walk_widths >= min_size & walk_widths <= max_size)
  # Check if valid walks are found
  if (length(valid_walk_indices) == 0) {
    stop("No walks found within the specified size range.")
  }
  # Subset to keep only valid walks
  filtered_walks <- gw[valid_walk_indices]
  nodes_list_mdm2 <- filtered_walks$grl
  # Flatten the GRangesList into a GRanges object
  gr_flat <- unlist(nodes_list_mdm2)
  # Calculate the median width
  median_width <- median(sum(width(nodes_list_mdm2)))
  # Set minoverlap to 90% of the average width
  minoverlap_value <- floor(0.90 * median_width)
  # Find overlaps among valid walks
  print("overlap check")
  overlaps <- findOverlaps(nodes_list_mdm2, nodes_list_mdm2, minoverlap = minoverlap_value, type = "any")
  # Convert the overlaps to a list of groups
  overlap_groups <- split(subjectHits(overlaps), queryHits(overlaps))
  # Function to find the largest walk
  get_largest_walk <- function(walk_indices) {
    if (length(walk_indices) == 0) return(NA)
    # Calculate the widths for valid indices only
    walk_sizes <- sapply(walk_indices, function(i) sum(width(nodes_list_mdm2[[i]])))
    # Return the index of the largest walk
    largest_walk_index <- walk_indices[which.max(walk_sizes)]
    return(largest_walk_index)
  }
  # Apply the function to all groups to retain only the largest Walk
  largest_walk_indices <- sapply(overlap_groups, get_largest_walk, USE.NAMES = FALSE)
  # Subset the gWalk object to keep only the largest Walks
  subset_gw <- gw[unique(largest_walk_indices)]
  saveRDS(subset_gw,paste0("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/Optimized_Cycles/",i,"_cycles.rds"))
  # Explicitly call garbage collection
  gc()
  # save the subsetted gWalk object
  
}

Subsets <- lapply(1:14, Get_Gwalk_SimilarGroups)d