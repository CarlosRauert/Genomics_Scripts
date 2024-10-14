# Required Libraries
library(gTrack)
library(GenomicRanges)
library(parallel)
library(dplyr)
#https://genome.ucsc.edu/cgi-bin/hgGene?db=hg19&hgg_gene=MDM2

# Function to determine the largest Walk from each group

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


Properties_List <- mclapply(fileList,function(gw_file){
    gw <- readRDS(gw_file)
    # Assuming `gw` is your gWalk object
    # Extract the segments/nodes of each walk as a GRangesList
    nodes_list <- gw$grl  # Extract the nodes from the gWalk object, a GRangesList
    nwalks_raw <- length(nodes_list)
    # Assuming `gw` is your gWalk object and MDM2 genomic coordinates are known
    mdm2_region <- GRanges(seqnames = "chr12", ranges = IRanges(start = 69202258, end = 69233629))
    seqlevelsStyle(mdm2_region) <- seqlevelsStyle(nodes_list)
    # Subset the gWalk object to keep only those walks that overlap with MDM2
    # Find overlaps between MDM2 and the GRangesList (nodes_list)
    mdm2_overlaps <- lapply(nodes_list, function(gr) {
    # Ensure each element of nodes_list is a GRanges
    if (is(gr, "GRanges")) {
        hits <- findOverlaps(gr, mdm2_region, type = "any")
        return(length(hits) > 0)  # Return TRUE if there are overlaps
    } else {
        return(FALSE)
    }
    })    # Filter the walks based on whether they overlap with MDM2
    mdm2_walks <- gw[unlist(mdm2_overlaps)] 
    nodes_list_mdm2 <- mdm2_walks$grl  # Extract the nodes from the gWalk object, a GRangesList
    # Calculate the widths of each walk
    nwalks_mdm2 <- length(nodes_list_mdm2)
    ratio=nwalks_mdm2/nwalks_raw
    return(list(nwalks_raw, ratio))
}, mc.cores=4)

ratio_list <- lapply(Properties_List, function(x){
    ratioX<-x[[2]]
    return(ratioX)
})
pdf(file="data/cephfs-1/home/users/rauertc_c/work/ratioPlot.pdf")
hist(ratio_list)
dev.off()

nwalks_list <- lapply(Properties_List, function(x){
    nwalkX<-x[[1]]
    return(nwalkX)
})

pdf(file="data/cephfs-1/home/users/rauertc_c/work/nwalksPlot.pdf")
hist(log10(nwalks))
dev.off()



ratios_list <- lapply
Get_Gwalk_SimilarGroups <- function(gw_file){
    gw <- readRDS(gw_file)
    # Assuming `gw` is your gWalk object
    # Extract the segments/nodes of each walk as a GRangesList
    nodes_list <- gw$grl  # Extract the nodes from the gWalk object, a GRangesList
    # Assuming `gw` is your gWalk object and MDM2 genomic coordinates are known
    mdm2_region <- GRanges(seqnames = "chr12", ranges = IRanges(start = 69202258, end = 69233629))
    seqlevelsStyle(mdm2_region) <- seqlevelsStyle(nodes_list)
    # Subset the gWalk object to keep only those walks that overlap with MDM2
    # Find overlaps between MDM2 and the GRangesList (nodes_list)
    mdm2_overlaps <- lapply(nodes_list, function(gr) {
    # Ensure each element of nodes_list is a GRanges
    if (is(gr, "GRanges")) {
        hits <- findOverlaps(gr, mdm2_region, type = "any")
        return(length(hits) > 0)  # Return TRUE if there are overlaps
    } else {
        return(FALSE)
    }
    })    # Filter the walks based on whether they overlap with MDM2
    mdm2_walks <- gw[unlist(mdm2_overlaps)] 
    nodes_list_mdm2 <- mdm2_walks$grl  # Extract the nodes from the gWalk object, a GRangesList
    # Calculate the widths of each walk
    walk_widths <- sapply(nodes_list_mdm2, function(gr) sum(width(gr)))
    # Define the size range in base pairs
    min_size <- 50000    # 50 kb
    max_size <- 5000000  # 5 mb
    # Identify walks that are within the specified size range
    valid_walk_indices <- which(walk_widths >= min_size & walk_widths <= max_size)
    # Subset the gWalk object to keep only those walks within the size range
    filtered_walks <- gw[valid_walk_indices]
    # Check if filtered_walks is empty
    if (length(filtered_walks) == 0) {
        stop("No walks found within the specified size range.")
    }
    # Flatten the GRangesList into a GRanges object
    gr_flat <- unlist(filtered_walks$grl)
    # Calculate the average width of the ranges in gr_flat
    median_width <- median(sum(width(nodes_list_mdm2)))
    # Set minoverlap to 90% of the average width
    minoverlap_value <- floor(0.90 * median_width)
    # We calculate 90% overlap by using the minoverlap argument in findOverlaps
    overlaps <- findOverlaps(nodes_list_mdm2, nodes_list_mdm2, minoverlap=minoverlap_value, type="any")

    # Convert the overlaps to a list of groups
    overlap_groups <- split(subjectHits(overlaps), queryHits(overlaps))

    get_largest_walk <- function(walk_indices) {
        # Check if the group is not empty and contains valid indices
        if (length(walk_indices) == 0) {
            return(NA)  # Return NA if the group is empty
        }
        69,202,258-69,233,629
        # Ensure that all indices are within bounds of the GRangesList
        valid_indices <- walk_indices[walk_indices <= length(nodes_list)]
  
        if (length(valid_indices) == 0) {
            return(NA)  # Return NA if no valid indices are found
        }
        # Find the largest walk based on the total width of the nodes
        walk_sizes <- sapply(valid_indices, function(i) sum(width(nodes_list[[i]])))
        largest_walk_index <- valid_indices[which.max(walk_sizes)]
  
        return(largest_walk_index)
    }

    overlap_groups_subset <- 
    # Apply the function to all groups to retain only the largest Walk
    largest_walk_indices <- sapply(overlap_groups, get_largest_walk)
    
    # Subset the gWalk object to keep only the largest Walks
    subset_gw <- gw[unique(largest_walk_indices)]

    # The resulting gWalk object is now subsetted with only the largest Walk from each overlapping group
    return(subset_gw)
}

Get_Gwalk_SimilarGroups("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1KU/A1KU_12_Cycles.rds")

getMDM2