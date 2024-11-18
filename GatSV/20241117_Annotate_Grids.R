system("export PKG_CONFIG_PATH=/data/cephfs-1/home/users/rauertc_c/work/miniforge3/envs/Gatsv/lib/pkgconfig:$PKG_CONFIG_PATH")

library(BiocGenerics)
library(caTools)
library(data.table)
library(e1071)
library(GenomeInfoDb)
library(GenomicRanges)
library(gUtils)
library(IRanges)
library(parallel)
library(rlang)
library(ROCR)
library(S4Vectors)
library(stats4)
library(stringr)
library(here)
library(optparse)
library(R.utils)

#ANNOTATE 
cat('Loading reference files...\n')
gnomad_hg38 = readRDS('/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/gnomAD.v4.hg38.rds')
gnomad_hg19 = readRDS('/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/gnomAD.v4.hg19.liftover.rds')
LINE_dt_hg38 = readRDS('/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/repeatmasker.hg38.LINE.bed')
SINE_dt_hg38 = readRDS('/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/repeatmasker.hg38.SINE.bed')
LINE_dt_hg19 = readRDS('/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/repeatmasker.hg19.LINE.bed')
SINE_dt_hg19 = readRDS('/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/repeatmasker.hg19.SINE.bed')
hg19_genes = readRDS('/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/gencode.genes.hg19.rds')
hg19_exons=readRDS('/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/gencode.exons.hg19.rds')
hg38_genes=readRDS('/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/gencode.genes.hg38.rds')
hg38_exons=readRDS('/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/gencode.exons.hg38.rds')
reptimedata_hg19 = readRDS('/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/reptime.hg19.rds')
reptimedata_hg38 = readRDS('/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/reptime.hg38.rds')
scaling_mat <- fread("/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/scalingmatrix.txt")

GaTSV <- readRDS("/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/svm/GaTSV.rda") #svmobject


#running the classifier on example data
metadata <- fread("/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/example_metadata.txt") #metadata file that contains the sample_ids (same as 'sample' input) and associated tp53_mutation_status
file_path <- "/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/example.sv.vcf" #replace with desired vcf path
sample <- "example"

run_GaTSV(file_path,sample,n_cores=16,genome='hg19',output_path = '/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out')

process_file(file_path,sample,n_cores=16,genome="hg19",output_path='/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/Example')



##Two output files are generated and stored in the output_path provided under the names:
#'filename'_processed.bedpe 
#'filename'_classified.bedpe
#"filename" is the basename of the vcf file provided without the extension

vcf_path_ex="/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/example.sv.vcf"
vcf_path="/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs/ACTN01020001T.purple.sv.vcf.gz"

vcf_to_dt <- function(vcf_path,sample_name) {
  cat(paste0(vcf_path, "\n"))
  if (!file.exists(vcf_path)) {
    print(paste("File does not exist",vcf_path))
  }
  
  cat("Reading file...\n")
  vcf_dt <- fread(cmd=paste("grep -v '^#'", vcf_path),sep='\t')
  
  # Set colnames of vcf_dt to standard...
  if (nrow(vcf_dt) == 0) {
    return (vcf_dt)
  }
  if (ncol(vcf_dt)==10) {
    setnames(vcf_dt, c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO","GENO","TUMOR"))
  } else {
    setnames(vcf_dt, c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO","GENO","NORMAL","TUMOR"), skip_absent=TRUE)
  }
  
  cat("Gathering Metadata...\n")
  if ("INFO" %in% colnames(vcf_dt) ) {
    vcf_dt[, SPAN := as.numeric(gsub(".*?SPAN=([-0-9]+).*","\\1",INFO))]
    vcf_dt$sample = sample_name
    vcf_dt[, uid := gsub("([0-9]+):(1|2)", "\\1", ID)]
    vcf_dt[, EVDNC := gsub(".*?EVDNC=([A-Z]+).*", "\\1", INFO)]
    vcf_dt[, MAPQ := as.integer(gsub(".*?;MAPQ=([0-9]+).*", "\\1", INFO))]
    vcf_dt[, HOMSEQ := gsub(".*?;HOMSEQ=([A-Z]+).*", "\\1", INFO)] 
    vcf_dt[, HOMSEQ := ifelse(grepl(";", HOMSEQ), "", HOMSEQ)] 
    vcf_dt[, INSERTION := gsub(".*?;INSERTION=([A-Z]+).*", "\\1", INFO)] 
    vcf_dt[, INSERTION := ifelse(grepl(";", INSERTION), "", INSERTION)]
    vcf_dt[, NDISC := as.numeric(gsub(".*?NDISC=([0-9]+).*", "\\1", INFO))]
    vcf_dt[, SVMETHOD := substr(INFO,regexpr("SVMETHOD=",INFO)+nchar("SVMETHOD="),regexpr(";NDISC",INFO)-1)]
    
  }
  
  # More extraction regexpr stuff...
  if ("TUMOR" %in% colnames(vcf_dt)) {
    vcf_dt[, TUMALT :=  as.integer(strsplit(TUMOR, ":")[[1]][2]) , by=uid]
    vcf_dt[, TUMCOV :=  as.integer(strsplit(TUMOR, ":")[[1]][3]) , by=uid]
    vcf_dt[, TUMLOD :=  as.numeric(strsplit(TUMOR, ":")[[1]][9]) , by=uid]
  }
  if ("NORMAL" %in% colnames(vcf_dt)) {
    vcf_dt[, NORMCOV :=  as.integer(strsplit(NORMAL, ":")[[1]][3]) , by=uid]
    vcf_dt[, NORMALT :=  as.integer(strsplit(NORMAL, ":")[[1]][2]) , by=uid]
    vcf_dt[, NORMLOD :=  as.numeric(strsplit(NORMAL, ":")[[1]][9]) , by=uid]
  }
  
  cat("Cleaning up...\n")
  vcf_dt[, strand := ifelse(grepl("^\\[", ALT) | grepl("^\\]", ALT), '-', '+')]
  vcf_dt[, inv := strand[1] == strand[2], by=uid]
  vcf_dt[, altstrand := rev(strand), by=uid]
  vcf_dt[, altpos := as.integer(gsub(".*?:([0-9]+).*", "\\1", ALT))]
  vcf_dt[, altchr := gsub(".*?(\\[|\\])(.*?):([0-9]+).*", "\\2", ALT)]
  vcf_dt[, end := start] 
  
  canonical_contigs <- c(c(1:24),c('X','Y'),paste0('chr',c(1:24)),paste0('chr',c('X','Y')))
  bad.ix <- vcf_dt[!(seqnames %in% canonical_contigs), uid] #modify this to only include canonical chromosomes
  vcf_dt <- vcf_dt[!uid %in% bad.ix]
  vcf_dt[, sid := sample_name]
  vcf_dt[, seqnames:= ifelse(grepl('chr',seqnames),seqnames,paste0("chr",seqnames)), by=uid]
  return(vcf_dt)
}

Example_dt <- vcf_to_dt(vcf_path_ex, "Example")

# Extract the first few rows
dt_head <- head(dt, 10)

# Create a table graphical object (grob) from the head of the data.table
table_grob <- tableGrob(dt_head)

# Specify the path for the PDF file
pdf("data_table_head.pdf", width = 8, height = 6)

# Draw the table on the PDF
grid.draw(table_grob)

# Close the PDF device to save the file
dev.off()



vcf_to_dt_Gridss <- function(vcf_path,sample_name) {
  cat(paste0(vcf_path, "\n"))
  if (!file.exists(vcf_path)) {
    print(paste("File does not exist",vcf_path))
  }
  # Read .gz file as data.table
  cat("Reading file...\n")
  con <- gzfile(vcf_path, "r")
  lines <- character()
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    if (!grepl("^#", line)) {
      lines <- c(lines, line)
    }
  }
  close(con)
  vcf_list <- strsplit(lines, "\t")  
  vcf_dt <- rbindlist(lapply(vcf_list, as.list), fill = TRUE, use.names=FALSE)  # Combine into a data.table

  # Set colnames of vcf_dt to standard...
  if (nrow(vcf_dt) == 0) {
    return (vcf_dt)
  }
  if (ncol(vcf_dt)==10) {
    setnames(vcf_dt, c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO","GENO","TUMOR"))
  } else {
    setnames(vcf_dt, c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO","GENO","NORMAL","TUMOR"), skip_absent=TRUE)
  }
  return(vcf_dt)
}

vcf_dt <- vcf_to_dt_Gridss("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs/ACTN01020001T.purple.sv.vcf.gz", "Test")
dt_example <- vcf_to_dt
i=1

get_metadata <- function(vcf_dt){
  cat("Gathering Metadata...\n")
  #calculate span
  spans <- lapply(1:nrow(vcf_dt), function(i){
    start_chrom=vcf_dt$seqnames[i]
    if (grepl("[", vcf_dt$ALT[i], fixed = TRUE) | grepl("]", vcf_dt$ALT[i], fixed = TRUE)){
      end_chrom=gsub(".*?([0-9XYMT]+):.*", "\\1", vcf_dt$ALT[i])
    }
    else{
      return(NA)
    }
    if (start_chrom==end_chrom){
      start_pos=as.numeric(vcf_dt$start[i])
      end_pos=as.numeric(gsub(".*:(\\d+).*", "\\1", vcf_dt$ALT[i]))
      span=abs(end_pos-start_pos)
      return(span)
    }
    else{
      return(-1)
    }
  })
  # get uid
  uids <- mclapply
  if ("INFO" %in% colnames(vcf_dt) ) {
    #vcf_dt[, SPAN := as.numeric(gsub(".*?SPAN=([-0-9]+).*","\\1",INFO))]
    vcf_dt[, SPAN := as.numeric(spans)]
    vcf_dt$sample = sample_name
    #vcf_dt[, uid := gsub("([0-9]+):(1|2)", "\\1", ID)]
    vcf_dt[, uid := gsub("([0-9]+):(o|h)", "\\1", ID)]
    vcf_dt[, EVDNC := gsub(".*?EVDNC=([A-Z]+).*", "\\1", INFO)]
    vcf_dt[, MAPQ := as.integer(gsub(".*?;MAPQ=([0-9]+).*", "\\1", INFO))]
    vcf_dt[, HOMSEQ := gsub(".*?;HOMSEQ=([A-Z]+).*", "\\1", INFO)] 
    vcf_dt[, HOMSEQ := ifelse(grepl(";", HOMSEQ), "", HOMSEQ)] 
    vcf_dt[, INSERTION := gsub(".*?;INSERTION=([A-Z]+).*", "\\1", INFO)] 
    vcf_dt[, INSERTION := ifelse(grepl(";", INSERTION), "", INSERTION)]
    vcf_dt[, NDISC := as.numeric(gsub(".*?NDISC=([0-9]+).*", "\\1", INFO))]
    vcf_dt[, SVMETHOD := substr(INFO,regexpr("SVMETHOD=",INFO)+nchar("SVMETHOD="),regexpr(";NDISC",INFO)-1)]
    
  }
  
  # More extraction regexpr stuff...
  if ("TUMOR" %in% colnames(vcf_dt)) {
    vcf_dt[, TUMALT :=  as.integer(strsplit(TUMOR, ":")[[1]][2]) , by=uid]
    vcf_dt[, TUMCOV :=  as.integer(strsplit(TUMOR, ":")[[1]][3]) , by=uid]
    vcf_dt[, TUMLOD :=  as.numeric(strsplit(TUMOR, ":")[[1]][9]) , by=uid]
  }
  if ("NORMAL" %in% colnames(vcf_dt)) {
    vcf_dt[, NORMCOV :=  as.integer(strsplit(NORMAL, ":")[[1]][3]) , by=uid]
    vcf_dt[, NORMALT :=  as.integer(strsplit(NORMAL, ":")[[1]][2]) , by=uid]
    vcf_dt[, NORMLOD :=  as.numeric(strsplit(NORMAL, ":")[[1]][9]) , by=uid]
  }
  
  cat("Cleaning up...\n")
  vcf_dt[, strand := ifelse(grepl("^\\[", ALT) | grepl("^\\]", ALT), '-', '+')]
  vcf_dt[, inv := strand[1] == strand[2], by=uid]
  vcf_dt[, altstrand := rev(strand), by=uid]
  vcf_dt[, altpos := as.integer(gsub(".*?:([0-9]+).*", "\\1", ALT))]
  vcf_dt[, altchr := gsub(".*?(\\[|\\])(.*?):([0-9]+).*", "\\2", ALT)]
  vcf_dt[, end := start] 
  
  canonical_contigs <- c(c(1:24),c('X','Y'),paste0('chr',c(1:24)),paste0('chr',c('X','Y')))
  bad.ix <- vcf_dt[!(seqnames %in% canonical_contigs), uid] #modify this to only include canonical chromosomes
  vcf_dt <- vcf_dt[!uid %in% bad.ix]
  vcf_dt[, sid := sample_name]
  vcf_dt[, seqnames:= ifelse(grepl('chr',seqnames),seqnames,paste0("chr",seqnames)), by=uid]
  return(vcf_dt)
}




process_file <- function(file,sample, n_cores,genome,output_path='./') {
  
  filt_file <- suppressWarnings(filter(file, sample))
  
  bedpe <- filt_file
  bedpe[, chrom1 := gsub('chr','',chrom1)]
  bedpe[, chrom2 := gsub('chr','',chrom2)]
  bedpe$start1 <- as.numeric(bedpe$start1)
  bedpe$start2 <- as.numeric(bedpe$start2)
  bedpe$end1 <- as.numeric(bedpe$end1)
  bedpe$end2 <- as.numeric(bedpe$end2)
  
  fuzzy <- closest_germline(bp = bedpe, cores = n_cores, genome = genome)
  
  line_sine <- closest_line_sine(bp = fuzzy, genome = genome, cores = n_cores)
  
  cat("Finding distance to closest SV... \n")
  nearest_sv_dist <- rbindlist(mclapply(1:nrow(line_sine), find_closest_sv, line_sine,mc.cores=n_cores))
  cat("done. \n")
  
  cat("Finding no of SVs in 5Mbp window... \n")
  sv_annotated <- rbindlist(mclapply(1:nrow(nearest_sv_dist), count_sv_5mbp, nearest_sv_dist, mc.cores = n_cores))
  cat("done. \n")
  
  cat("Adding replication timing info... \n")
  reptime_added <- rbindlist(mclapply(1:nrow(sv_annotated), rep_time, sv_annotated, genome=genome, mc.cores = n_cores))
  cat("done. \n")
  
  chroms <-  c(as.character(c(1:22)),'X','Y')
  reptime_added$chrom1 <- as.character(reptime_added$chrom1)
  reptime_added$chrom2 <- as.character(reptime_added$chrom2)
  store_indices = which(!reptime_added$chrom1 %in% chroms| !reptime_added$chrom2 %in% chroms) #done to address the observation that some SVs were mapping to hpv
  bedpe_clean <- reptime_added[!store_indices,]
  bedpe_clean$start1 <- as.numeric(bedpe_clean$start1)
  bedpe_clean$start2 <- as.numeric(bedpe_clean$start2)
  
  cat("Performing gene/exon annotation...")
  bedpe_annot <- rbindlist(mclapply(1:nrow(bedpe_clean), annot_geneexon, bedpe_clean, genome=genome,mc.cores = n_cores)) 
  cat("done. \n")
  
  cat("Checking TP53 status...")
  tp53_added <- check_tp53(bedpe_annot)
  cat("done. \n")
  
  filename <- tp53_added$sample[1]
  write.table(tp53_added, paste0(output_path,filename,'_processed.bedpe'), row.names = F, col.names = T, sep = "\t", quote = F)
  
  return(tp53_added)
}