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

run_GaTSV(file_path,sample,n_cores=32,genome='hg19',output_path = '/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out')
process_file(file_path,sample,n_cores=16,genome="hg19",output_path='/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/Example')
vcf_path_ex="/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/example.sv.vcf"
vcf_path="/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs/ACTN01020001T.purple.sv.vcf.gz"
vcf_dt <- vcf_to_dt_Gridss("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs/ACTN01020001T.purple.sv.vcf.gz", "Test")
vcf_dt <- vcf_to_dt_Gridss("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs/ACTN01020001T.purple.sv.vcf.gz", "Test")
Example_dt <- vcf_to_dt(vcf_path_ex, "Example")
sample_name="Testo"

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

get_metadata <- function(vcf_dt, sample_name){
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
  uids <- lapply(1:nrow(vcf_dt), function(i){
    if (grepl("gridss",vcf_dt$ID[i],fixed=TRUE)){
      uid <- sub(".*_", "", vcf_dt$ID[i]) # Remove everything up to and including the underscore
      uid <- substr(uid, 1, nchar(uid) - 1) # Remove the last character
      return(uid)
    }
    else{
      return(NA)
    }
  })
  # get evidence - TSI ?
  evdnc <- lapply(1:nrow(vcf_dt), function(i){
    #AS <- as.numeric(sub(".*AS=([^;]*);.*", "\\1", vcf_dt$INFO[i]))
    AS <- ifelse(
      startsWith(vcf_dt$INFO[i], "AS="),  # Check if the string starts with "AS="
      as.numeric(sub("^AS=([^;]*);.*$", "\\1", vcf_dt$INFO[i])),  # Extract substring after "AS=" until the first ";"
      NA  # Set to NA if the string does not start with "AS="
    )
    RP <- as.numeric(sub(".*;RP=([^;]*);.*", "\\1", vcf_dt$INFO[i]))
    SR <- as.numeric(sub(".*;SR=([^;]*);.*", "\\1", vcf_dt$INFO[i]))
    if(is.na(SR)|is.na(RP)|is.na(AS)){
      return(NA)
    }
    if(AS>0){
      if(RP>0){
        return("ASDIS")
      }
      else{
        if(SR>0){
          return("TSI_L")
        }
        else{
          return("ASSMB")
        }
      }  
    }
    else{
      if (RP>0){
        return("DSCRD")
      }
      else{
        if (SR>0){
          return("TSI_G")
        }
        else{
          return(NA)
        }
      }
    }
  })
  # get inserted sequence
  ins <- lapply(1:nrow(vcf_dt), function(i){
    if (grepl("\\.", vcf_dt$ALT[i])){
      return(NA)
    }
    before <- sub("\\[.*|\\].*", "", vcf_dt$ALT[i])
    revd <- paste(rev(strsplit(vcf_dt$ALT[i], NULL)[[1]]), collapse = "")
    after_rev <- sub("\\[.*|\\].*", "", revd)
    after <- paste(rev(strsplit(after_rev, NULL)[[1]]), collapse = "")
    if (nchar(before)>1){
      return(before)
    }
    else{
      if(nchar(after)>1){
        return(after)
      }
      else{
        return(NA)
      }
    }
  })
  if ("INFO" %in% colnames(vcf_dt) ) {
    #vcf_dt[, SPAN := as.numeric(gsub(".*?SPAN=([-0-9]+).*","\\1",INFO))]
    vcf_dt[, SPAN := as.numeric(spans)]
    vcf_dt$sample = sample_name
    #vcf_dt[, uid := gsub("([0-9]+):(1|2)", "\\1", ID)]
    vcf_dt[, uid := as.numeric(uids)]
    #vcf_dt[, EVDNC := gsub(".*?EVDNC=([A-Z]+).*", "\\1", INFO)]
    vcf_dt[, EVDNC := as.character(evdnc)]
    #vcf_dt[, MAPQ := as.integer(gsub(".*?;MAPQ=([0-9]+).*", "\\1", INFO))] # No direct equivalent only mean Quality score
    vcf_dt[, MAPQ := as.integer(gsub(".*?;BMQ=([0-9]+).*", "\\1", INFO))]
    vcf_dt[, HOMSEQ := gsub(".*?;HOMSEQ=([A-Z]+).*", "\\1", INFO)]
    vcf_dt[, HOMSEQ := ifelse(grepl(";", HOMSEQ), "", HOMSEQ)] 
    #vcf_dt[, INSERTION := gsub(".*?;INSERTION=([A-Z]+).*", "\\1", INFO)] 
    #vcf_dt[, INSERTION := ifelse(grepl(";", INSERTION), "", INSERTION)]
    vcf_dt[, INSERTION := as.character(ins)]
    # vcf_dt[, NDISC := as.numeric(gsub(".*?NDISC=([0-9]+).*", "\\1", INFO))]
    vcf_dt[, NDISC := NA] # if it does not run, redo
    #vcf_dt[, SVMETHOD := substr(INFO,regexpr("SVMETHOD=",INFO)+nchar("SVMETHOD="),regexpr(";NDISC",INFO)-1)]
    vcf_dt[, SVMETHOD := NA] # if it does not run, redo
  }
  
  # More extraction regexpr stuff...
  if ("TUMOR" %in% colnames(vcf_dt)) {
    #vcf_dt[, TUMALT :=  as.integer(strsplit(TUMOR, ":")[[1]][2]) , by=uid] # might not have that
    vcf_dt[, TUMALT := NA]
    #vcf_dt[, TUMCOV :=  as.integer(strsplit(TUMOR, ":")[[1]][3]) , by=uid] # might not have that
    vcf_dt[, TUMCOV := NA]
    #vcf_dt[, TUMLOD :=  as.numeric(strsplit(TUMOR, ":")[[1]][9]) , by=uid] # does not have that
    vcf_dt[, TUMLOD := NA]
  }
  if ("NORMAL" %in% colnames(vcf_dt)) {
    #vcf_dt[, NORMCOV :=  as.integer(strsplit(NORMAL, ":")[[1]][3]) , by=uid] # might not have that
    vcf_dt[, NORMCOV := NA]
    #vcf_dt[, NORMALT :=  as.integer(strsplit(NORMAL, ":")[[1]][2]) , by=uid] # might not have that
    vcf_dt[, NORMALT := NA]
    #vcf_dt[, NORMLOD :=  as.numeric(strsplit(NORMAL, ":")[[1]][9]) , by=uid] # does not have that
    vcf_dt[, NORMLOD := NA]
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

build_bedpe_with_metadata <- function(merged_dt) {
  cat("Building bedpe...\n")
  
  ### get mate indexes
  #merged_dt[, mates_idx := unlist(strsplit(ID, ":"))[1], by = "uid"]
  merged_dt[, mates_idx := merged_dt$uid]
  #merged_dt[, which_mate := unlist(strsplit(ID, ":"))[2], by = "uid"]
  merged_dt[, which_mate := seq_len(.N), by = uid]
  temp_bedpe <- NULL
  removed_bnd <- NULL
  for(i in 1:length(unique(merged_dt$mates_idx))){
    foo <- merged_dt[mates_idx == unique(merged_dt$mates_idx)[i]]
    
    if(!(nrow(foo)== 2)) {
      mes <- paste0("Breakpoint ",  unique(merged_dt$mates_idx)[i], " has incorrect number of mates for ", foo$sample, " It has been removed.")
      # excludes this breakpoint from 
      continue = FALSE
      removed_bnd <- rbind(removed_bnd, foo)
      warning(mes[1])
    }
    else {
      continue = TRUE
    }
    if(continue) {
      #### build bedpe
      foo[,split_ID := c(1:length(foo$seqnames))] 
      #foo[, which_mate := unlist(strsplit(ID, ":"))[2], by = "split_ID"] 
      foo1 <- foo[which_mate == 1]
      foo2 <- foo[which_mate == 2]
      
      bedpe_base <- as.data.frame(cbind(foo1$seqnames, foo1$start, foo1$end,
                                        foo2$seqnames, foo2$start, foo2$end))
      colnames(bedpe_base) <- c("chrom1", "start1", "end1", "chrom2","start2","end2")
      
      if(!(foo1$sid == foo2$sid)){
        stop("Multiple samples are being processed, one at a time please...")
      }
      
      bedpe_base <- cbind(bedpe_base, paste0(foo$sample[1],"_", foo$mates_idx[1]))
      colnames(bedpe_base)[7] <- "name"
      
      bedpe_base <- cbind(bedpe_base, foo$QUAL[1])
      colnames(bedpe_base)[8] <- "score"
      
      bedpe_base <- cbind(bedpe_base, foo1$strand[1])
      bedpe_base <- cbind(bedpe_base, foo2$strand[1])
      colnames(bedpe_base)[9:10] <- c("strand1", "strand2")
      
      refs_alts <- as.data.frame(cbind(foo1$REF[1], foo1$ALT[1],
                                       foo2$REF[1], foo2$ALT[1]))
      colnames(refs_alts) <- c("REF_1","ALT_1","REF_2","ALT_2")
      bedpe_base <- cbind(bedpe_base, refs_alts)
      bedpe_base <- cbind(bedpe_base, foo1[,c("SPAN", "HOMSEQ","INSERTION","NDISC","FILTER","sample", "EVDNC","TUMALT", 'GENO','TUMOR')])
      mapqs <- as.data.frame(cbind(foo1$MAPQ[1], foo2$MAPQ[1]))
      colnames(mapqs) <- c("MAPQ_1","MAPQ_2")
      bedpe_base <- cbind(bedpe_base, mapqs)
      temp_bedpe <- rbind(temp_bedpe, bedpe_base)
    }
  }
  
  bedpe <- as.data.table(temp_bedpe)
  return(bedpe)
}

filter_gridss <- function(lof_pth, sample) {
  
  vcf_dt <- vcf_to_dt_Gridss(lof_pth, sample)
  dt_with_mets <- get_metadata(vcf_dt, sample)
  vcf_bedpe <- build_bedpe_with_metadata(dt_with_mets)
  
  # as substitute for coverage: use number of split reads + number of spanning reads + number of soft clips
  #vcf_bedpe[, NALT_SR := unlist(strsplit(TUMOR, ":"))[3], by = 'TUMOR']
  cov_est <- lapply(1:nrow(vcf_bedpe), function(i){
    SR=as.numeric(unlist(strsplit(vcf_bedpe$TUMOR[i], ":"))[27])
    REF=as.numeric(unlist(strsplit(vcf_bedpe$TUMOR[i], ":"))[23])
    BSC=as.numeric(unlist(strsplit(vcf_bedpe$TUMOR[i], ":"))[13])
    cov_est=SR+REF+BSC
  })
  vcf_bedpe[, NALT_SR := as.numeric(cov_est)]
  #vcf_bedpe[, NALT := unlist(strsplit(TUMOR, ":"))[1], by = 'TUMOR']
  
  tmp <- vcf_bedpe[MAPQ_1 == 60 | MAPQ_2 == 60] #takes only reads that have a mapq score that are 60 or higher since MAPQ scores are capped at 60
  tmp <- tmp[!(EVDNC == 'DSCRD')] #removes reads that are discordant
  tmp <- tmp[NALT_SR > 1] #the number of reads covering the site/depth of coverage must be greater than 1
  tmp <- tmp[SPAN > 49 | SPAN== -1] #SPAN = -1 refers to a translocation. Span shorter than 50bp are considered simple indels
  return(tmp)
  #}
}
andre3000 <- process_file_gridss("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs/ACTN01020001T.purple.sv.vcf.gz","Testo",4,"hg19","/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/Testo/")
file="/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs/ACTN01020001T.purple.sv.vcf.gz"
sample="Testo"
n_cores=32
genome="hg19"
output_path="/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/Testo/"

process_file_gridss <- function(file,sample, n_cores,genome,output_path='./') {
  
  filt_file <- suppressWarnings(filter_gridss(file, sample))
  
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

andre3000 <- process_file_gridss("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs/ACTN01020001T.purple.sv.vcf.gz","Testo",4,"hg19","/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/Testo/")


file_path="/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs/ACTN01020001T.purple.sv.vcf.gz"
sample="Testo"
n_cores=32
genome="hg19"
output_path="/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/Testo/"

run_GaTSV_gridss <- function(file_path,sample,n_cores=1,genome='hg19',output_path = './'){
  cat(paste0('Reference genome: ',genome,'\n'))
  cat(paste0('Writing outputs to: ', output_path,'\n'))
  cat(paste0('Using ',n_cores,' core(s) for parallel processing \n'))
  tmp <- process_file_gridss(file = file_path, sample=sample, n_cores =n_cores,genome=genome,output_path=output_path)
  cutoff_prob <-  0.2684 #optimal tpr+ppv cutoff
  
  test <- add_last_feat(tmp)
  test <- test[SPAN>=1e3|SPAN==-1,]
  test[, log_SPAN := log(SPAN, base = 10)] 
  test$log_SPAN[which(is.na(test$log_SPAN))] <- max(subset(test$log_SPAN, !is.na(test$log_SPAN)))
  test$log_SPAN[which(is.infinite(test$log_SPAN))] <- max(subset(test$log_SPAN, !is.infinite(test$log_SPAN)))
  test_sub_prelog <- test[, .SD,.SDcols = features_tolog]
  test_sub_log <- data.table(apply(test_sub_prelog,2,FUN = log_feat))
  colnames(test_sub_log) <- paste0('log_', colnames(test_sub_log))
  
  test <- cbind(test, test_sub_log)
  test_df <- as.data.frame(test[,.SD,.SDcols=features_toscale]) #select the features to scale
  gtest_scaled <- data.table()
  for (i in colnames(test_df)){
    row_l <- scaling_mat[which(scaling_mat$feature== i),]
    feature_col <- test_df[grepl(i,colnames(test_df))]
    scaled_feature <- (feature_col -(row_l$mean))/row_l$sd
    gtest_scaled <- cbind(gtest_scaled,scaled_feature)
  }
  cat("Performing classification...\n")
  y_pred_radial <- predict(GaTSV, newdata = gtest_scaled, decision.values = T, probability = T)
  probabilities_radial <-data.table(attr(y_pred_radial, 'probabilities'))
  setcolorder(probabilities_radial, c('0', '1'))
  
  probabilities_radial[,"pred_class"] <- lapply(1:length(probabilities_radial$`1`),function(i){
    return (ifelse(probabilities_radial$`1`[i]>= cutoff_prob,'SOMATIC','GERMLINE'))})
  test <- cbind(test,as.character(probabilities_radial$pred_class))
  colnames(test)[ncol(test)] <- 'predicted_class'
  filename <- test$sample[1]
  write.table(test, paste0(output_path,filename,'_classified.bedpe'), row.names = F, col.names = T, sep = "\t", quote = F)
  cat('done.')
  #return(test)
}

run_GaTSV_gridss(file,sample,n_cores, genome, output_path)




#running the classifier on example data
metadata <- fread("/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/example_metadata.txt") #metadata file that contains the sample_ids (same as 'sample' input) and associated tp53_mutation_status
file_path <- "/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/example.sv.vcf" #replace with desired vcf path
sample <- "example"

#make metadatafile
# Define the directory and output file path
input_dir <- "/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs"
output_file <- "/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/samples_tp53_status.txt"

# List all files in the directory
file_list <- list.files(input_dir)

# Extract substrings until the first '.' and remove duplicates
sample_list <- unique(sub("\\..*", "", file_list))

# Create a data frame with 'sample' and 'tp53_mutation_status' columns
result_df <- data.frame(
  sample = sample_list,
  tp53_mutation_status = NA  # Fill with NA
)

# Write the data frame to a text file in a tab-delimited format for fread compatibility
write.table(
  result_df, 
  file = output_file, 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)

cat("File written to:", output_file, "\n")

run_GaTSV(file_path,sample,n_cores=16,genome='hg19',output_path = '/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out')
process_file(file_path,sample,n_cores=16,genome="hg19",output_path='/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/Example')
vcf_path_ex="/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/example.sv.vcf"
vcf_path="/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs/ACTN01020001T.purple.sv.vcf.gz"
vcf_dt <- vcf_to_dt_Gridss("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs/ACTN01020001T.purple.sv.vcf.gz", "Test")
vcf_dt <- vcf_to_dt_Gridss("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/WGS/HMF_GRIDSS_vcfs/ACTN01020001T.purple.sv.vcf.gz", "Test")
Example_dt <- vcf_to_dt(vcf_path_ex, "Example")
sample_name="Testo"

file_path="/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/data/example.sv.vcf"
sample="Example"
n_cores=32
genome="hg19"
output_path='/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/Example'

run_GaTSV <- function(file_path,sample,n_cores=1,genome='hg19',output_path = './'){
  cat(paste0('Reference genome: ',genome,'\n'))
  cat(paste0('Writing outputs to: ', output_path,'\n'))
  cat(paste0('Using ',n_cores,' core(s) for parallel processing \n'))
  tmp <- process_file(file = file_path, sample=sample, n_cores =n_cores,genome=genome,output_path=output_path)
  cutoff_prob <-  0.2684 #optimal tpr+ppv cutoff
  
  test <- add_last_feat(tmp)
  test <- test[SPAN>=1e3|SPAN==-1,]
  test[, log_SPAN := log(SPAN, base = 10)] 
  test$log_SPAN[which(is.na(test$log_SPAN))] <- max(subset(test$log_SPAN, !is.na(test$log_SPAN)))
  test$log_SPAN[which(is.infinite(test$log_SPAN))] <- max(subset(test$log_SPAN, !is.infinite(test$log_SPAN)))
  test_sub_prelog <- test[, .SD,.SDcols = features_tolog]
  test_sub_log <- data.table(apply(test_sub_prelog,2,FUN = log_feat))
  colnames(test_sub_log) <- paste0('log_', colnames(test_sub_log))
  
  test <- cbind(test, test_sub_log)
  test_df <- as.data.frame(test[,.SD,.SDcols=features_toscale]) #select the features to scale
  test_scaled <- data.table()
  for (i in colnames(test_df)){
    row_l <- scaling_mat[which(scaling_mat$feature== i),]
    feature_col <- test_df[grepl(i,colnames(test_df))]
    scaled_feature <- (feature_col -(row_l$mean))/row_l$sd
    test_scaled <- cbind(test_scaled,scaled_feature)
  }
  cat("Performing classification...\n")
  y_pred_radial <- predict(GaTSV, newdata = test_scaled, decision.values = T, probability = T)
  probabilities_radial <-data.table(attr(y_pred_radial, 'probabilities'))
  setcolorder(probabilities_radial, c('0', '1'))
  
  probabilities_radial[,"pred_class"] <- lapply(1:length(probabilities_radial$`1`),function(i){
    return (ifelse(probabilities_radial$`1`[i]>= cutoff_prob,'SOMATIC','GERMLINE'))})
  test <- cbind(test,as.character(probabilities_radial$pred_class))
  colnames(test)[ncol(test)] <- 'predicted_class'
  filename <- test$sample[1]
  write.table(test, paste0(output_path,filename,'_classified.bedpe'), row.names = F, col.names = T, sep = "\t", quote = F)
  cat('done.')
  #return(test)
}
