library(BSgenome.Hsapiens.UCSC.hg19)

suppressPackageStartupMessages(require(BiocGenerics))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(caTools))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(e1071))
suppressPackageStartupMessages(require(GenomeInfoDb))
suppressPackageStartupMessages(require(GenomicRanges))
suppressPackageStartupMessages(require(gUtils))
suppressPackageStartupMessages(require(IRanges))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(rlang))
suppressPackageStartupMessages(require(ROCR))
suppressPackageStartupMessages(library(rstudioapi))
suppressPackageStartupMessages(require(S4Vectors))
suppressPackageStartupMessages(require(stats4))
suppressPackageStartupMessages(require(stringr))
suppressPackageStartupMessages(require(here))

################################################################################
####### LOAD FUNCTIONS #########################################################

# THIS ASSUMES THAT THE CHROMOSOMES YOU HAVE ARE 1-22, XY
# DIFFERENT ANNOTATIONS ARE NECESSARY IF YOU HAVE 1-24
# THIS ALSO ASSUMES YOU DONT HAVE chr IN THE seqnames ANT ALT
# IF YOU DO, YOU CAN SKIP SOME OF THESE STEPS

filter_adequate_reads <- function(vcf_dt) {

  check_tumor_valid <- function(i, cur_dt) {
    cur_row <- cur_dt[i,]
    format <- cur_row$FORMAT
    evidence <- cur_row$TUMOR

    if(format == "PR:SR") {
      supporting_reads <- unlist(strsplit(evidence, "[:|,]"))
      pr_alt <- as.numeric(supporting_reads[2])
      sr_alt <- as.numeric(supporting_reads[4])
      if(pr_alt > 0 & sr_alt > 0) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      return(FALSE)
    }
  }

  check_normal_valid <- function(i, cur_dt) {
    cur_row <- cur_dt[i,]

    evidence <- cur_row$NORMAL
    evidence_split <- unlist(strsplit(evidence, "[:|,]"))
    pr_alt <- as.numeric(evidence_split[length(evidence_split) - 2])
    sr_alt <- as.numeric(evidence_split[length(evidence_split)])

    if(pr_alt > 0 & sr_alt > 0) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }

  if (ncol(vcf_dt) == 12) {
    vcf_dt$read_supported <- lapply(1:nrow(vcf_dt), check_tumor_valid, vcf_dt)
    vcf_dt <- vcf_dt[read_supported == TRUE]


  } else if (ncol(vcf_dt) == 11) {
    vcf_dt$read_supported <- lapply(1:nrow(vcf_dt), check_normal_valid, vcf_dt)
    vcf_dt <- vcf_dt[read_supported == TRUE]
  } else {
    stop("Unexpected number of columns")
  }

  return(vcf_dt)
}

preprocess_manta_vcfs <- function(vcf_dt) {
  # FOR SOMATIC
  colnames12 <- c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO", "FORMAT","NORMAL","TUMOR","Sample")
  # FOR DIPLOID/GERMLINE
  colnames11 <- c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO", "FORMAT","NORMAL","Sample")

  if (ncol(vcf_dt) == 12) {
    colnames(vcf_dt) <- colnames12
  } else if (ncol(vcf_dt) == 11) {
    colnames(vcf_dt) <- colnames11
  } else {
    stop("Unexpected number of columns")
  }

  vcf_dt <- filter_adequate_reads(vcf_dt)

  sv.calls <- subset(vcf_dt, select = c("seqnames","start","ID","REF","ALT","QUAL","FILTER", "INFO","Sample"))

  canonical_chr <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                     "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                     "21", "22", "X", "Y")

  # Modify seqnames column
  sv.calls <- sv.calls[seqnames %in% canonical_chr]
  sv.calls$seqnames <- paste0("chr", sv.calls$seqnames)

  # Check which of ALT has non-canonical chromosomes
  noncanon_rows <- function(cell) {
    if(grepl("\\[|\\]", cell) & grepl(":", cell)) {
      pre_colon <- unlist(strsplit(cell, ":", fixed = TRUE))[1]
      part2 <- unlist(strsplit(pre_colon, "\\[|\\]", perl = TRUE))[2]

      mod_string <- part2 %in% canonical_chr

      return(mod_string)
    } else {
      return(TRUE)
    }
  }
  to_keep <- unlist(lapply(sv.calls$ALT, noncanon_rows))
  sv.calls <- sv.calls[to_keep]

  add_chr <- function(cell) {
    if(grepl("\\[|\\]", cell) & grepl(":", cell)) {
      pre_colon <- unlist(strsplit(cell, ":", fixed = TRUE))[1]
      post_colon <- paste0(":", unlist(strsplit(cell, ":", fixed = TRUE))[2])

      part1 <- unlist(strsplit(pre_colon, "\\[|\\]", perl = TRUE))[1]
      part1 <- paste0(part1, substr(cell, nchar(part1) + 1, nchar(part1) + 1))
      part2 <- paste0("chr", unlist(strsplit(pre_colon, "\\[|\\]", perl = TRUE))[2])

      mod_string <- paste0(part1, part2, post_colon)

      return(mod_string)
    } else {
      return(cell)
    }
  }
  sv.calls$ALT <- lapply(sv.calls$ALT, add_chr)

  sv.calls <- sv.calls[sv.calls$FILTER == "PASS", ]

  return(sv.calls)
}


reformat_manta_vcfs <- function(sv_calls){
  # processing manta deletion and duplication calls
  # use ID column instead of alt column bc some deletion calls, the alt column doesn't have the mantaDEL tag
  # manta output contains deletions > 50bp that are not annotated with <DEL> in the ALT column
  # so ID column is best
  # sv.calls[grepl('DEL',ID) & ALT != '<DEL>',]
  dup_del_sv <- sv_calls[grepl("TANDEM|DEL", ID)]
  dup_del_sv_2nd <- copy(dup_del_sv)
  if(nrow(dup_del_sv) > 0){
    dup_del_sv[, c("ALT", "ID") := list(ifelse(
      grepl("MantaDEL", ID),
      # for the first breakend, for deletions, the lower genomic breakpoint is the + end
      # therefore they should be given t[p[ notation
      paste0(REF, "[", seqnames, ":", gsub("END=", "", str_extract(INFO, "END=[0-9]*")), "["),
      paste0("]", seqnames, ":", gsub("END=", "", str_extract(INFO, "END=[0-9]*")), "]", REF)
    ), paste0(ID, ":1"))]
    # extract the second breakend genomic coordinate for the same set of calls
    dup_del_sv_2nd[, REF := as.character(getSeq(
      Hsapiens, seqnames, as.numeric(gsub("END=", "", str_extract(INFO, "END=[0-9]*"))),
      as.numeric(gsub("END=", "", str_extract(INFO, "END=[0-9]*")))
    ))][, c("ALT", "ID") := list(ifelse(
      grepl("MantaDEL", ID),
      # for the second breakend, for deletions, the upper genomic breakpoint is the - end
      # therefore they should be given t[p[ notation
      paste0("]", seqnames, ":", start, "]", REF),
      paste0(REF, "[", seqnames, ":", start, "[")
    ), paste0(ID, ":2"))]
    dup_del_sv_2nd$start <- as.numeric(gsub("END=", "", str_extract(dup_del_sv$INFO, "END=[0-9]*")))
  }

  # Processing manta inversion calls (this is mainly for the 1KG data where `convertInversion.py` has been run)
  inv_calls <- sv_calls[grepl("MantaINV", ID)]
  inv_call_2nd <- copy(inv_calls)
  if(nrow(inv_calls) > 0){
    # for inversions, the calls with INV5 in info are -/- and INV3 are +/+
    inv_calls[, c("ALT", "ID") := list(
      ifelse(
        grepl("INV5", INFO),
        paste0("[", seqnames, ":", gsub('END=', '', str_extract(INFO, "END=[0-9]*")), "[", REF),
        paste0(REF, "]", seqnames, ":", gsub('END=', '', str_extract(INFO, "END=[0-9]*")), "]")
      ), paste0(ID, ":1")
    )]
    inv_call_2nd[, 'REF' := as.character(getSeq(
      Hsapiens, seqnames, as.numeric(gsub('END=', '', str_extract(INFO, 'END=[0-9]*'))),
      as.numeric(gsub('END=', '', str_extract(INFO, 'END=[0-9]*')))
    ))][, c("ALT", "ID",'INFO') := list(
      ifelse(
        grepl("INV5", INFO),
        paste0("[", seqnames, ":", start, "[", REF),
        paste0(REF, "]", seqnames, ":", start, "]")
      ), paste0(ID, ":2"),
      # swap out the SVINSSEQ string in inv_call_2nd 'SVINSSEQ=[AGCT]*' with the revercomplement of the string
      ifelse(
        grepl("SVINSSEQ", INFO),
        gsub('SVINSSEQ=[AGCT]*', paste0(
          'SVINSSEQ=', reverseComplement(DNAString(
            gsub('SVINSSEQ=', '', str_extract(INFO, 'SVINSSEQ=[AGCT]*'))
          ))
        ), INFO), INFO
      )
    ),.(ID, Sample)]
    inv_call_2nd$start <- as.numeric(gsub('END=', '', str_extract(inv_calls$INFO, 'END=[0-9]*')))
  }

  rbind(sv_calls[!grepl("TANDEM|DEL", ID)], dup_del_sv, dup_del_sv_2nd, inv_calls, inv_call_2nd)
}

# filter for major chromosomes
# filter for minimal qual score (svaba), FILTER == 'PASS' (manta)
# expand and reformat 1 row SVs (mantaDEL/DUP/INV) into 2 rows per SV
# make sure that every SV have both ends passing the above criteria
filter_and_reformat_sv_calls <- function(sv.calls, caller){
  # filter for only the major chromosomes
  sv.calls <- sv.calls[seqnames %in% c(c(1:22, "X", "Y"), paste0("chr", c(1:22, "X", "Y"))), ]

  # filter for minimal qual score (svaba), FILTER == 'PASS' (manta)
  # expand and reformat 1 row SVs (mantaDEL/DUP/INV) into 2 rows per SV
  # if (caller %in% c("consensus")) {
  # sv.calls$callers <- gsub("CALLER=", "", str_extract(sv.calls$INFO, "CALLER=[SD]*"))
  # print(table(sv.calls$callers))
  # print(paste0("taking out DS or SD calls which makes up ", round(sum(sv.calls$callers == "DS" | sv.calls$callers == "SD") / nrow(sv.calls), 3), " of the data"))
  # sv.calls <- sv.calls[grepl("CALLER=SD|CALLER=DS", INFO)]
  # } else
  if (caller %in% c("manta")) {
    sv.calls <- sv.calls[FILTER == 'PASS']
    sv.calls <- reformat_manta_vcfs(sv.calls)
  } else {
    # with svaba we use QUAL score cut off
    sv.calls <- sv.calls[as.numeric(QUAL) > 15]
  }

  # making sure every breakend will have its pair
  sv.calls[, breakend_ID := paste(Sample, seqnames, ID, sep = "__")]
  sv.calls[, both_end_pass_filter := ifelse(.N == 2, TRUE, FALSE), , by = .(SV_pair_ID = gsub(":[0-2]$", "", ID), Sample)]
  sv.calls[, breakend_num_per_sv_ID := .N, , by = .(SV_pair_ID = gsub(":[0-2]$", "", ID), Sample)]
  print(paste0('Total number of breakends: ',nrow(sv.calls)))
  print(paste0("Both end pass filter: ", nrow(sv.calls[both_end_pass_filter == TRUE])))
  sv.calls <- sv.calls[both_end_pass_filter == TRUE]
  return(sv.calls)
}


# Extracting info from INFO column =====================================================================================
extract_info <- function(sv.calls){
  # extract insertion sequences from the info column
  sv.calls$ins_seq <- gsub("INSERTION=|FORSEQ=|SVINSSEQ=", "", str_extract(sv.calls$INFO, "INSERTION=[ACGT]*|FORSEQ=[ACGT]*|SVINSSEQ=[ACGT]*"))
  sv.calls$ins_len <- nchar(sv.calls$ins_seq)
  # extracting MH
  sv.calls$mh_seq <- gsub("HOMSEQ=", "", str_extract(sv.calls$INFO, "HOMSEQ=[ACGT]*"))
  sv.calls$mh_len <- nchar(sv.calls$mh_seq)
  # extracting evidence type
  sv.calls$evdnc <- gsub("EVDNC=", "", str_extract(sv.calls$INFO, "EVDNC=[A-Z\\_]*"))
  # extracting MAPQ
  sv.calls$mapq <- as.numeric(gsub("MAPQ=", "", str_extract(sv.calls$INFO, "MAPQ=[0-9]*")))
  sv.calls$disc_mapq <- as.numeric(gsub("DISC_MAPQ=", "", str_extract(sv.calls$INFO, "DISC_MAPQ=[0-9]*")))
  # extracting +/- breakpoint
  sv.calls$cnt_type <- unlist(lapply(grepl("^[AGCTN]", sv.calls$ALT), function(x) {
    ifelse(x, "+", "-")
  }))
  return(sv.calls)
}


collapse_VCF_to_bedpe <- function(vcf_data){
  # minimal columns needed -- Sample, seqnames, start, cnt_type, ID
  vcf_data[, c("SV_ID", "breakend_num") :=
             list(
               paste0(Sample, "__", gsub(":[0-2]$", "", ID)),
               gsub(".*:", "", ID)
             )]
  print(paste0('There are ',nrow(vcf_data[,.N,SV_ID][N == 2]),' SVs with two breakend.'))
  # there are some SVs without two breakends, we will discard them
  vcf_data <- vcf_data[SV_ID %in% vcf_data[,.N,SV_ID][N == 2]$SV_ID]
  # per sv pair, we sort the breakend numbers, sort them, and give them a unified order
  vcf_data[order(breakend_num), breakend_order := c("breakend1", "breakend2"), SV_ID]
  # SV_config_combo
  vcf_data[,'SV_config_combo' := paste0(cnt_type,collapse = ''),SV_ID]
  # dcast into the short format with breakend 1 and 2 identified by SV_ID

  print(vcf_data)

  bedpe_data <- dcast(
    # melt into long format
    melt(
      vcf_data,
      id.vars = c(
        "SV_ID", "Sample", "breakend_order", "SV_config_combo"
      )
    )[, variable := paste0(variable, "_", breakend_order)][, breakend_order := NULL],
    SV_ID + Sample + SV_config_combo ~ variable,
    value.var = "value"
  )
  # print(bedpe_data[
  #     ,.N,.(
  #         svtype = gsub('.*__|:.*','',SV_ID),
  #         SV_config_combo
  #     )
  # ][order(svtype)])
  # there are some bugs from old code that we will ignore here
  if(nrow(bedpe_data[grepl('DUP|DEL',SV_ID) & SV_config_combo %in% c('--','++')])/nrow(bedpe_data) == 0){
    # nice, nothing is wrong
    print('Connection Type Validated by SV Types')
  }else if(nrow(bedpe_data[grepl('DUP|DEL',SV_ID) & SV_config_combo %in% c('--','++')])/nrow(bedpe_data) < 0.0002){
    print(paste0(
      "There are some bugs from old code that we will ignore here. ",
      "The number of deletions and duplications with SV_config_combo == '--' or '++' is ",
      nrow(bedpe_data[grepl('DUP|DEL',SV_ID) & SV_config_combo %in% c('--','++')])/nrow(bedpe_data)*100,
      "% of the total SVs."
    ))
    bedpe_data <- bedpe_data[!(grepl('DUP|DEL',SV_ID) & SV_config_combo %in% c('--','++'))]
  }else{
    stop('too many bugs from old code')
  }
  bedpe_data[,svclass := ifelse(
    as.character(seqnames_breakend1) == as.character(seqnames_breakend2),
    ifelse(
      as.character(cnt_type_breakend1) == as.character(cnt_type_breakend2),
      ifelse(as.character(cnt_type_breakend1) == '+', "h2hINV", "t2tINV"),
      ifelse(as.character(cnt_type_breakend1) == "+",  "DEL", "DUP")
    ),
    "TRA"
  )]
  bedpe_data$svclass <- factor(bedpe_data$svclass, levels = c("DEL", "DUP", "h2hINV", "t2tINV", "TRA"))
  bedpe_data
}


make_bedpe_clean <- function(sv_bedpe) {
  # Input the bedpe output from collapse_VCF_to_bedpe and outputs a bedpe that matches the format of what we will need
  # First convert the SVs to the annotations we use

  sv_bedpe_df <- as.data.frame(sv_bedpe)

  sv_bedpe_df$chrom1 <- as.character(ifelse(grepl('chr', sv_bedpe$seqnames_breakend1), sv_bedpe$seqnames_breakend1, paste0("chr", sv_bedpe$seqnames_breakend1)))
  sv_bedpe_df$start1 <- as.numeric(sv_bedpe$start_breakend1)
  sv_bedpe_df$end1 <- as.numeric(sv_bedpe$start_breakend1)
  sv_bedpe_df$chrom2 <- as.character(ifelse(grepl('chr', sv_bedpe$seqnames_breakend2), sv_bedpe$seqnames_breakend2, paste0("chr", sv_bedpe$seqnames_breakend2)))
  sv_bedpe_df$start2 <- as.numeric(sv_bedpe$start_breakend2)
  sv_bedpe_df$end2 <- as.numeric(sv_bedpe$start_breakend2)
  sv_bedpe_df$name <- as.character(sv_bedpe$SV_ID)
  sv_bedpe_df$score <- NA
  sv_bedpe_df$strand1 <- as.character(sv_bedpe$cnt_type_breakend1)
  sv_bedpe_df$strand2 <- as.character(sv_bedpe$cnt_type_breakend2)
  sv_bedpe_df$REF_1 <- as.character(sv_bedpe$REF_breakend1)
  sv_bedpe_df$ALT_1 <- as.character(sv_bedpe$ALT_breakend1)
  sv_bedpe_df$REF_2 <- as.character(sv_bedpe$REF_breakend2)
  sv_bedpe_df$ALT_2 <- as.character(sv_bedpe$ALT_breakend2)


  sv_bedpe_df$SPAN <- ifelse(sv_bedpe_df$chrom1 != sv_bedpe_df$chrom2, -1, abs(sv_bedpe_df$start1 - sv_bedpe_df$start2))
  sv_bedpe_df$HOMSEQ <- sv_bedpe$mh_seq_breakend1
  sv_bedpe_df$INSERTION <- sv_bedpe$ins_seq_breakend1
  sv_bedpe_df$sample <- sv_bedpe$Sample
  sv_bedpe_df$SVTYPE_MANTA <- sv_bedpe$svclass
  sv_bedpe_df$INFO1 <- sv_bedpe$INFO_breakend1
  sv_bedpe_df$INFO2 <- sv_bedpe$INFO_breakend2

  sv_bedpe_df <- subset(sv_bedpe_df, select = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2",
                                                "REF_1", "ALT_1", "REF_2", "ALT_2", "SPAN", "HOMSEQ", "INSERTION", "sample", "SVTYPE_MANTA",
                                                "INFO1", "INFO2"))
  sv_bedpe_df$sample <- gsub("(_somaticSV\\.vcf|_diploidSV\\.vcf)$", "", sv_bedpe_df$sample)

  return(sv_bedpe_df)
}

################################################################################
## GATSV ANNOTATION SCRIPTS ##
source("/Volumes/xchip_beroukhimlab/siyun/germline_classifier/scripts/post_preprint/20241004_reviewer_comments/gatsv_bedpe_annotation_scripts.R")

################################################################################
################################################################################

# read in manta VCFs
setwd("/Volumes/xchip_beroukhimlab/siyun/germline_classifier/scripts/post_preprint/20241004_reviewer_comments/phgg_manta/vcfs/")

somatic_manta_vcf_paths <- list.files("./somatic", full.names = TRUE)
germline_manta_vcf_paths <- list.files("./diploid", full.names = TRUE)

somatic_manta_vcfs <- do.call(rbind, lapply(somatic_manta_vcf_paths, function(x){
  tmp <- data.table::fread(cmd=paste("grep -v '^#'", x),sep='\t')
  tmp[, V12:=gsub("./somatic/", "", x)]
  return(tmp)
}))

germline_manta_vcfs <- do.call(rbind, lapply(germline_manta_vcf_paths, function(x){
  tmp <- data.table::fread(cmd=paste("grep -v '^#'", x),sep='\t')
  tmp[, V12:=gsub("./diploid/", "", x)]
  return(tmp)
}))

processing_workflow <- function(sv.calls) {
  # Input vcf as a datatable, returns bedpe file
  sv.calls <- preprocess_manta_vcfs(sv.calls)
  sv.calls <- reformat_manta_vcfs(sv.calls)
  sv.calls <- filter_and_reformat_sv_calls(sv.calls, "manta")
  sv.calls <- extract_info(sv.calls)
  # This will have a warning message, but it's ok
  sv_bedpe <- collapse_VCF_to_bedpe(sv.calls)
  sv_bedpe <- make_bedpe_clean(sv_bedpe)
}

somatic_manta_bedpe <- processing_workflow(somatic_manta_vcfs)
germline_manta_bedpe <- processing_workflow(germline_manta_vcfs)

################  PHGG TP53 CALLS  #############################################
################################################################################
phGG_mutectcalls <- fread("/Volumes/xchip_beroukhimlab/wolu/testing_svaba/scripts/phGG_casestudy/20190714M2maf.maf")
case_ids <- c("00-001R", "00-004R", "00-019P", "00-030P", "10-417-574_S11_T", "10-417-2218_S12_T", "10-417-6804-A17_T")
high_sv_phGGM2 <- subset(phGG_mutectcalls,Tumor_Sample_Barcode %in% case_ids) #we know that all cases have mutation calls so we know that it's either 1 or -1
high_sv_phGGM2_filt <- high_sv_phGGM2[!(which(Variant_Classification %in% c("3'UTR", "5'UTR", "incRNA", "Intron", "Silent", "IGR", "5'Flank", "RNA"))),]
high_sv_phGGM2_filt_tp53 <- subset(high_sv_phGGM2_filt, Hugo_Symbol == 'TP53')
phGG_tp53_status <- as.data.table(case_ids)
colnames(phGG_tp53_status) <- 'sample_name'
phGG_tp53_status[,tp53_mut_status := ifelse(sample_name %in% high_sv_phGGM2_filt_tp53$Tumor_Sample_Barcode,1,-1)]
################################################################################
################################################################################

# set metadata and
metadata <- phGG_tp53_status
colnames(metadata) <- c("sample", "tp53_status")
metadata$sample <- c("00_001R", "00-004R", "00-019P", "00_030P", "10-417-574", "10_417_2218", "10_417_6804")

somatic_manta_bedpe <- merge(somatic_manta_bedpe, metadata)
germline_manta_bedpe <- merge(germline_manta_bedpe, metadata)

# Reorder columns

somatic_manta_bedpe <- subset(somatic_manta_bedpe, select = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2",
                                              "REF_1", "ALT_1", "REF_2", "ALT_2", "SPAN", "HOMSEQ", "INSERTION", "sample", "SVTYPE_MANTA",
                                              "INFO1", "INFO2", "tp53_status"))
germline_manta_bedpe <- subset(germline_manta_bedpe, select = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2",
                                                              "REF_1", "ALT_1", "REF_2", "ALT_2", "SPAN", "HOMSEQ", "INSERTION", "sample", "SVTYPE_MANTA",
                                                              "INFO1", "INFO2", "tp53_status"))

saveRDS(somatic_manta_bedpe, paste0(cur_dir, "20241101_somatic_manta_bedpe_preannotation_READFILT.rds"))
saveRDS(germline_manta_bedpe, paste0(cur_dir, "20241101_germline_manta_bedpe_preannotation_READFILT.rds"))

cur_dir <- "/Volumes/xchip_beroukhimlab/siyun/germline_classifier/scripts/post_preprint/20241004_reviewer_comments/phgg_manta/"
gnomad_hg19 <- readRDS(paste0(cur_dir, "gnomAD.v4.hg19.liftover.rds"))
LINE_dt_hg19 <- readRDS(paste0(cur_dir, "repeatmasker.hg19.LINE.bed"))
SINE_dt_hg19 <- readRDS(paste0(cur_dir, "repeatmasker.hg19.SINE.bed"))
reptimedata_hg19 <- readRDS(paste0(cur_dir, "reptime.hg19.rds"))
hg19_genes <- readRDS(paste0(cur_dir, "gencode.genes.hg19.rds"))
hg19_exons <- readRDS(paste0(cur_dir, "gencode.exons.hg19.rds"))
scaling_mat <- fread(paste0(cur_dir, "scalingmatrix.txt"))

somatic_manta_bedpe_processed <- process_file(somatic_manta_bedpe, n_cores = 2, genome='hg19', filename = '20241029_phgg_somatic_partAnnotated', output_path = cur_dir)
germline_manta_bedpe_processed <- process_file(germline_manta_bedpe, n_cores = 6, genome='hg19', filename = '20241029_phgg_germline_partAnnotated', output_path = cur_dir)

######################### POST MOST PROCESSING #################################
################################################################################

somatic_manta_bedpe_processed <- readRDS(paste0(cur_dir, "20241101_phgg_somatic_partAnnotated_READFILT.rds"))
germline_manta_bedpe_processed <- readRDS(paste0(cur_dir, "20241101_phgg_germline_partAnnotated_READFILT.rds"))

somatic_manta_bedpe_processed %>% group_by(SVTYPE_MANTA, svtype) %>% summarise(counts = n())
germline_manta_bedpe_processed %>% group_by(SVTYPE_MANTA, svtype) %>% summarise(counts = n())

somatic_manta_bedpe_processed <- add_last_feat(somatic_manta_bedpe_processed)
germline_manta_bedpe_processed <- add_last_feat(germline_manta_bedpe_processed)

somatic_manta_test_scaled <- final_processing(somatic_manta_bedpe_processed, 1)
germline_manta_test_scaled <- final_processing(germline_manta_bedpe_processed, 0)
manta_test_scaled <- rbind(somatic_manta_test_scaled, germline_manta_test_scaled)

#saveRDS(manta_test_scaled, paste0(cur_dir, "20241101_phgg_manta_test_scaled.rds"))
saveRDS(manta_test_scaled, paste0(cur_dir, "phgg_manta_test_bedpe_UNFILTERED1k.rds"))

somatic_manta_test_bedpe <- final_processing(somatic_manta_bedpe_processed, 1)
germline_manta_test_bedpe <- final_processing(germline_manta_bedpe_processed, 0)
manta_test_bedpe <- rbind(somatic_manta_test_bedpe, germline_manta_test_bedpe)

saveRDS(manta_test_bedpe, paste0(cur_dir, "phgg_manta_test_bedpe.rds"))
