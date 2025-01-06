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

## TEST SVM ####################################################################
################################################################################
cur_dir <- "/Volumes/xchip_beroukhimlab/siyun/germline_classifier/scripts/post_preprint/20241004_reviewer_comments/phgg_manta/"
GaTSV <- readRDS(paste0(cur_dir, "GaTSV.rda"))

#test_scaled <- readRDS(paste0(cur_dir, "/20241118_phgg_scaled_LOWFILT.rds"))
#test_scaled <- readRDS(paste0(cur_dir, "/20241118_phgg_scaled_MIDFILT.rds"))
test_scaled <- readRDS(paste0(cur_dir, "/20241118_phgg_scaled_HIFILT.rds"))

y_pred_radial <- predict(GaTSV, newdata = test_scaled, decision.values = T, probability = T)

probabilities_radial <-data.table(attr(y_pred_radial, 'probabilities'))
setcolorder(probabilities_radial, c('0', '1'))

cutoff_prob <-  0.2684
probabilities_radial[,"pred_class"] <- lapply(1:length(probabilities_radial$`1`),function(i){
  return (ifelse(probabilities_radial$`1`[i]>= cutoff_prob,'SOMATIC','GERMLINE'))})

#test <- readRDS(paste0(cur_dir, "/20241118_phgg_1ktestbedpe_LOWFILT.rds"))
#test <- readRDS(paste0(cur_dir, "/20241118_phgg_1ktestbedpe_MIDFILT.rds"))
test <- readRDS(paste0(cur_dir, "/20241118_phgg_1ktestbedpe_HIFILT.rds"))

test <- cbind(test,as.character(probabilities_radial$pred_class))

test_scaled$sv_class <- test$sv_class
colnames(test)[ncol(test)] <- 'predicted_class'

## ANALYZE PERFORMANCE #########################################################
################################################################################

# Give the option of running all of the code themselves
plot_roc_curve <- function(pr_values, title, output_path) {
  #ROC Curve
  perf <- performance(pr_values, "tpr", "fpr")
  auc_pred <- performance(pr_values, measure = "auc")
  auc_values <- auc_pred@y.values[[1]]

  # Plotting
  pdf(output_path)
  par(mgp=c(2.5,1,0))
  gg <- plot(perf, main = paste0(title), colorize = F,cex.lab=1.8, cex.main=1.8, ylim=c(0,1), xlim=c(0,1)) +
    text(0.8,0.1, labels = paste0("AUC:", substr(auc_values, 1,5)), cex=1.5) +
    abline(a = 0, b = 1)
  gg
  dev.off()
}

plot_pr_curve <- function(pr_values, title, output_path) {
  # PR Curve
  perf <- performance(pr_values, "prec", "rec")
  auc_pred <- performance(pr_values, measure = "aucpr")
  auc_values <- auc_pred@y.values[[1]]

  # Plotting
  pdf(output_path)
  par(mgp=c(2.5,1,0))
  gg <- plot(perf, main = paste0(title), colorize = F,cex.lab=1.8, cex.main=1.8, ylim=c(0,1), xlim=c(0,1)) +
    text(0.2,0.1, labels = paste0("AUC:", substr(auc_values, 1,5)), cex=1.5) +
    abline(a = 1, b = -1)
  gg
  dev.off()
}

test_cutoffs <- function(probabilities_radial, test_scaled) {
  max_sum <- 0
  max_sens <- 0
  max_ppv <- 0
  max_cutoff <- 0

  for (cutoff_prob in seq(from = 0.06, to = 0.94, by = 0.002)) {
    probabilities_radial$pred_class <- ifelse(probabilities_radial$`1` >= cutoff_prob, 1, 0)

    truth <- as.factor(test_scaled$sv_class)
    estimate <- as.factor(probabilities_radial$pred_class)

    # PRECISION RECALL
    cur_sens <- sens_vec(truth, estimate, event_level = "second")
    cur_ppv <- ppv_vec(truth, estimate, event_level = "second")

    cur_max_sum <- cur_sens + cur_ppv

    if (cur_ppv > max_ppv) {
      max_sum <- cur_max_sum
      max_sens <- cur_sens
      max_ppv <- cur_ppv
      max_cutoff <- cutoff_prob
    }
  }
  print(paste0("BEST SENSITIVITY: ", max_sens))
  print(paste0("BEST PPV: ", max_ppv))
  print(paste0("BEST CUTOFF: ", max_cutoff))

  return(max_cutoff)
}

pr_values <- prediction(as.numeric(probabilities_radial$`1`), as.numeric(test_scaled$sv_class))

probabilities_table <- data.table(cbind("prob"=probabilities_radial$`1`,"class"=test_scaled$sv_class))

plot_roc_curve(pr_values, paste0("GaTSV pHGG Manta ROC HIFILT"), paste0(cur_dir, "HIFILT_phgg_manta", "_roc_performance.pdf"))
plot_pr_curve(pr_values, paste0("GaTSV pHGG Manta PR HIFILT"), paste0(cur_dir, "HIFILT_phgg_manta", "_pr_performance.pdf"))


library(yardstick)

probabilities_radial$pred_class <- ifelse(probabilities_radial$`1` >= 0.2684, 1, 0)

truth <- as.factor(test_scaled$sv_class)
estimate <- as.factor(probabilities_radial$pred_class)

sens_vec(truth, estimate, event_level = "second")
spec_vec(truth, estimate, event_level = "second")
ppv_vec(truth, estimate, event_level = "second")


############################

test_cutoffs(probabilities_radial, test_scaled)
