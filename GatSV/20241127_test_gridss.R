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

# get all processed bed files and put into one DF
# Set the working directory
setwd("/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/GetScaled")
# Get a list of all subdirectories
subdirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)

# Define a function to process a single subdirectory
process_subdir <- function(subdir) {
  # Get the list of files in the subdirectory
  files <- list.files(path = subdir, full.names = TRUE)

  # Filter for the file containing "processed" in its name
  processed_file <- files[grep("processed", files)]

  # If a processed file exists, read it and return the data
  if (length(processed_file) > 0) {
    return(fread(processed_file[1]))  # Assuming only one "processed" file per subdirectory
  } else {
    return(NULL)  # Return NULL if no file found
  }
}

# Use mclapply to process subdirectories in parallel
results <- mclapply(subdirs, process_subdir, mc.cores = detectCores())

# Combine all non-NULL results into a single data.table
all_proc_ <- rbindlist(results, fill = TRUE)

# Result: all_proc contains concatenated data from all "processed" files

# Assuming 'GaTSV' is your trained SVM model, 'test_scaled' is your test data
# Number of cores to use for parallel processing
num_cores <- detectCores() - 1  # Use one less than total cores to avoid overloading

# Step 1: Split the test data into chunks
chunk_size <- ceiling(nrow(test_scaled) / num_cores)
chunks <- split(test_scaled, rep(1:num_cores, each = chunk_size, length.out = nrow(test_scaled)))

# Step 2: Run predictions in parallel using mclapply
y_pred_radial_list <- pbmclapply(chunks, function(chunk) {
  predict(GaTSV, newdata = chunk, decision.values = TRUE, probability = TRUE)
}, mc.cores = num_cores)

# Step 3: Extract and process probabilities for each prediction
probabilities_list <- lapply(y_pred_radial_list, function(y_pred) {
  # Extract the probabilities and convert to data.table
  probabilities_radial <- data.table(attr(y_pred, 'probabilities'))
  
  # Ensure the columns are in the correct order
  setcolorder(probabilities_radial, c('0', '1'))
  
  return(probabilities_radial)
})

# Step 4: Concatenate all the data.tables into one
final_probabilities <- rbindlist(probabilities_list)




## TEST SVM ####################################################################
################################################################################
GaTSV <- readRDS("/data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/GatSV/svm/GaTSV.rda") #svmobject
common <- colnames(GaTSV$SV)

test_scaled <- all_proc[, ..common]

y_pred_radial <- predict(GaTSV, newdata = test_scaled, decision.values = T, probability = T)

probabilities_radial <-data.table(attr(y_pred_radial, 'probabilities'))
setcolorder(probabilities_radial, c('0', '1'))

probabilities_radial<-final_probabilities

cutoff_prob <-  0.2684
probabilities_radial[,"pred_class"] <- lapply(1:length(probabilities_radial$`1`),function(i){
  return (ifelse(probabilities_radial$`1`[i]>= cutoff_prob,'GERMLINE','SOMATIC'))})

#test <- readRDS(paste0(cur_dir, "/20241118_phgg_1ktestbedpe_LOWFILT.rds"))
#test <- readRDS(paste0(cur_dir, "/20241118_phgg_1ktestbedpe_MIDFILT.rds"))
test <- readRDS(paste0(cur_dir, "/20241118_phgg_1ktestbedpe_HIFILT.rds"))

test <- all_proc

test <- cbind(test,as.character(probabilities_radial$pred_class))
test[, sv_class := ifelse(FILTER == "PON", 1, 0)]
test[, sv_class_pred := ifelse(predicted_class == "GERMLINE", 1, 0)]

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
    print(cutoff_prob)
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

plot_roc_curve(pr_values, paste0("GaTSV HMF GRIDSS ROC"), paste0("/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/Plots/HMF_roc_performance_inv.pdf"))
plot_pr_curve(pr_values, paste0("GaTSV HMF GRIDSS PR"), paste0("/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/Plots/HMF_pr_performance_inv.pdf"))


library(yardstick)

probabilities_radial$pred_class <- ifelse(probabilities_radial$`1` >= 0.2684, 1, 0)

truth <- as.factor(test$sv_class)
estimate <- as.factor(test$sv_class_pred)

sens_vec(truth, estimate, event_level = "second")
spec_vec(truth, estimate, event_level = "second")
ppv_vec(truth, estimate, event_level = "second")


############################

test_cutoffs(probabilities_radial, test_scaled)
ewfef=test_cutoffs_(probabilities_radial, test_scaled, "/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/Plots/")

test_cutoffs_ <- function(probabilities_radial, test_scaled, save_dir = ".") {
  # Initialize an empty data frame to store results
  results <- data.frame(cutoff = numeric(0), sensitivity = numeric(0), ppv = numeric(0))
  
  max_sum <- 0
  max_sens <- 0
  max_ppv <- 0
  max_cutoff <- 0

  for (cutoff_prob in seq(from = 0.06, to = 0.99, by = 0.001)) {
    print(cutoff_prob)
    probabilities_radial$pred_class <- ifelse(probabilities_radial$`1` >= cutoff_prob, 1, 0)

    truth <- as.factor(test_scaled$sv_class)
    estimate <- as.factor(probabilities_radial$pred_class)

    # Compute Precision and Recall
    cur_sens <- sens_vec(truth, estimate, event_level = "second")
    cur_ppv <- ppv_vec(truth, estimate, event_level = "second")

    # Save the results for each cutoff value
    results <- rbind(results, data.frame(cutoff = cutoff_prob, sensitivity = cur_sens, ppv = cur_ppv))

    cur_max_sum <- cur_sens + cur_ppv

    if (cur_ppv > max_ppv) {
      max_sum <- cur_max_sum
      max_sens <- cur_sens
      max_ppv <- cur_ppv
      max_cutoff <- cutoff_prob
    }
  }

  # Print the best cutoff values
  print(paste0("BEST SENSITIVITY: ", max_sens))
  print(paste0("BEST PPV: ", max_ppv))
  print(paste0("BEST CUTOFF: ", max_cutoff))

  # Save plots to specified directory
  # Ensure directory exists
  if (!dir.exists(save_dir)) {
    dir.create(save_dir)
  }

  # Plot Sensitivity vs Cutoff and save to PDF
  pdf(file.path(save_dir, "sensitivity_vs_cutoff.pdf"))
  ggplot(results, aes(x = cutoff, y = sensitivity)) +
    geom_line() +
    labs(title = "Sensitivity vs Cutoff", x = "Cutoff", y = "Sensitivity") +
    theme_minimal()
  dev.off()

  # Plot PPV vs Cutoff and save to PDF
  pdf(file.path(save_dir, "ppv_vs_cutoff.pdf"))
  ggplot(results, aes(x = cutoff, y = ppv)) +
    geom_line() +
    labs(title = "PPV vs Cutoff", x = "Cutoff", y = "PPV") +
    theme_minimal()
  dev.off()

  # Return the best cutoff value
  return(max_cutoff)
}


# Step 1: Count occurrences of "PON" and all other values
value_counts <- all_proc_[, .(count = .N), by = FILTER]
# Create a new column for categorization: "PON" or "Others"
value_counts[, category := ifelse(FILTER == "PON", "PON", "PASS")]

# Step 2: Summarize the counts for "PON" and "Others"
summary_counts <- value_counts[, .(count = sum(count)), by = category]

# Step 3: Calculate proportions
summary_counts[, proportion := count / sum(count) * 100]

# Step 4: Create the pie chart with labels for the proportions
pie_plot <- ggplot(summary_counts, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +  # Convert bar chart to pie chart
  labs(title = "PON vs PASS after preprocessing") +
  theme_void() +  # Removes background grid and axes
  scale_fill_manual(values = c("PON" = "blue", "PASS" = "gray")) +
  geom_text(aes(label = paste0(round(proportion, 1), "%")), position = position_stack(vjust = 0.5))  # Add proportions as text labels

# Step 5: Save the plot as a PDF to the specified directory
pdf_path <- "/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/Plots/pie_chart_proportion_.pdf"
ggsave(pdf_path, plot = pie_plot, width = 8, height = 6)

# Confirmation message
cat("Pie chart saved to:", pdf_path, "\n")

test_cutoffs_ <- function(probabilities_radial, test_scaled, save_dir = ".", num_cores = 32) {
  # Initialize empty variables to store the best results
  max_sum <- 0
  max_sens <- 0
  max_ppv <- 0
  max_cutoff <- 0

  # Define the function for each cutoff computation
  cutoff_function <- function(cutoff_prob) {
    # Create prediction class based on cutoff
    probabilities_radial$pred_class <- ifelse(probabilities_radial$`1` >= cutoff_prob, 1, 0)
    
    # Get true and predicted values
    truth <- as.factor(test_scaled$sv_class)
    estimate <- as.factor(probabilities_radial$pred_class)

    # Compute Precision and Recall & specificity
    cur_spec <- spec_vec(truth, estimate, event_level="second")
    cur_sens <- sens_vec(truth, estimate, event_level = "second")
    cur_ppv <- ppv_vec(truth, estimate, event_level = "second")

    # Return results as a list
    return(data.frame(cutoff = cutoff_prob, sensitivity = cur_sens, ppv = cur_ppv, specificity=cur_spec))
  }

  # Parallelize the loop using mclapply
  cutoffs <- seq(from = 0.312, to = 0.314, by = 0.00001)
  results_list <- mclapply(cutoffs, cutoff_function, mc.cores = num_cores)

  # Combine the results into a single data frame
  results <- do.call(rbind, results_list)

  # Find the best cutoff based on the highest PPV
  best_row <- results[which.max(results$ppv),]
  max_sens <- best_row$sensitivity
  max_ppv <- best_row$ppv
  max_cutoff <- best_row$cutoff

  # Print the best cutoff values
  print(paste0("BEST SENSITIVITY: ", max_sens))
  print(paste0("BEST PPV: ", max_ppv))
  print(paste0("BEST CUTOFF: ", max_cutoff))

  # Save plots to the specified directory
  if (!dir.exists(save_dir)) {
    dir.create(save_dir)
  }

  # Plot Sensitivity vs Cutoff and save to PDF
  pdf(file.path(save_dir, "sensitivity_vs_cutoff.pdf"))
  print(ggplot(results, aes(x = cutoff, y = sensitivity)) +
    geom_line() +
    geom_vline(xintercept = 0.2864, color = "red", linetype = "dashed", size = 1) +  # Add red vertical line
    labs(title = "Sensitivity vs Cutoff", x = "Cutoff", y = "Sensitivity") +
    theme_minimal())
  dev.off()

  # Plot PPV vs Cutoff and save to PDF
  pdf(file.path(save_dir, "ppv_vs_cutoff.pdf"))
  print(ggplot(results, aes(x = cutoff, y = ppv)) +
    geom_line() +
    geom_vline(xintercept = 0.2864, color = "red", linetype = "dashed", size = 1) +  # Add red vertical line
    labs(title = "PPV vs Cutoff", x = "Cutoff", y = "PPV") +
    theme_minimal())
  dev.off()

  # Plot Spec vs Cutoff and save to PDF
  pdf(file.path(save_dir, "specificity_vs_cutoff.pdf"))
  print(ggplot(results, aes(x = cutoff, y = specificity)) +
    geom_line() +
    geom_vline(xintercept = 0.2864, color = "red", linetype = "dashed", size = 1) +  # Add red vertical line
    labs(title = "Specificity vs Cutoff", x = "Cutoff", y = "Specificity") +
    theme_minimal())
  dev.off()

  # Return the best cutoff value
  return(results)
}

best_cutoff <- test_cutoffs_(probabilities_radial, test_scaled, save_dir = "/data/cephfs-1/home/users/rauertc_c/work/GatSV/Out/Plots/", num_cores = 32)