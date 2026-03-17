#!/usr/bin/env Rscript
# MintTea Multi-omics Integration Analysis (Fixed Version)
# 修复了模块提取问题，添加了更健壮的错误处理

library(MintTea)
library(dplyr)
library(readr)
library(readxl)
library(ggplot2)
library(pheatmap)

cat("Step 1: Loading data...\n")

# Load microbiome data
mic_raw <- read.csv("C:/Users/ASUS/Desktop/Gut/Gut microbiota.csv",
                    row.names = 1,
                    check.names = FALSE,
                    encoding = "UTF-8")

# Extract species-level features
extract_species <- function(tax_string) {
  if (is.na(tax_string)) return(NA)
  if (grepl("\\|s__", tax_string)) {
    parts <- strsplit(tax_string, "\\|")[[1]]
    species_part <- parts[grepl("^s__", parts)]
    if (length(species_part) > 0) {
      return(gsub("^s__", "", species_part[1]))
    }
  }
  return(NA)
}

species_names <- sapply(rownames(mic_raw), extract_species)
mic_species <- mic_raw[!is.na(species_names), ]
rownames(mic_species) <- species_names[!is.na(species_names)]

# Aggregate duplicates
mic_species <- as.data.frame(t(mic_species))
mic_species <- mic_species %>%
  group_by(row.names(.)) %>%
  summarise(across(everything(), sum)) %>%
  as.data.frame()
rownames(mic_species) <- mic_species[,1]
mic_species <- mic_species[,-1]
mic_species <- t(mic_species)

cat(sprintf("Microbiome: %d species\n", nrow(mic_species)))

# Load metabolomics
metab_neg <- read_excel("C:/Users/ASUS/Desktop/Gut/Metabolites quantification profiles in negative ion modes.xlsx")
metab_pos <- read_excel("C:/Users/ASUS/Desktop/Gut/Metabolites quantification profiles in positive ion modes.xlsx")

# Process metabolomics
process_metab <- function(df) {
  df$CompoundName <- ifelse(is.na(df$Name), df$Compound_ID, df$Name)
  rownames(df) <- df$CompoundName
  numeric_cols <- sapply(df, is.numeric)
  df_num <- df[, numeric_cols]
  return(as.data.frame(t(df_num)))
}

metab_neg_t <- process_metab(metab_neg)
metab_pos_t <- process_metab(metab_pos)

metab_combined <- cbind(metab_neg_t, metab_pos_t)
metab_dedup <- metab_combined[, !duplicated(colnames(metab_combined))]

# Filter metabolites
missing_pct <- colSums(metab_dedup == 0) / nrow(metab_dedup)
metab_filtered <- metab_dedup[, missing_pct <= 0.5]

detection_rate <- colSums(metab_filtered > 0) / nrow(metab_filtered)
metab_filtered <- metab_filtered[, detection_rate >= 0.1]

cat(sprintf("Metabolome: %d compounds\n", ncol(metab_filtered)))

# Impute missing values
for (col in colnames(metab_filtered)) {
  min_val <- min(metab_filtered[metab_filtered[,col] > 0, col], na.rm = TRUE)
  if (is.finite(min_val) && min_val > 0) {
    metab_filtered[metab_filtered[,col] == 0, col] <- min_val / 2
  } else {
    metab_filtered[metab_filtered[,col] == 0, col] <- 1e-6
  }
}

# Match samples
mic_samples <- gsub(" ", "_", colnames(mic_species))
metab_samples <- gsub(" ", "_", rownames(metab_filtered))

common_samples <- intersect(mic_samples, metab_samples)
cat(sprintf("Common samples: %d\n", length(common_samples)))

# Align data
mic_aligned <- mic_species[, match(common_samples, mic_samples)]
colnames(mic_aligned) <- common_samples
mic_aligned <- as.data.frame(t(mic_aligned))

metab_aligned <- metab_filtered[match(common_samples, metab_samples), ]
rownames(metab_aligned) <- common_samples

cat("\nStep 2: Data transformation...\n")

# CLR transformation for microbiome
clr_transform <- function(x) {
  x_pos <- x + 1e-6  # Add pseudocount
  gm <- exp(mean(log(x_pos)))
  return(log(x_pos / gm))
}

mic_clr <- as.data.frame(t(apply(mic_aligned, 1, clr_transform)))

# Log2 transformation for metabolomics
metab_log <- log2(metab_aligned + 1)

# Remove zero variance features
mic_var <- apply(mic_clr, 2, var)
mic_clr <- mic_clr[, mic_var > 0]

metab_var <- apply(metab_log, 2, var)
metab_log <- metab_log[, metab_var > 0]

cat(sprintf("After variance filtering: %d microbiome, %d metabolome features\n",
            ncol(mic_clr), ncol(metab_log)))

# Create study groups
study_group <- ifelse(grepl("NCD", common_samples), "disease", "healthy")

cat("\nStep 3: Preparing MintTea input with relaxed parameters...\n")

# Add view prefixes
colnames(mic_clr) <- paste0("T__", colnames(mic_clr))
colnames(metab_log) <- paste0("M__", colnames(metab_log))

# Create input dataframe
minttea_input <- data.frame(
  SampleID = common_samples,
  StudyGroup = study_group,
  mic_clr,
  metab_log,
  check.names = FALSE
)

cat(sprintf("Input dimensions: %d samples x %d features\n",
            nrow(minttea_input), ncol(minttea_input) - 2))

cat("\nStep 4: Running MintTea with relaxed parameters...\n")

# 使用更宽松的参数
minttea_result <- tryCatch({
  MintTea(
    proc_data = minttea_input,
    study_group_column = "StudyGroup",
    control_group_name = "healthy",
    case_group_name = "disease",
    sample_id_column = "SampleID",
    view_prefixes = c("T", "M"),
    param_edge_thresholds = 0.5,  # 降低阈值从 0.7 到 0.5
    param_repeats = 10,            # 减少重复次数
    param_num_cores = 1
  )
}, error = function(e) {
  cat("MintTea error:", e$message, "\n")
  return(NULL)
})

cat("\n=== MintTea Result Structure ===\n")
if (is.null(minttea_result)) {
  cat("ERROR: MintTea returned NULL\n")
} else {
  cat("Result class:", class(minttea_result), "\n")
  cat("Result is a list:", is.list(minttea_result), "\n")
  cat("Result names:", paste(names(minttea_result), collapse = ", "), "\n")
  cat("Result length:", length(minttea_result), "\n")
}

# 尝试多种方式提取模块
modules_df <- NULL

if (!is.null(minttea_result)) {
  # 方式1: 直接提取 modules
  if ("modules" %in% names(minttea_result)) {
    modules_df <- minttea_result$modules
    cat("Found modules in minttea_result$modules\n")
  }

  # 方式2: 检查 results 中的 modules
  if (is.null(modules_df) && "results" %in% names(minttea_result)) {
    if ("modules" %in% names(minttea_result$results)) {
      modules_df <- minttea_result$results$modules
      cat("Found modules in minttea_result$results$modules\n")
    }
  }

  # 方式3: 检查第一个元素
  if (is.null(modules_df) && length(minttea_result) > 0) {
    first_elem <- minttea_result[[1]]
    if (is.data.frame(first_elem) || is.matrix(first_elem)) {
      modules_df <- as.data.frame(first_elem)
      cat("Found data frame/matrix in first element\n")
    }
  }
}

cat("\n=== Module Extraction Result ===\n")
if (is.null(modules_df)) {
  cat("WARNING: No modules found in MintTea result!\n")
  cat("This may be due to:\n")
  cat("  1. Small sample size (n=", length(common_samples), ")\n")
  cat("  2. Weak multi-omics associations\n")
  cat("  3. Strict parameter thresholds\n")
  cat("\nProceeding with correlation-based analysis instead...\n")
} else if (nrow(modules_df) == 0) {
  cat("WARNING: modules_df has 0 rows!\n")
  cat("MintTea found no significant modules.\n")
  modules_df <- NULL
} else {
  cat("SUCCESS: Found", nrow(modules_df), "modules\n")
  cat("Module columns:", paste(colnames(modules_df), collapse = ", "), "\n")
  print(head(modules_df))
}

cat("\nStep 5: Generating visualizations...\n")

# 生成可视化（即使没有模块也能运行）

# 1. 相关性热图
cat("Generating correlation heatmap...\n")

# 选择关键特征（方差最大的前50个）
mic_top <- mic_clr[, order(apply(mic_clr, 2, var), decreasing = TRUE)[1:min(50, ncol(mic_clr))]]
metab_top <- metab_log[, order(apply(metab_log, 2, var), decreasing = TRUE)[1:min(50, ncol(metab_log))]]

# 计算交叉相关性
cor_matrix <- cor(mic_top, metab_top, use = "pairwise.complete.obs")

# 清理特征名（移除前缀）
rownames(cor_matrix) <- gsub("^T__", "", rownames(cor_matrix))
colnames(cor_matrix) <- gsub("^M__", "", colnames(cor_matrix))

pdf("C:/Users/ASUS/Desktop/Gut/minttea_correlation_heatmap.pdf", width = 12, height = 10)
pheatmap(cor_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         main = "Microbiome-Metabolome Correlation (Top 50 features each)",
         fontsize = 6)
dev.off()

cat("Saved: minttea_correlation_heatmap.pdf\n")

# 2. 如果有模块，绘制模块特异性图
if (!is.null(modules_df) && nrow(modules_df) > 0) {
  cat("Generating module-specific visualizations...\n")

  pdf("C:/Users/ASUS/Desktop/Gut/minttea_modules.pdf", width = 10, height = 8)

  for (i in 1:min(5, nrow(modules_df))) {  # 最多展示5个模块
    cat(sprintf("Processing module %d/%d\n", i, nrow(modules_df)))

    module_id <- modules_df[i, "ModuleID"]
    module_features <- strsplit(as.character(modules_df[i, "Features"]), ";")[[1]]

    # 区分微生物和代谢物特征
    mic_features <- grep("^T__", module_features, value = TRUE)
    metab_features <- grep("^M__", module_features, value = TRUE)

    cat(sprintf("  Module %s: %d microbiome + %d metabolome features\n",
                module_id, length(mic_features), length(metab_features)))

    # 绘制该模块内部相关性
    if (length(mic_features) > 0 && length(metab_features) > 0) {
      module_mic <- mic_clr[, mic_features, drop = FALSE]
      module_metab <- metab_log[, metab_features, drop = FALSE]

      module_cor <- cor(module_mic, module_metab, use = "pairwise.complete.obs")

      pheatmap(module_cor,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               color = colorRampPalette(c("blue", "white", "red"))(100),
               main = paste("Module", module_id, "Correlation"),
               fontsize = 8)
    }
  }

  dev.off()
  cat("Saved: minttea_modules.pdf\n")
}

# 3. PCA 可视化
cat("Generating PCA plots...\n")

library(ggfortify)

combined_data <- cbind(mic_clr, metab_log)
pca_result <- prcomp(combined_data, scale. = TRUE)

pdf("C:/Users/ASUS/Desktop/Gut/minttea_pca.pdf", width = 10, height = 8)

pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Group = study_group,
  Sample = common_samples
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() +
  labs(title = "PCA of Combined Multi-omics Data",
       x = sprintf("PC1 (%.1f%%)", summary(pca_result)$importance[2,1] * 100),
       y = sprintf("PC2 (%.1f%%)", summary(pca_result)$importance[2,2] * 100)) +
  scale_color_manual(values = c("healthy" = "blue", "disease" = "red"))

dev.off()
cat("Saved: minttea_pca.pdf\n")

# 4. 检查目标特征
cat("\n=== Checking Target Features ===\n")

target_mic <- "T__Ruminococcus gnavus"
target_metab <- "M__Phenethylamine"

if (target_mic %in% colnames(mic_clr)) {
  cat("✓ Found Ruminococcus gnavus\n")

  # 找到与它相关性最强的代谢物
  rum_cor <- cor(mic_clr[, target_mic], metab_log, use = "pairwise.complete.obs")
  top_metab <- sort(abs(rum_cor[1,]), decreasing = TRUE)[1:10]

  cat("  Top 10 correlated metabolites:\n")
  for (j in 1:length(top_metab)) {
    metab_name <- gsub("^M__", "", names(top_metab)[j])
    cat(sprintf("    %d. %s (r = %.3f)\n", j, metab_name, rum_cor[1, names(top_metab)[j]]))
  }

  if (target_metab %in% names(rum_cor)) {
    cat(sprintf("\n  Correlation with Phenethylamine: r = %.3f\n",
                rum_cor[1, target_metab]))
  }
} else {
  cat("✗ Ruminococcus gnavus not found in microbiome data\n")

  # 列出所有 Ruminococcus 物种
  rum_species <- grep("T__Ruminococcus", colnames(mic_clr), value = TRUE)
  if (length(rum_species) > 0) {
    cat("  Available Ruminococcus species:\n")
    for (sp in rum_species) {
      cat(sprintf("    - %s\n", gsub("^T__", "", sp)))
    }
  }
}

if (target_metab %in% colnames(metab_log)) {
  cat("\n✓ Found Phenethylamine\n")

  # 找到与它相关性最强的微生物
  pea_cor <- cor(metab_log[, target_metab], mic_clr, use = "pairwise.complete.obs")
  top_mic <- sort(abs(pea_cor[1,]), decreasing = TRUE)[1:10]

  cat("  Top 10 correlated microbes:\n")
  for (j in 1:length(top_mic)) {
    mic_name <- gsub("^T__", "", names(top_mic)[j])
    cat(sprintf("    %d. %s (r = %.3f)\n", j, mic_name, pea_cor[1, names(top_mic)[j]]))
  }
} else {
  cat("\n✗ Phenethylamine not found in metabolome data\n")

  # 模糊搜索
  pea_like <- grep("Phenethyl", colnames(metab_log), value = TRUE, ignore.case = TRUE)
  if (length(pea_like) > 0) {
    cat("  Similar metabolites found:\n")
    for (metab in pea_like) {
      cat(sprintf("    - %s\n", gsub("^M__", "", metab)))
    }
  }
}

cat("\n=== Analysis Complete ===\n")
cat("Output files:\n")
cat("  - minttea_correlation_heatmap.pdf\n")
if (!is.null(modules_df) && nrow(modules_df) > 0) {
  cat("  - minttea_modules.pdf\n")
}
cat("  - minttea_pca.pdf\n")
