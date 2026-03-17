#!/usr/bin/env Rscript
# 诊断脚本：查找 Phenethylamine 在数据处理过程中哪一步丢失

library(readxl)
library(dplyr)

cat("=== Phenethylamine 诊断报告 ===\n\n")

# 1. 检查原始 Excel 文件
cat("Step 1: 检查原始 Excel 文件...\n")

metab_neg <- read_excel("C:/Users/ASUS/Desktop/Gut/Metabolites quantification profiles in negative ion modes.xlsx")
metab_pos <- read_excel("C:/Users/ASUS/Desktop/Gut/Metabolites quantification profiles in positive ion modes.xlsx")

# 搜索 Phenethylamine（不区分大小写，模糊匹配）
search_phenethyl <- function(df, mode_name) {
  cat(sprintf("\n在 %s 中搜索...\n", mode_name))

  # 在 Name 列中搜索
  if ("Name" %in% colnames(df)) {
    matches_name <- grep("phenethyl", df$Name, ignore.case = TRUE, value = TRUE)
    if (length(matches_name) > 0) {
      cat(sprintf("  在 Name 列找到 %d 个匹配:\n", length(matches_name)))
      print(unique(matches_name))
    } else {
      cat("  Name 列中未找到\n")
    }
  }

  # 在 Compound_ID 列中搜索
  if ("Compound_ID" %in% colnames(df)) {
    matches_id <- grep("phenethyl", df$Compound_ID, ignore.case = TRUE, value = TRUE)
    if (length(matches_id) > 0) {
      cat(sprintf("  在 Compound_ID 列找到 %d 个匹配:\n", length(matches_id)))
      print(unique(matches_id))
    } else {
      cat("  Compound_ID 列中未找到\n")
    }
  }

  # 精确匹配 "Phenethylamine"
  exact_match <- FALSE
  if ("Name" %in% colnames(df)) {
    exact_match <- any(df$Name == "Phenethylamine", na.rm = TRUE)
  }

  if (exact_match) {
    cat("  ✓ 找到精确匹配 'Phenethylamine'\n")
    idx <- which(df$Name == "Phenethylamine")
    pea_row <- df[idx[1], ]

    # 检查数据质量
    numeric_cols <- sapply(pea_row, is.numeric)
    pea_values <- as.numeric(pea_row[numeric_cols])

    cat(sprintf("    样本数: %d\n", sum(numeric_cols)))
    cat(sprintf("    非零值数量: %d (%.1f%%)\n",
                sum(pea_values > 0, na.rm = TRUE),
                100 * sum(pea_values > 0, na.rm = TRUE) / sum(numeric_cols)))
    cat(sprintf("    缺失值数量: %d (%.1f%%)\n",
                sum(pea_values == 0, na.rm = TRUE),
                100 * sum(pea_values == 0, na.rm = TRUE) / sum(numeric_cols)))
    cat(sprintf("    均值: %.2f\n", mean(pea_values, na.rm = TRUE)))
    cat(sprintf("    中位数: %.2f\n", median(pea_values, na.rm = TRUE)))

    return(TRUE)
  }

  return(FALSE)
}

found_neg <- search_phenethyl(metab_neg, "Negative mode")
found_pos <- search_phenethyl(metab_pos, "Positive mode")

if (!found_neg && !found_pos) {
  cat("\n⚠️  警告: 在两个原始文件中都未找到 Phenethylamine!\n")
  cat("    可能的原因:\n")
  cat("    1. 该代谢物未被检测到\n")
  cat("    2. 使用了不同的名称（同义词）\n")
  cat("    3. 在数据预处理阶段已被移除\n\n")

  # 列出所有包含 "amine" 的代谢物
  cat("列出所有包含 'amine' 的代谢物作为参考:\n")

  all_amines <- c()
  if ("Name" %in% colnames(metab_neg)) {
    amines_neg <- grep("amine", metab_neg$Name, ignore.case = TRUE, value = TRUE)
    all_amines <- c(all_amines, paste0(amines_neg, " (Neg)"))
  }
  if ("Name" %in% colnames(metab_pos)) {
    amines_pos <- grep("amine", metab_pos$Name, ignore.case = TRUE, value = TRUE)
    all_amines <- c(all_amines, paste0(amines_pos, " (Pos)"))
  }

  all_amines <- unique(all_amines)
  if (length(all_amines) > 0) {
    cat(sprintf("  找到 %d 个包含 'amine' 的代谢物:\n", length(all_amines)))
    for (i in 1:min(20, length(all_amines))) {
      cat(sprintf("    %d. %s\n", i, all_amines[i]))
    }
    if (length(all_amines) > 20) {
      cat(sprintf("    ... 还有 %d 个\n", length(all_amines) - 20))
    }
  }
}

# 2. 模拟数据处理流程，追踪 Phenethylamine
cat("\n\nStep 2: 模拟数据处理流程...\n")

process_metab <- function(df) {
  df$CompoundName <- ifelse(is.na(df$Name), df$Compound_ID, df$Name)
  rownames(df) <- df$CompoundName
  numeric_cols <- sapply(df, is.numeric)
  df_num <- df[, numeric_cols]
  return(as.data.frame(t(df_num)))
}

metab_neg_t <- process_metab(metab_neg)
metab_pos_t <- process_metab(metab_pos)

cat(sprintf("处理后: Negative mode %d 特征, Positive mode %d 特征\n",
            ncol(metab_neg_t), ncol(metab_pos_t)))

# 检查是否存在
pea_in_neg <- "Phenethylamine" %in% colnames(metab_neg_t)
pea_in_pos <- "Phenethylamine" %in% colnames(metab_pos_t)

cat(sprintf("  Phenethylamine 在 Negative mode: %s\n", pea_in_neg))
cat(sprintf("  Phenethylamine 在 Positive mode: %s\n", pea_in_pos))

# 合并
metab_combined <- cbind(metab_neg_t, metab_pos_t)
cat(sprintf("\n合并后: %d 特征\n", ncol(metab_combined)))

pea_in_combined <- "Phenethylamine" %in% colnames(metab_combined)
cat(sprintf("  Phenethylamine 在合并数据: %s\n", pea_in_combined))

if (pea_in_combined) {
  # 检查是否会被去重删除
  is_dup <- duplicated(colnames(metab_combined))
  pea_indices <- which(colnames(metab_combined) == "Phenethylamine")

  if (length(pea_indices) > 1) {
    cat(sprintf("    ⚠️  Phenethylamine 出现 %d 次（重复）\n", length(pea_indices)))
    cat(sprintf("    去重后会保留第 %d 个出现的\n", pea_indices[1]))
  }
}

# 去重
metab_dedup <- metab_combined[, !duplicated(colnames(metab_combined))]
cat(sprintf("\n去重后: %d 特征\n", ncol(metab_dedup)))

pea_in_dedup <- "Phenethylamine" %in% colnames(metab_dedup)
cat(sprintf("  Phenethylamine 在去重数据: %s\n", pea_in_dedup))

if (pea_in_dedup) {
  # 检查缺失值比例
  pea_values <- metab_dedup[, "Phenethylamine"]
  missing_pct <- sum(pea_values == 0) / length(pea_values)

  cat(sprintf("    缺失值比例: %.1f%%\n", missing_pct * 100))

  if (missing_pct > 0.5) {
    cat("    ✗ 会被 >50% 缺失值过滤器移除!\n")
  } else {
    cat("    ✓ 通过 >50% 缺失值过滤\n")

    # 检查检出率
    detection_rate <- sum(pea_values > 0) / length(pea_values)
    cat(sprintf("    检出率: %.1f%%\n", detection_rate * 100))

    if (detection_rate < 0.1) {
      cat("    ✗ 会被 <10% 检出率过滤器移除!\n")
    } else {
      cat("    ✓ 通过 <10% 检出率过滤\n")
    }
  }
}

# 应用过滤器
missing_pct <- colSums(metab_dedup == 0) / nrow(metab_dedup)
metab_filtered <- metab_dedup[, missing_pct <= 0.5]

cat(sprintf("\n缺失值过滤后 (≤50%%): %d 特征\n", ncol(metab_filtered)))

pea_after_missing <- "Phenethylamine" %in% colnames(metab_filtered)
cat(sprintf("  Phenethylamine 存在: %s\n", pea_after_missing))

detection_rate <- colSums(metab_filtered > 0) / nrow(metab_filtered)
metab_filtered <- metab_filtered[, detection_rate >= 0.1]

cat(sprintf("\n检出率过滤后 (≥10%%): %d 特征\n", ncol(metab_filtered)))

pea_after_detection <- "Phenethylamine" %in% colnames(metab_filtered)
cat(sprintf("  Phenethylamine 存在: %s\n", pea_after_detection))

# 最终诊断
cat("\n\n=== 诊断结论 ===\n")

if (!found_neg && !found_pos) {
  cat("❌ Phenethylamine 在原始 Excel 文件中不存在\n")
  cat("   建议: 检查是否使用了不同的化合物名称或 ID\n")
} else if (!pea_in_combined) {
  cat("❌ Phenethylamine 在数据转置/合并过程中丢失\n")
  cat("   建议: 检查 Name 和 Compound_ID 列是否有空值\n")
} else if (!pea_in_dedup) {
  cat("❌ Phenethylamine 在去重过程中被删除（不太可能）\n")
} else if (!pea_after_missing) {
  cat("❌ Phenethylamine 因缺失值 >50% 被过滤\n")
  cat("   建议: 放宽缺失值阈值或手动保留该特征\n")
} else if (!pea_after_detection) {
  cat("❌ Phenethylamine 因检出率 <10% 被过滤\n")
  cat("   建议: 放宽检出率阈值或手动保留该特征\n")
} else {
  cat("✓ Phenethylamine 应该在最终数据中\n")
  cat("  可能是名称不完全匹配（空格、大小写等）\n")
}

# 列出最终数据中所有可能相关的代谢物
cat("\n最终数据中包含 'amine' 或 'phenyl' 的代谢物:\n")
related <- grep("amine|phenyl", colnames(metab_filtered), ignore.case = TRUE, value = TRUE)
if (length(related) > 0) {
  for (i in 1:length(related)) {
    cat(sprintf("  %d. %s\n", i, related[i]))
  }
} else {
  cat("  （无）\n")
}

cat("\n诊断完成!\n")
