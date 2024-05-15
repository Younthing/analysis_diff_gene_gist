# ---------------------------------------------------------
#
# Title: Difference gene analysis
# Description: Differential gene analysis and visualization
# Author: fanxi
# Date: 2024-05-15
# Version:
# Change Log:
# - 2024-05-14: Initial version by fanxi
# - 2024-05-15: Add formatting
#
# Dependencies:
# - dplyr, magrittr, rlang, forecast
# - purrr::keep, stringr::str_detect
# Instructions:
# - need config_argparse_yaml.R, config_my_theme.R, param.yaml
# Usage:
# - Rscript analysis_diff_gene.R
#
# TODO: usegalaxy
# ---------------------------------------------------------

# Libraries & Dependencies --------------------------------
# 加载所需的 R 包和其他依赖
library(dplyr)
library(forcats)
library(rlang) # %||% : if x is NULL, return y
library(magrittr) # %T>%

library(limma)
rm(list = ls())

# ---------------------------------------------------------

# 传参------------------------------------------------------------

dir <- "Diff_analysis"
fs::dir_create(dir)

param_parse_script <- "config_argparse_yaml.R"
param <- argparser::include(param_parse_script) # 类似soures但是可以取得局部作用域最后的值
method_select <- "limma"

## 本地变量
local_params <- param$value$analysis$diff_analysis %T>% str()

# 函数------------------------------------------------------------

#' @param method_select 选择分析方法
#' @return name为limma的参数列表
diff_params <- purrr::keep(local_params, ~ stringr::str_detect(.$name, pattern = method_select)) %>%
  unlist(recursive = TRUE) %>%
  as.list() %T>%
  jsonlite::write_json("diff_params.json", row.names = FALSE) %T>%
  str()

## 处理矩阵
#' @param expr_file 输入表达矩阵文件
#' @param group_file 输入分组文件
expr_file <- diff_params$expr_file %||% "./2-GEO校正/batch_normal_matrix.csv"
expr_data <- read.csv(expr_file, row.names = 1)

## 处理分组
group_file <- diff_params$group_file %||% "./2-GEO校正/combined_group.csv"
group_data <- read.csv(group_file, row.names = NULL)

group_data$sample <- make.names(group_data$sample) # 从字符向量中生成语法上有效的名称。

## 确保矩阵样本和分组顺序对应
expr_data <- expr_data[, group_data$sample]

## 设置对比
#' @param contrast_vector 输入对比组
contrasts_vector <- eval(parse(text = diff_params$contrasts)) %||% c("Disorder-Normal") # 默认为“Disorder-Normal”


# 差异表达分析------------------------------------------------------------
## 前面已经做了log2处理,系统效应和方差分量通常被假定为在对数尺度上相加
## （Oberg AL. et al JPR 2008；Hill EG. et al JPR 2008）
boxplot(expr_data, las = 2, ylab = "Intensity")

## 实验设计矩阵

cond <- group_data$group

design <- model.matrix(~ 0 + cond) # 0 不考虑截距
colnames(design) <- gsub("cond", "", colnames(design))

## 设置对比
contrast <- makeContrasts(contrasts = contrasts_vector, levels = design)

fit1 <- lmFit(expr_data, design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit_res <- eBayes(fit2)
plotSA(fit_res)


# 提取结果------------------------------------------------------------

dif <- topTable(fit_res, coef = 1, n = Inf)
dif <- na.omit(dif)
dif$gene <- row.names(dif)

## 循环保存
diff_res_list <- list()
for (i in seq_along(colnames(fit_res$coefficients))) {
  choose_trt <- colnames(fit_res$coefficients)[i]
  deg_results <- topTable(
    fit_res,
    coef = i,
    n = Inf
  )
  deg_results$gene <- row.names(deg_results)
  diff_res_list[[i]] <- deg_results
  names(diff_res_list)[i] <- choose_trt
  write.csv(deg_results, fs::path(dir, paste0("Diff_result_", choose_trt, ".csv"))) # nolint

  ## 提醒
  cat(paste0("Diff_result_", choose_trt, ".csv"), "\n")
}
saveRDS(diff_res_list, fs::path(dir, "Diff_result_list.rds"))

# 火山图------------------------------------------------------------

## limma默认用FDR
## FDR 校正的 p 值与原始值相比丢失了一些信息
## Benjamini-Hochberg和p是1:1的关系，没有权重等其他的话
## 所以火山图推荐用原始的p值
library(ggplot2)
library(ggrepel)
library(dplyr)
logfc_cut <- diff_params$logfc_threshold |> as.numeric()
argparser::include(diff_params$theme_file)

for (i in seq_along(diff_res_list)) {
  choose_trt <- colnames(fit_res$coefficients)[i]
  deg_results <- diff_res_list[[i]]
  deg_results <- deg_results %>%
    mutate(label = case_when(
      P.Value > 0.05 ~ "no sig",
      logFC > abs(logfc_cut) ~ "up gene",
      logFC < -abs(logfc_cut) ~ "down gene",
      TRUE ~ "no sig"
    ))

  top_results <- deg_results |>
    dplyr::filter(P.Value < 0.05) |>
    dplyr::arrange(logFC) |>
    slice_head(n = 5)
  tail_results <- deg_results |>
    dplyr::filter(P.Value < 0.05) |>
    dplyr::arrange(logFC) |>
    slice_tail(n = 5)
  tmp <- rbind(top_results, tail_results)

  max_abs_x <- max(abs(min(deg_results$logFC)), abs(max(deg_results$logFC)))

  ggplot(deg_results, aes(x = logFC, y = -log10(P.Value), color = label)) +
    geom_point(size = 0.5) +
    xlab(expression("log2FC")) + # x-axis label
    ylab(expression(" -log10(P-value)")) + # y-axis label
    # geom_vline(xintercept = c(-0.5, 0.5), colour = "red") + # Add fold change cutoffs
    geom_hline(yintercept = -log10(0.05), colour = "red", linetype = "dashed", size = 0.25) +
    geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
    scale_colour_manual(values = c("down gene" = "blue", "up gene" = "red", "no sig" = "grey")) + # Manual color scale
    geom_text_repel(
      data = tmp,
      aes(logFC, -log10(P.Value), label = gene),
      size = 1,
      color = "black", # 设置字体颜色为黑色,
      max.overlaps = 10,
    ) +
    xlim(-max_abs_x, max_abs_x) +
    labs(color = "") +
    my_theme

  ggsave(fs::path(dir, paste0("Volcano_plot_", choose_trt, ".pdf")), width = 6, height = 6, units = "cm")
  cat(paste0("Volcano_plot_", choose_trt, ".pdf"), "\n")

  top_results <- deg_results |>
    dplyr::filter(P.Value < 0.05) |>
    dplyr::arrange(logFC) |>
    slice_head(n = 20)
  tail_results <- deg_results |>
    dplyr::filter(P.Value < 0.05) |>
    dplyr::arrange(logFC) |>
    slice_tail(n = 20)
  tmp <- rbind(top_results, tail_results)
  write.csv(tmp$gene, fs::path(dir, paste0("top_sig_gene_list_", choose_trt, ".csv")), row.names = FALSE)
}

# 保存显著差异基因

sig_gene_list <- lapply(diff_res_list, function(x) {
  x <- x[x$P.Val < 0.05, ]
  x <- x[x$logFC > abs(logfc_cut) | x$logFC < -abs(logfc_cut), ]
  return(x$gene)
})

library(fs)

# Save RDS file with fs::path()
saveRDS(diff_res_list, fs::path(dir, "diff_res_list.rds"))

## Loop to save gene lists
lapply(seq_along(sig_gene_list), function(i) {
  choose_trt <- names(sig_gene_list)[i]

  # Use fs::path() for constructing the file path with dir
  file_path <- fs::path(dir, paste0("sig_gene_list_", choose_trt, ".csv"))

  # Save the gene list as CSV
  write.csv(sig_gene_list[[i]], file_path, row.names = FALSE)

  # Print the file path
  cat(file_path, "\n")
})
