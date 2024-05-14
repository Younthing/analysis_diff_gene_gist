fs::dir_create("Diff_analysis")
rm(list = ls())
# 数据预处理------------------------------------------------------------

## 处理矩阵
expr_file <- "./2-GEO校正/batch_normal_matrix.csv"

expr_data <- read.csv(expr_file, row.names = 1)

## 前面已经做了log2处理,系统效应和方差分量通常被假定为在对数尺度上相加
## （Oberg AL. et al JPR 2008；Hill EG. et al JPR 2008）
boxplot(expr_data, las = 2, ylab = "Intensity")

## 处理分组
group_file <- "./2-GEO校正/combined_group.csv"
group_data <- read.csv(group_file, row.names = NULL)

group_data$sample <- make.names(group_data$sample) # 从字符向量中生成语法上有效的名称。

## 确保矩阵样本和分组顺序对应
expr_data <- expr_data[, group_data$sample]

# 差异表达分析------------------------------------------------------------

## 实验设计矩阵

cond <- group_data$group

design <- model.matrix(~ 0 + cond) # 0 不考虑截距
colnames(design) <- gsub("cond", "", colnames(design))

## 设置对比
contrasts_vector <- c("Disorder-Normal")
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
  write.csv(deg_results, paste0("Diff_analysis/Diff_result_", choose_trt, ".csv")) # nolint

  ## 提醒
  cat(paste0("Diff_result_", choose_trt, ".csv"), "\n")
}

saveRDS(diff_res_list, "Diff_analysis/Diff_result_list.rds")

# 火山图------------------------------------------------------------

## limma默认用FDR
## FDR 校正的 p 值与原始值相比丢失了一些信息
## Benjamini-Hochberg和p是1:1的关系，没有权重等其他的话
## 所以火山图推荐用原始的p值
library(ggplot2)
library(ggrepel)
library(dplyr)
logfc_cut <- 0
argparser::include("./0-my_theme.R")

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

  ggsave(paste0("Diff_analysis/Volcano_plot_", choose_trt, ".pdf"), width = 6, height = 6, units = "cm")
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
  write.csv(tmp$gene, paste0("Diff_analysis/top_sig_gene_list_", choose_trt, ".csv"), row.names = FALSE)
}

# 保存显著差异基因

sig_gene_list <- lapply(diff_res_list, function(x) {
  x <- x[x$P.Val < 0.05, ]
  x <- x[x$logFC > abs(logfc_cut) | x$logFC < -abs(logfc_cut), ]
  return(x$gene)
})

saveRDS(sig_gene_list, "Diff_analysis/sig_gene_list.rds")


## 循环保存基因
lapply(seq_along(sig_gene_list), function(i) {
  choose_trt <- names(sig_gene_list)[i]
  write.csv(
    sig_gene_list[[i]],
    paste0("Diff_analysis/sig_gene_list_", choose_trt, ".csv"),
    row.names = FALSE
  )

  cat(paste0("sig_gene_list_", choose_trt, ".csv"), "\n")
})
