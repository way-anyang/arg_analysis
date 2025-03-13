# 从命令行参数获取 ARG（替代交互输入）
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("请通过命令行参数指定目标ARG，例如：Rscript arg_analysis.R 'tetA'", call. = FALSE)
}
arg <- args[1]

# 读取数据文件（确保Linux路径正确）
rgi_data <- read.delim("/data/way/gggenes_rgi_preparation.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
mge_data <- read.delim("/data/way/gggenes_mge_preparation.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# 筛选包含目标ARG的行
rgi_matches <- rgi_data[grep(arg, do.call(paste, rgi_data), ignore.case = TRUE), ]

# 创建临时数据存储
temp_data <- data.frame(
  contig = character(),
  gene = character(),
  start = numeric(),
  end = numeric(),
  strand = character(),
  direction = numeric(),
  stringsAsFactors = FALSE
)

# 处理匹配的RGI行
if (nrow(rgi_matches) > 0) {
  for (i in 1:nrow(rgi_matches)) {
    # 获取当前RGI行
    rgi_row <- rgi_matches[i, ]
    
    # 创建基础行（来自RGI数据）
    base_row <- data.frame(
      contig = rgi_row$V3,
      gene = rgi_row$V7,
      start = rgi_row$V4,
      end = rgi_row$V5,
      strand = rgi_row$V6,
      direction = ifelse(rgi_row$V6 == "+", 1, -1),
      stringsAsFactors = FALSE
    )
    
    # 在MGE数据中查找匹配行
    mge_match <- subset(mge_data,
                        V1 == rgi_row$V1 &
                          V2 == rgi_row$V2 &
                          V3 == rgi_row$V3)
    
    # 处理MGE匹配结果
    if (nrow(mge_match) > 0) {
      for (j in 1:nrow(mge_match)) {
        # 使用MGE数据覆盖前五列
        temp_data <- rbind(temp_data, data.frame(
          contig = mge_match[j, ]$V3,
          gene = mge_match[j, ]$V4,
          start = mge_match[j, ]$V5,
          end = mge_match[j, ]$V6,
          strand = mge_match[j, ]$V7,
          direction = ifelse(mge_match[j, ]$V7 == "+", 1, -1),
          stringsAsFactors = FALSE
        ))
      }
    } else {
      # 保留RGI数据
      temp_data <- rbind(temp_data, base_row)
    }
    
    # 新增：在RGI数据中查找前三列匹配的行
    rgi_matching_rows <- subset(rgi_data,
                                V1 == rgi_row$V1 &
                                  V2 == rgi_row$V2 &
                                  V3 == rgi_row$V3)
    
    # 新增：将匹配的行添加到临时数据中
    if (nrow(rgi_matching_rows) > 0) {
      for (j in 1:nrow(rgi_matching_rows)) {
        temp_data <- rbind(temp_data, data.frame(
          contig = rgi_matching_rows[j, ]$V3,  # 第三列
          gene = rgi_matching_rows[j, ]$V7,    # 第七列
          start = rgi_matching_rows[j, ]$V4,   # 第四列
          end = rgi_matching_rows[j, ]$V5,     # 第五列
          strand = rgi_matching_rows[j, ]$V6,  # 第六列
          direction = ifelse(rgi_matching_rows[j, ]$V6 == "+", 1, -1),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
} else {
  cat("No matching ARG found in RGI data.\n")
}

write.table(temp_data, 
            "/data/way/temp.csv",
            sep = ",",
            row.names = FALSE,
            quote = TRUE,
            col.names = TRUE)
            
            
library(ggplot2)
library(gggenes)
library(GenomicRanges)
library(ape)
library(ggtree)
library(patchwork)

# 读取基因组数据
data <- read.csv("/data/way/temp.csv")

# 创建 GenomicRanges 对象
gr <- GRanges(
  seqnames = data$contig,
  ranges = IRanges(start = data$start, end = data$end),
  strand = data$strand,
  gene_id = data$gene
)

# 将 GenomicRanges 对象转换为 data.frame
gr_df <- as.data.frame(gr)

# 标记冷色基因和暖色基因
cold_genes <- c("In_pC52_003-CP042548", "ISFinder_ISKox2", "Tn6369-CP017672", 
                "Tn125-JN872328", "Tn6924-CP082952", "ISCR21-MF072961", "In37-AY259086")
gr_df$gene_type <- ifelse(gr_df$gene_id %in% cold_genes, "cold", "warm")

# 处理 arg 对齐位置
ndm1_positions <- gr_df[gr_df$gene_id == arg, c("seqnames", "start", "strand")]
colnames(ndm1_positions) <- c("contig", "ndm1_start", "ndm1_strand")
gr_df <- merge(gr_df, ndm1_positions, by.x = "seqnames", by.y = "contig", all.x = TRUE)

# 调整基因相对位置
gr_df$aligned_start <- ifelse(
  gr_df$ndm1_strand == "-",
  gr_df$ndm1_start - gr_df$end,
  gr_df$start - gr_df$ndm1_start
)
gr_df$aligned_end <- ifelse(
  gr_df$ndm1_strand == "-",
  gr_df$ndm1_start - gr_df$start,
  gr_df$end - gr_df$ndm1_start
)

# 生成颜色映射
cold_colors <- colorRampPalette(c("#2c7bb6", "#00a6ca", "#00ccbc", "#90eb9d"))(7)
names(cold_colors) <- cold_genes
other_genes <- setdiff(unique(gr_df$gene_id), cold_genes)
warm_colors <- colorRampPalette(c("#d73027", "#f46d43", "#fdae61", "#fee090"))(length(other_genes))
names(warm_colors) <- other_genes
color_mapping <- c(cold_colors, warm_colors)

# 生成基因组图
p_genes <- ggplot() +
  geom_gene_arrow(
    data = subset(gr_df, gene_type == "cold"),
    aes(xmin = aligned_start, xmax = aligned_end, y = seqnames, fill = gene_id,
        forward = strand == "+"),
    arrowhead_height = unit(3, "mm"),
    arrowhead_width = unit(1.5, "mm")
  ) +
  geom_gene_arrow(
    data = subset(gr_df, gene_type == "warm"),
    aes(xmin = aligned_start, xmax = aligned_end, y = seqnames, fill = gene_id,
        forward = strand == "+"),
    arrowhead_height = unit(3, "mm"),
    arrowhead_width = unit(1.5, "mm")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = color_mapping, breaks = c(cold_genes, other_genes)) +
  theme_genes() +
  labs(x = "Relative Position to NDM-1", y = "Contig") +  # 修改 y 轴标签
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),  # 显示 y 轴文本
    axis.ticks.y = element_line(),  # 显示 y 轴刻度
    legend.key.size = unit(0.4, "cm")
  )

# 调整 y 轴间距
p_genes <- p_genes + scale_y_discrete(expand = expansion(add = 0.9))

# 合并图形
combined_plot <- p_genes + 
  plot_layout(widths = c(3, 0.1), heights = c(13, 2)) +
  plot_annotation(theme = theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))

# 保存结果（放宽尺寸限制）
ggsave(
  "/data/way/combined_plot.pdf",
  combined_plot,
  width = 14, 
  height = 6,
  limitsize = FALSE  # 添加此行以绕过尺寸检查
)