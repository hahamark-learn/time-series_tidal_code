############################################################
## Persistent vs transient HGT analysis at phylum level
## Input: WAAFLE HGT result table
############################################################

library(tidyverse)
library(ggalluvial)

############################################################
## 1. 读入数据
############################################################

hgt <- read.table(
  "all_hgt_time",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

############################################################
## 2. 拆分 UniRef90 注释
############################################################

hgt2 <- hgt %>%
  separate_rows(`ANNOTATIONS:UNIREF90`, sep = "\\|") %>%
  filter(
    !is.na(`ANNOTATIONS:UNIREF90`),
    `ANNOTATIONS:UNIREF90` != "--",
    `ANNOTATIONS:UNIREF90` != ""
  )

############################################################
## 3. 从 TAXONOMY_A / TAXONOMY_B 提取 phylum
############################################################

hgt2 <- hgt2 %>%
  mutate(
    phylum_A = str_extract(TAXONOMY_A, "p__[^|]+"),
    phylum_B = str_extract(TAXONOMY_B, "p__[^|]+")
  ) %>%
  filter(!is.na(phylum_A))

############################################################
## 4. 构建 phylum-level HGT ID
## 不判断方向，只表示 TAXONOMY_A 对应谱系中的 HGT-associated gene
############################################################

hgt2 <- hgt2 %>%
  mutate(
    phylum_hgt = paste(
      phylum_A,
      `ANNOTATIONS:UNIREF90`,
      sep = " | "
    )
  )

############################################################
## 5. 统计每个 phylum_hgt 出现于多少个时间点
############################################################

hgt_persistence <- hgt2 %>%
  group_by(phylum_hgt, phylum_A, `ANNOTATIONS:UNIREF90`) %>%
  summarise(
    n_time = n_distinct(TIME),
    time_points = paste(sort(unique(TIME)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_time))

write.table(
  hgt_persistence,
  "01_HGT_persistence_all_phylumA.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
## 6. 定义 Persistent / Transient
## Persistent: n_time >= 3
## Transient: n_time < 3
############################################################

hgt_class <- hgt_persistence %>%
  mutate(
    HGT_type = ifelse(
      n_time >= 3,
      phylum_A,
      "Transient (<3 timepoints)"
    ),
    persistence_group = ifelse(
      n_time >= 3,
      "Persistent",
      "Transient"
    )
  ) %>%
  select(
    phylum_hgt,
    phylum_A,
    `ANNOTATIONS:UNIREF90`,
    n_time,
    time_points,
    HGT_type,
    persistence_group
  )

write.table(
  hgt_class,
  "02_HGT_class_persistent_vs_transient.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
## 7. 合并分类信息回原始表
############################################################

plot_data <- hgt2 %>%
  left_join(
    hgt_class %>% select(phylum_hgt, HGT_type, persistence_group, n_time),
    by = "phylum_hgt"
  )

############################################################
## 8. 时间顺序
############################################################

time_levels <- c(
  "2020_08",
  "2020_10",
  "2020_12",
  "2021_02",
  "2021_04",
  "2021_06"
)

plot_data$TIME <- factor(plot_data$TIME, levels = time_levels)

############################################################
## 9. Barplot 数据
## Y轴为每个时间点内的比例
############################################################

bar_data <- plot_data %>%
  distinct(TIME, phylum_hgt, HGT_type) %>%
  group_by(TIME, HGT_type) %>%
  summarise(
    HGT_number = n(),
    .groups = "drop"
  ) %>%
  group_by(TIME) %>%
  mutate(
    proportion = HGT_number / sum(HGT_number)
  ) %>%
  ungroup()

write.table(
  bar_data,
  "03_barplot_data_persistent_vs_transient.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

wohai############################################################
## 10. 设置颜色
############################################################

my_colors <- c(
  "p__Euryarchaeota" = "#E64B35",
  "p__Thaumarchaeota" = "#4DBBD5",
  "p__Crenarchaeota" = "#00A087",
  "p__Asgardarchaeota" = "#3C5488",
  "p__Nanoarchaeota" = "#F39B7F",
  "p__Altiarchaeota" = "#8491B4",
  "Transient (<3 timepoints)" = "black"
)

############################################################
## 11. Proportion barplot
############################################################

p_bar <- ggplot(
  bar_data,
  aes(
    x = TIME,
    y = proportion,
    fill = HGT_type
  )
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = my_colors) +
  theme_bw() +
  labs(
    x = "Time point",
    y = "Relative proportion of HGT-associated genes",
    fill = "HGT category"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

p_bar

ggsave(
  "04_persistent_vs_transient_HGT_proportion_barplot.pdf",
  p_bar,
  width = 7,
  height = 5
)

############################################################
## 12. Alluvial plot 数据
############################################################

alluvial_data <- plot_data %>%
  distinct(TIME, phylum_hgt, HGT_type) %>%
  filter(!is.na(TIME))

write.table(
  alluvial_data,
  "05_alluvial_data_persistent_vs_transient.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
## 13. Alluvial plot
############################################################

p_alluvial <- ggplot(
  alluvial_data,
  aes(
    x = TIME,
    stratum = HGT_type,
    alluvium = phylum_hgt,
    fill = HGT_type
  )
) +
  geom_flow(alpha = 0.6) +
  geom_stratum(width = 0.3) +
  scale_fill_manual(values = my_colors) +
  theme_bw() +
  labs(
    x = "Time point",
    y = "HGT-associated genes",
    fill = "HGT category"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

p_alluvial

ggsave(
  "06_persistent_vs_transient_HGT_alluvial.pdf",
  p_alluvial,
  width = 8,
  height = 5
)



library(patchwork)

p_combined <- p_bar / p_alluvial

ggsave(
  "Persistent_vs_Transient_HGT_barplot_alluvial.pdf",
  p_combined,
  width = 8,
  height = 10
)



############################################################
## 14. 输出 summary
############################################################

summary_table <- hgt_class %>%
  group_by(HGT_type, persistence_group) %>%
  summarise(
    HGT_number = n(),
    mean_n_time = mean(n_time),
    max_n_time = max(n_time),
    .groups = "drop"
  ) %>%
  arrange(persistence_group, desc(HGT_number))

write.table(
  summary_table,
  "07_summary_persistent_vs_transient_HGT.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
