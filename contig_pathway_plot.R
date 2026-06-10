# 加载必要的包
library(tidyverse)
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
#读入文件夹位置
setwd("~/Desktop/article3/metagenome/contig/coverm_result/")
#读入氮循环的相关途径：
#读入氮循环中不同代谢途径：
read1 <- read.table("methane_m00357",header = T,check.names = F)
read1$Gene <- rownames(read1)
data_long <- read1 %>%
  melt(id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
  

#读入相关的时间信息
samples <- read.table("sample_information",header = T)
colnames(samples) <- c("Sample","group","class")
data_m00531 <- merge(data_long,samples,by="Sample",all.x = T)

time_summary <- data_m00531 %>%
  group_by(Gene,group,class) %>%
  summarise(
    Mean_Expression = mean(Expression),
    SE_Expression = sd(Expression) / sqrt(n()),
    .groups = 'drop'
  )
ggplot(time_summary, aes(x = class, y = Mean_Expression, group = Gene, color = Gene)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  # 修改x轴标签
  scale_x_continuous(
    breaks = 1:6,
    labels = c("Aug 2020", "Oct 2020", "Dec 2020", "Feb 2021", "Apr 2021", "Apr 2021")
  )+# 添加主题
  theme_bw() +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  # 添加坐标轴标签和标题
  labs(
    x = "Time",
    y = "Mean Expression",
    title = "Gene Expression Over Time",
    color = "Gene"
  )




#根据特定的途径来研究基因的时间动态：
time_summary_special <- time_summary %>% filter(Gene %in% c("K00193", "K00194","K00197"))

ggplot(time_summary_special, aes(x = class, y = Mean_Expression, group = Gene, color = Gene)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  # 修改x轴标签
  scale_x_continuous(
    breaks = 1:6,
    labels = c("Aug 2020", "Oct 2020", "Dec 2020", "Feb 2021", "Apr 2021", "Apr 2021")
  )+# 添加主题
  theme_bw() +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  # 添加坐标轴标签和标题
  labs(
    x = "Time",
    y = "Mean Expression",
    title = "Gene Expression Over Time",
    color = "Gene"
  )

