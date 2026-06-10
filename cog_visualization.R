library(reshape2)
library(dplyr)
library(tidyr)
library(circlize)
library(RColorBrewer)


#COG  数量
library(ggplot2)
library(scales)
library(cowplot)

#hgtector 可视化分析
#Thaumarchaeota
setwd("/Users/hahamark/Desktop/article3/metagenome/hgt/hgtector")
data2 <- read.table('Thaumarchaeota_cog_count_raw',header = T,sep = " ")
data3 <- read.table('Euryarchaeota_cog_count_raw',header = T,sep = " ")
hgtector_all<- rbind(data2,data3)



#绘制contig_level的图片
contig_all <- read.table('/Users/hahamark/Desktop/article3/contig_hgt/all_time_hgt_COG',header = T,sep = "\t")
#修改做图代码，修改为散点图

hgtector_all$group <- "hgtector"
contig_all$group <- "waafle"
hgt_all<- rbind(hgtector_all,contig_all)

#设置渐变性颜色：
color1 <- brewer.pal(10,"Paired")
color2 <- colorRampPalette(color1)(24)
# 绘制散点图
p1 <- ggplot(hgt_all, aes(x = phylum, y = type, size = value, color = type)) +
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1)) +  # 添加抖动避免点重叠
  scale_size_continuous(range = c(3, 18))+
  labs(x = "Category", y = "Value", size = "Size", color = "Color") +  # 设置轴标签和图例标题
  facet_grid(~group,scales = "free_x")+
  scale_color_manual(values=color2)


ggsave(filename = "/Users/hahamark/Desktop/article3/contig_hgt/cog_visualization.pdf",width = 6,height = 10)
