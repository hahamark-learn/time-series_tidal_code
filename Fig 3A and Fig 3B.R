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






#metachip_archaea_bacteria
data5 <- read.table("/Users/hahamark/Desktop/article3/metagenome/hgt/hgt_result/Euryarchaeota_cog_count_raw",sep = " ",header = T)
#metachip_archaea_to_archaea
data6 <- read.table("/Users/hahamark/Desktop/article3/metagenome/hgt/hgt_result/Thermoplasmatota_cog_count_raw",sep = " ",header = T)


#绘制分面图
data2$phylum <- "Thaumarchaeota"
data3$phylum <- "Euryarchaeota"

hgtector_all<- rbind(data2,data3)



#绘制metachip的图片

data5$phylum <- "Euryarchaeota"
data6$phylum <- "Thermoplasmatota"
metachip_all<- rbind(data5,data6)

#修改做图代码，修改为散点图

hgtector_all$group <- "hgtector"
metachip_all$group <- "metachip"
hgt_all<- rbind(hgtector_all,metachip_all)

#设置渐变性颜色：
color1 <- brewer.pal(10,"Paired")
color2 <- colorRampPalette(color1)(20)
# 绘制散点图
p1 <- ggplot(hgt_all, aes(x = phylum, y = type, size = value, color = type)) +
  geom_point(alpha = 0.7, position = position_jitter(width = 0.15)) +  # 添加抖动避免点重叠
  scale_size_continuous(range = c(3, 18))+
  labs(x = "Category", y = "Value", size = "Size", color = "Color") +  # 设置轴标签和图例标题
  facet_grid(~group,scales = "free_x")+
  scale_color_manual(values=color2)

ggsave("/Users/hahamark/Desktop/article3/metagenome/hgt/hgt_plot_25.6.17.pdf", plot = p1, width = 210 / 25.4, height = 297 / 25.4)
