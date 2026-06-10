#NMDS
#读入相关数据
library(vegan)
library(ggplot2)
library(patchwork)
#读入古菌功能
data3 <- read.table("~/Desktop/article3/metagenome/contig/coverm_result/all_kegg_count_standardized_result.txt",header = T,row.names = 1,sep = "\t",check.names = F)
time_group_meta <- read.table("~/Desktop/article3/metagenome/contig/coverm_result/sample_information_nmds",header = T,row.names = 1,sep = "\t",check.names = F)
tdata3 <- t(data3)
df_nmds <- metaMDS(tdata3, distance = "bray",k = 2,trymax = 100)


#计算r2 

original_distance <- vegdist(tdata3, method = "bray")

nmds_distance <- dist(df_nmds$points)

r_squared <- cor(original_distance, nmds_distance)^2

print(r_squared)
#提取作图数据
df_points <- as.data.frame(df_nmds$points)
#添加samp1es变量
df_points$samples <- row.names(df_points)
#修改列名
names(df_points)[1:2] <- c('NMDS1', 'NMDS2')

#修改列名
time_group_meta$samples <- row.names(time_group_meta)
#将绘图数据和分组合并
df <- merge(df_points,time_group_meta,by="samples")
df$group <- factor(df$group)
temperature_colors <- c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C")
p1 <- ggplot(data=df,aes(x=NMDS1,y=NMDS2))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(aes(color = temprature,shape = group), size=5)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  scale_color_gradientn(
    colours = temperature_colors,
    limits = c(-1,30)
  )+
  scale_shape_manual(
    values = c(11, 17, 15, 18, 16,9)
  )+
  # X轴设置
  scale_x_continuous(
    limits = c(-0.4, 1.2),
    breaks = seq(-0.4, 1.2, by = 0.4)
  ) +
  
  # Y轴设置
  scale_y_continuous(
    limits = c(-0.25, 0.75),
    breaks = seq(-0.25, 0.75, by = 0.25)
  )


p1

ggsave("~/Desktop/article3/metagenome/contig/coverm_result/function_nmds_plot_26_1_26.pdf", plot = p1, width = 10, height = 6)





#读入古菌分类
data1 <- read.table("~/Desktop/article3/archaea-16s-6-18/ASVs_counts-sliva-raw-cuta-taxonomy.txt",header = T,row.names = 1,sep = "\t",check.names = F)
time_group <-read.table("~/Desktop/article3/archaea-16s-6-18/sample_information_nmds",header = T,row.names = 1,sep = "\t",check.names = F)

tdata1 <- t(data1)
df_nmds_asv <- metaMDS(tdata1, distance = "bray",k = 2,trymax = 100)



#计算r2 

original_distance_asv <- vegdist(tdata1, method = "bray")

nmds_distance_asv <- dist(df_nmds_asv$points)

r_squared_asv <- cor(original_distance_asv, nmds_distance_asv)^2

print(r_squared_asv)



#提取作图数据
df_points_asv <- as.data.frame(df_nmds_asv$points)
#添加samp1es变量
df_points_asv$samples <- row.names(df_points_asv)
#修改列名
names(df_points_asv)[1:2] <- c('NMDS1', 'NMDS2')
head(df_points_asv)
p <- ggplot(df_points_asv,aes(x=NMDS1, y=NMDS2))+#指定数据、X轴、Y轴
  geom_point(size=3)+#绘制点图并设定大小
  theme_bw()#主题
p

#修改列名
time_group$samples <- row.names(time_group)
#将绘图数据和分组合并
df_asv <- merge(df_points_asv,time_group,by="samples")
df_asv$group <- factor(df_asv$group)
temperature_colors <- c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C")

p2<- ggplot(data=df_asv,aes(x=NMDS1,y=NMDS2))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(aes(color = tempature,shape = group), size=5)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  scale_color_gradientn(
    colours = temperature_colors,
    limits = c(-1,30)
  )+
  scale_shape_manual(
    values = c(11, 17, 15, 18, 16,9)
  )

p2

ggsave("~/Desktop/article3/archaea-16s-6-18/nmds_plot_26_1_26.pdf", plot = p2, width = 10, height = 6)


#为了美观，合并作图

# 合并两个数据框的坐标
all_x <- c(df$NMDS1, df_asv$NMDS1)
all_y <- c(df$NMDS2, df_asv$NMDS2)

# 计算统一范围
x_range <- range(all_x)
y_range <- range(all_y)

x_range
y_range

x_limits <- c(-0.8, 1.2)
y_limits <- c(-0.5, 0.75)

x_breaks <- seq(-0.8, 1.2, by = 0.4)
y_breaks <- seq(-0.5, 0.75, by = 0.25)

p1 <- p1 +
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks)

p2 <- p2 +
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks)


p_combined <- p2+ p1 +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")

p_combined


df_nmds$stress

df_nmds_asv$stress
ggsave("~/Desktop/article3/archaea-16s-6-18/archaea_function_taxonomy_nmds_plot_26_4_18.pdf", plot = p_combined, width = 16, height = 6)

