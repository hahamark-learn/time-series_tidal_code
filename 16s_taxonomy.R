#16s柱状堆积图
library(tidyverse)
library(ggprism)
library(vegan)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(tibble)
#读入数据，这里要排除
setwd("~/Desktop/article3/archaea-16s-6-18/")
#读入原始的count数据
asv_raw <- read.table("ASVs_counts-sliva-raw-cuta.txt",header = T,row.names = 1,sep = "\t",check.names = F)
taxa_archaea <- read.table("archaea_taxonomy_archaea",header = T, row.names = 1,sep = "\t",check.names = F)
#先进行相关抽平的工作
colSums(asv_raw)
#使用该代码抽平
asv_Flattening = as.data.frame(t(rrarefy(t(asv_raw), min(colSums(asv_raw)))))
tasv_Flattening <- t(asv_Flattening)
#查看抽平后的每个样本的和
colSums(asv_Flattening)
#将抽平后的otu表保存到该工作目录，准备后面的多样性分析
#排出按照行进行累加为0的行。
row_sums <- rowSums(asv_Flattening)
rows_to_remove <- which(row_sums==0)
asv_clean_flattening <- asv_Flattening[-rows_to_remove,]
write.table(asv_clean_flattening,"ASVs_counts_clean_Flattening-sliva.txt",quote=F,sep = "\t")



#提取相关的archaea的count数据
asv_count_archaea_taxa <- merge(taxa_archaea,asv_clean_flattening,by="row.names")
asv_count_archaea <- asv_count_archaea_taxa[,-c(2:7)]
asv_count_archaea <- column_to_rownames(asv_count_archaea,var="Row.names")
#写入只含有archaea count数量
write.table(file = "ASVs_counts-sliva-raw-cuta-taxonomy_Flattening.txt",asv_count_archaea,quote = F,sep = "\t")
group<- read.table("sample_information",header = T,row.names = 1,check.names = F)
dat <- merge(x=asv_count_archaea,y=taxa_archaea,by='row.names')
dat=dplyr::rename(dat,OTUID=Row.names)
num_a <- ncol(asv_count_archaea)+1
aa<-aggregate(dat[,2:num_a],by=list(dat$Phylum),FUN=sum)
row.names(aa)=aa$Group.1   
head(aa)
aa<-dplyr::select(aa,-Group.1)
#根据行求和结果对数据排序
order<-sort(rowSums(aa[,1:ncol(aa)]),index.return=TRUE,decreasing=T)   
#根据列求和结果对表格排序
cc<-aa[order$ix,]


##只展示排名前10的物种，之后的算作Others(根据需求改数字),发现古菌只有6个
##dd<-rbind(colSums(cc[11:as.numeric(length(rownames(cc))),]),cc[10:1,])
##head(dd, n = 3)
##rownames(dd)[1]<-"Others"
##head(dd, n = 3)
#
bb<-merge(t(cc),group,by = "row.names")
kk<-tidyr::gather(bb,Phylum,Abundance,-c(group,Row.names,class))
hh <- kk %>%group_by(group,Phylum) %>%dplyr :: summarise(Abundance=sum(Abundance))

color1 <- rev(c("#BC80BD","#CCEBC5","#B3DE69","#BEBADA","#DADF00","#8DD3C7","#FDB462","#FCCDE5","#D9D9D9"))
hh$Phylum = factor(hh$Phylum,order = T,levels = rev(c("Thaumarchaeota","Euryarchaeota","Crenarchaeota","Asgardaeota","Nanoarchaeaeota","Altiarchaeota","Hydrothermarchaeota","Diapherotrites","Hadesarchaeaeota")))
#hh$Location <-factor(hh$Location,levels= c('Lianyungang','Zhuhai','Sanya','Dandong','Dongying','Qingdao','Xiamen','Shantou','Ningbo','Wenzhou','Beihai','Yancheng'))
#c('Lianyungang','Zhuhai','Sanya','Dandong','Dongying','Qingdao','Xiamen','Shantou','Ningbo','Wenzhou','Beihai','Yancheng')
hh$group <-factor(hh$group,levels=c("2008" ,"2010","2012","2102","2104","2106"))
#绘制柱状堆积图

p1<- ggplot(hh)+
  geom_bar(aes(x =group,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.93)+
  scale_fill_manual(values = color1) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 90))
p1

ggsave("time_16s_archaea_taxonomy.pdf",p1,width = 8,height = 6)


#例子1-1as
plot_16s <- ggplot(hh3)+
  geom_bar(aes(x =id,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.9)+
  scale_fill_manual(values = rev(color1)) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))+
  guides(fill = FALSE)

#提取为古菌的部分进行 alpha 研究
asv2 <- merge(x=asv,y=taxa1,by='row.names')
row.names(asv2) <- asv2$Row.names
#去除读取的列
asv3 <- asv2[,-c(1,86,87,88,89,90,91)]
#生产新的asv表
asv<- asv3
#数据抽平
colSums(asv)
#使用该代码抽平
asv_Flattening = as.data.frame(t(rrarefy(t(asv), min(colSums(asv)))))

tasv_Flattening <- t(asv_Flattening)
#查看抽平后的每个样本的和
colSums(asv_Flattening)
#将抽平后的otu表保存到该工作目录，准备后面的多样性分析
asv <- asv_Flattening
list2 <- c(rowSums(asv)>0)
list3 <- data.frame(list2)
#通过统计每一行之和，来取不为0的那一行。
asv1 <- asv[rowSums(asv)>0,,drop=F]
write.table(asv1,"clean_flatting_sum_asv",sep = "\t",quote=F)
#排除到综合为0的行之后，计算相关的多样性
richness1<- specnumber(asv1,MARGIN = 2)
data1 <- data.frame(richness1)

#合并列名，进行统计
rownames(group)<-group$id

data2 <- merge(x=data1,y=group,by='row.names')
#按照时间线为x轴，进行相关性的拟合

#重新写一下 时间的相关脚本


#这里按照季节进行分类：
#环境因子和古菌otu的相关性
#读入otu表
#读入环参


#这里按照季节进行相关的分析：
#Winter: January–March; Spring: April–June; Summer: July–September; Fall: October–December) 
#读入数据
#16s柱状堆积图
library(tidyverse)
library(ggprism)
library(vegan)
library(RColorBrewer)
library(reshape2)
library(dplyr)
#读入数据，这里要排除
setwd("~/Desktop/article3/16s/")
asv<- read.table("ASVs_counts_flatting.txt",sep = '\t',header = T,row.names = 1)
taxa1<- read.table("ASVs_taxonomy-rdp.txt", sep="\t",header = T,row.names = 1)
group<- readxl::read_xlsx("sample_df_1.xlsx")
rownames(group) <- group$id 
dat <- merge(x=asv,y=taxa1,by='row.names')
dat=dplyr::rename(dat,OTUID=Row.names)
num_a <- ncol(asv)+1
aa<-aggregate(dat[,2:num_a],by=list(dat$Phylum),FUN=sum)
row.names(aa)=aa$Group.1   
head(aa)
aa<-dplyr::select(aa,-Group.1)
head(aa, n = 3)
#根据行求和结果对数据排序
order<-sort(rowSums(aa[,2:ncol(aa)]),index.return=TRUE,decreasing=T)   
#根据列求和结果对表格排序
cc<-aa[order$ix,]
head(cc, n = 10)
cc1 <- cc[1:6,]

##只展示排名前10的物种，之后的算作Others(根据需求改数字),发现古菌只有6个
##dd<-rbind(colSums(cc[11:as.numeric(length(rownames(cc))),]),cc[10:1,])
##head(dd, n = 3)
##rownames(dd)[1]<-"Others"
##head(dd, n = 3)
#
bb<-merge(t(cc1),dplyr::select(group,id,season),
          by.x = "row.names",by.y ="id")
kk<-tidyr::gather(bb,Phylum,Abundance,-c(season,Row.names))
kk$Phylum<-ifelse(kk$Phylum=='Unassigned','Others',kk$Phylum)
hh <- kk %>%group_by(season,Phylum) %>%dplyr :: summarise(Abundance=sum(Abundance))

color1 <- brewer.pal(12,"Set3")
color1 <- rev(c("#BC80BD","#CCEBC5","#B3DE69","#BEBADA","#DADF00","#8DD3C7"))
hh$Phylum = factor(hh$Phylum,order = T,levels = rev(row.names(cc1)))
#hh$Location <-factor(hh$Location,levels= c('Lianyungang','Zhuhai','Sanya','Dandong','Dongying','Qingdao','Xiamen','Shantou','Ningbo','Wenzhou','Beihai','Yancheng'))
#c('Lianyungang','Zhuhai','Sanya','Dandong','Dongying','Qingdao','Xiamen','Shantou','Ningbo','Wenzhou','Beihai','Yancheng')
hh$season <-factor(hh$season,levels=c("spring","summer","fall","winter"))


#color2 <- c("#FDB462" ,"#FCCDE5", "#D9D9D9" ,"#B3DE69", "#CCEBC5","#BC80BD" )


hh2 <- data.frame(hh1,id1)
colnames(hh2)<-(c("Location","id"))
hh3 <- merge(hh,hh2,by="Location")

#绘制柱状堆积图

p1<- ggplot(hh)+
  geom_bar(aes(x =season,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.93)+
  scale_fill_manual(values = color1) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 90))

ggsave("time_16s_season_taxonomy.pdf",p1,width = 8,height = 6)






