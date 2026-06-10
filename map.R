package.list=c("geoviz","tidyverse","sf","terra","rasterVis","ggspatial",
               "rgdal","rnaturalearth","rnaturalearthdata","raster")
for (package in package.list) {
  if(!require(package,character.only=T,quietly=T)){
    install.packages(package)
    library(package,character.only = T)
  }
  
}


#即墨地图绘制
library(ggplot2)
library(tidyverse)
library(sf)
library(geoviz)
library(RColorBrewer)
library(terra)
setwd("~/Desktop/article3/")
shandong <- read_sf("shandong.json")
p1<-ggplot()+
  geom_sf(data=shandong,aes(fill=NULL))+
  labs(x=NULL,y=NULL)+
  geom_sf(data=shandong,fill="white",size=0.4,color="grey")+
  xlim(120.6,121.2)+
  ylim(36.3,36.55)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#CADEF0",color="gray71",size = 0.5),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)
  )

ggsave("jimo.pdf",p1,width = 10,height = 6)

#中国地图的绘制
china <- read_sf("china.json")
p2<-ggplot()+
  geom_sf(data=china,aes(fill=NULL))+
  labs(x=NULL,y=NULL)+
  geom_sf(data=china,fill="white",size=0.4,color="grey")+
  xlim(75,137)+
  ylim(15,55)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#CADEF0",color="gray71",size = 0.5),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)
  )

ggsave("china.pdf",p2,width = 10,height = 6)


