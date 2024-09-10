############## 使用下载好的数据
#两种方式，或者使用rda，或者query但是不down，第二种虽然可以但是还是麻烦，就选择使用rda的

load("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/TCGA-BLCA.rda")
ls()
data
assayNames(data)
geneexp <- assay(data,i = "unstranded")
samplesTP <- TCGAquery_SampleTypes(colnames(geneexp), typesample = c("TP"))
PGA5 <- geneexp['ENSG00000168883.20',samplesTP]
View(PGA5)

###转换出ID
gene_id<-rowData(data)[rowData(data)$gene_name=="USP39",'gene_id']


samplesNT <- TCGAquery_SampleTypes(colnames(geneexp), typesample = c("NT"))
PGA5 <- geneexp['ENSG00000168883.20',samplesNT]
View(PGA5)

Pair_sample = TCGAquery_MatchedCoupledSampleTypes(colnames(geneexp),c("NT","TP"))
PGA5 <- geneexp['ENSG00000168883.20',Pair_sample]
View(PGA5)
typeof(PGA5)
class(PGA5)
class(geneexp)

pd<-as.data.frame(PGA5)
class(pd)

gd<-as.data.frame(geneexp)
class(gd)
rownames(gd)
colnames(gd)
PGA5 <- gd['ENSG00000168883.20',Pair_sample]



### 只有一个样本
gd<-as.data.frame(geneexp)
PGA5 <- gd['ENSG00000168883.20','TCGA-GC-A3WC-11A-11R-A22U-07',drop=FALSE]
View(PGA5)
typeof(PGA5)
class(PGA5)
class(geneexp)
tPGA5<-t(PGA5)

pd<-as.data.frame(PGA5)
class(pd)

colnames(PGA5)
#####

####################患者信息

class(colData(data))

patients<-colData(data)[Pair_sample,c('days_to_death','days_to_last_follow_up')]
typeof(patients)

patients<-colData(data)[Pair_sample,c('barcode','days_to_death')]
typeof(patients)
pad<-as.data.frame(patients)
class(pad)



dd<-as.data.frame(colData(data))
class(dd)
colnames(dd)

patients<-dd[Pair_sample,'days_to_death']

class(patients)

patients<-dd[Pair_sample,c('days_to_last_follow_up'),drop=FALSE]
patients
class(patients)


df3<-merge(tPGA5, patients, by = "row.names", all = F)

#####1.去掉B，合并T N ，调用pairline函数画图；
###文件格式
###
###______sample_______ | ____USP39_____ | _____type_____ | ____group____
###                                           N or T



###单独测试
load("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/TCGA-CHOL.rda")
ls()
geneexp <- assay(data,i = "unstranded")
#TCGAquery_SampleTypes(barcode, typesample)
##  colnames()=  vector
#### 选择是A的冰冻样本 
barcode <- colnames(geneexp) %>% 
	as.data.frame() %>%
  filter(grepl(pattern="-[0-9]+A-",colnames(geneexp)))
#取出一列就是向量类型
b<-barcode[,1]
#### 选取样本
 samplesNT=TCGAquery_SampleTypes(b, 'NT')

Pair_sample=TCGAquery_MatchedCoupledSampleTypes(b,c("NT","TP"))

#### 基因名转换
gene_id<-rowData(data)[rowData(data)$gene_name=="USP39",'gene_id']

####
gd<-as.data.frame(geneexp)
Pair_sample_express <- gd[gene_id,Pair_sample,drop=FALSE]
##转置并加一列type
Pair_sample_express<-t(Pair_sample_express)
##样本信息
sample_type_information<-colData(data)[ Pair_sample  ,c('sample_type','patient'),drop=FALSE]

###合并
df3<-merge(Pair_sample_express, sample_type_information, by = "row.names", all = F)

#######画图
library(ggplot2)
      ggplot(df3,aes(sample_type,ENSG00000168883.20)) +
      geom_line(aes(group = patient),
              color="grey"
              )+
    geom_point(aes(color=sample_type),size=3,shape=1) +
       theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y=element_blank())+
  labs(title = i)
    write.csv(dfu,csv_file)
  ggsave(fig_file)



###### 没办法直接让线在圈外，采用图层叠加的方式，线，白色圈，空心圈

ggplot(df3,aes(sample_type,ENSG00000168883.20))+
  geom_line(aes(group = patient),color="grey")+
  geom_point(size=4,shape=16,color='white')+
  geom_point(size=4,shape=1,aes(color=sample_type))+
       theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y=element_blank())+
  labs(title = i)
    write.csv(dfu,csv_file)
  ggsave(fig_file)


################################### 循环  正式使用 （颜色 注释 线粗细都合适
projects <- getGDCprojects()

library(dplyr)
projects <- projects %>% 
as.data.frame() %>% 
select(project_id,tumor) %>% 
  filter(grepl(pattern="TCGA",project_id))
###加载变量
projects

for (i in 1:nrow(projects)) {
## 0.运行信息
  print(paste0(" number ",i,",project name: ",projects$project_id[i]))

###
fig_file=paste('/home/wus/2023-2-23_008_TCGAbiolinks_download/pairline/',projects$project_id[i],"_paired_points.png")
####### 配对数据构建
TCGA_rda<-paste0("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/",projects$project_id[i],".rda")

load(TCGA_rda)

geneexp <- assay(data,i = "tpm_unstrand")
#TCGAquery_SampleTypes(barcode, typesample)
##  colnames()=  vector
#### 选择是A的冰冻样本 
barcode <- colnames(geneexp) %>% 
  as.data.frame() %>%
  filter(grepl(pattern="-[0-9]+A-",colnames(geneexp)))
#取出一列就是向量类型
b<-barcode[,1]
#### 选取样本
#samplesNT=TCGAquery_SampleTypes(b, 'NT')

Pair_sample=TCGAquery_MatchedCoupledSampleTypes(b,c("NT","TP"))
print(paste0('TCGAquery_MatchedCoupledSampleTypes:',Pair_sample,'   length(Pair_sample):',length(Pair_sample)))
#####
if(Pair_sample[1]=='Error message: there exist no matched samples'){
  print(paste0(projects$project_id[i],'have no matched sample'))
  next;}
else{
  print(Pair_sample)

#### 基因名转换
gene_id<-rowData(data)[rowData(data)$gene_name=="USP39",'gene_id']

####
gd<-as.data.frame(geneexp)
Pair_sample_express <- gd[gene_id,Pair_sample,drop=FALSE]
##转置并加一列type
Pair_sample_express<-t(Pair_sample_express)
##样本信息
sample_type_information<-colData(data)[ Pair_sample  ,c('sample_type','patient'),drop=FALSE]

###合并
df3<-merge(Pair_sample_express, sample_type_information, by = "row.names", all = F)

#### 添加NT顺序
df3$sample_type <-as.factor(df3$sample_type)
df3$sample_type <- factor(df3$sample_type, levels=rev(levels(df3$sample_type)))

###### 没办法直接让线在圈外，采用图层叠加的方式，线，白色圈，空心圈
## T: #9a0000  N: #0b5d0a  #line: #6c6c6d

ggplot(df3,aes(sample_type,ENSG00000168883.20))+
  geom_line(aes(group = patient),color="#6c6c6d",size=1)+
  geom_point(size=4,shape=16,color='white')+
  geom_point(size=5,shape=1,aes(color=sample_type),stroke = 2)+
    scale_color_manual(values = c("#9a0000", "#0b5d0a"))+ #加粗边线，修改颜色
       theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y=element_blank(),
    text = element_text(size = 20))+#改变全部字体大小)
    ylab("USP39_TPM")+
  #x轴字体大小加大
  #y轴用自己写的标注
  #图注用自己写的
  labs(title = paste0(projects$project_id[i],'  (N=T=',length(Pair_sample)/2,')'))
    #write.csv(dfu,csv_file)
  ggsave(fig_file)
}
}

################################### 循环  正式使用 （颜色 注释 线粗细都合适） 保存两个版本的图片
projects <- getGDCprojects()

library(dplyr)
projects <- projects %>% 
as.data.frame() %>% 
select(project_id,tumor) %>% 
  filter(grepl(pattern="TCGA",project_id))
###加载变量
projects

for (i in 1:nrow(projects)) {
## 0.运行信息
  print(paste0(" number ",i,",project name: ",projects$project_id[i]))

###
fig_file=paste('/home/wus/2023-2-23_008_TCGAbiolinks_download/pairline/',projects$project_id[i],"_paired_points-NT.png")
###
fig_file_txt=paste('/home/wus/2023-2-23_008_TCGAbiolinks_download/pairline/',projects$project_id[i],"_paired_points_txt-NT.png")
####### 配对数据构建
TCGA_rda<-paste0("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/",projects$project_id[i],".rda")

load(TCGA_rda)

geneexp <- assay(data,i = "tpm_unstrand")
#TCGAquery_SampleTypes(barcode, typesample)
##  colnames()=  vector
#### 选择是A的冰冻样本 
barcode <- colnames(geneexp) %>% 
  as.data.frame() %>%
  filter(grepl(pattern="-[0-9]+A-",colnames(geneexp)))
#取出一列就是向量类型
b<-barcode[,1]
#### 选取样本
#samplesNT=TCGAquery_SampleTypes(b, 'NT')

Pair_sample=TCGAquery_MatchedCoupledSampleTypes(b,c("NT","TP"))
print(paste0('TCGAquery_MatchedCoupledSampleTypes:',Pair_sample,'   length(Pair_sample):',length(Pair_sample)))
#####
if(Pair_sample[1]=='Error message: there exist no matched samples'){
  print(paste0(projects$project_id[i],'have no matched sample'))
  next;}
else{
  print(Pair_sample)

#### 基因名转换
gene_id<-rowData(data)[rowData(data)$gene_name=="USP39",'gene_id']

####
gd<-as.data.frame(geneexp)
Pair_sample_express <- gd[gene_id,Pair_sample,drop=FALSE]
##转置并加一列type
Pair_sample_express<-t(Pair_sample_express)
##样本信息
sample_type_information<-colData(data)[ Pair_sample  ,c('sample_type','patient'),drop=FALSE]

###合并
df3<-merge(Pair_sample_express, sample_type_information, by = "row.names", all = F)
#### 添加NT顺序
df3$sample_type <-as.factor(df3$sample_type)
df3$sample_type <- factor(df3$sample_type, levels=rev(levels(df3$sample_type)))

###### 没办法直接让线在圈外，采用图层叠加的方式，线，白色圈，空心圈
## T: #9a0000  N: #0b5d0a  #line: #6c6c6d

ggplot(df3,aes(sample_type,ENSG00000168883.20))+
  geom_line(aes(group = patient),color="#6c6c6d",size=1)+
  geom_point(size=4,shape=16,color='white')+
  geom_point(size=5,shape=1,aes(color=sample_type),stroke = 2)+
    scale_color_manual(values = c("#9a0000", "#0b5d0a"))+ #加粗边线，修改颜色
       theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y=element_blank(),
    text = element_text(size = 20))+#改变全部字体大小)
    ylab("USP39_TPM")+
  #x轴字体大小加大
  #y轴用自己写的标注
  #图注用自己写的
  labs(title = paste0(projects$project_id[i],'  (N=T=',length(Pair_sample)/2,')'))
    #write.csv(dfu,csv_file)
  ggsave(fig_file,dpi=300,width=7,height=6)
}



ggplot(df3,aes(sample_type,ENSG00000168883.20))+
  geom_line(aes(group = patient),color="#6c6c6d",size=1)+
  geom_point(size=4,shape=16,color='white')+
  geom_point(size=5,shape=1,aes(color=sample_type),stroke = 2)+
    scale_color_manual(values = c("#9a0000", "#0b5d0a"))+ #加粗边线，修改颜色
       theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    text = element_text(size = 20))+#改变全部字体大小)
    ylab("USP39_TPM")+
  #x轴字体大小加大
  #y轴用自己写的标注
  #图注用自己写的
  labs(title = paste0(projects$project_id[i],'  (N=T=',length(Pair_sample)/2,')'))
  ggsave(fig_file_txt,dpi=300,width=7.5,height=7)
}
















#####2.