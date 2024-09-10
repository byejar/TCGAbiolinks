library(ggpubr)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
#library(ggstatsplot)



################# 单独测试
setwd('/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data')
load("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/TCGA-PAAD.rda")
ls()
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
 samples=TCGAquery_SampleTypes(b, c("TP"))

#R=TCGAquery_MatchedCoupledSampleTypes(b,c("NT","TP"))

#### 基因名转换
EI_gene_list<-vector()
EI_gene=c('ANGPT2','CDH5','ESAM','ESM1','ERG','ICAM2','TIE1')
for (i in EI_gene){
  print(i)
  gene_id<-rowData(data)[rowData(data)$gene_name==i,'gene_id']
  print(gene_id)
  EI_gene_list<-append(EI_gene_list,gene_id)
}
EI_gene_list
####
gd<-as.data.frame(geneexp)
samples_express <- gd[EI_gene_list,samples,drop=FALSE]
##转置并加一列type
samples_express<-t(samples_express)
samples_express
colnames(samples_express) <- EI_gene

csv_file='TCGA-PAAD.csv'
write.csv( samples_express, csv_file,row.names=T)

##########USP33 USP39 FBXO22 DDIT4
cor_gene_list<-vector()
cor_gene=c('USP33','USP39','FBXO22','DDIT4')
for (i in cor_gene){
  print(i)
  gene_id<-rowData(data)[rowData(data)$gene_name==i,'gene_id']
  print(gene_id)
  cor_gene_list<-append(cor_gene_list,gene_id)
}
cor_gene_list
####
gd<-as.data.frame(geneexp)
samples_express <- gd[cor_gene_list,samples,drop=FALSE]
##转置并加一列type
samples_express<-t(samples_express)
samples_express[1:4,1:4]
colnames(samples_express) <- cor_gene

csv_file='TCGA-PAAD-cor_gene.csv'
write.csv( samples_express, csv_file,row.names=T)



./calculate_ei -o /home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/output_ei/TCGA-PAAD-ei.csv  /home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data/TCGA-PAAD.csv


a<-read.csv('/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/output_ei/TCGA-PAAD-ei.csv.csv')
b<-read.csv('/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data/TCGA-PAAD-cor_gene.csv')
a[1:4,1:4]
b[1:4,1:4]

m=merge(a,b,by.x='Unnamed..0',by.y='X')
m
write.csv(m,'/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data/TCGA-PAAD-merge.csv')
################
mtcars<-read.csv('TCGA-PAAD-merge.csv')
x <- mtcars$FBXO22
y <- mtcars$EI_Score
plot(x, y, main = "CHOL",
     xlab = "FBXO22", ylab = "EI_Score",
     pch = 19, frame = T)
# 添加回归线
mtcars2<-mtcars[,grepl('FBXO22 | EI_Score', colnames(mtcars))] #提取列
abline(lm(y ~ x, data = mtcars2), col = "red")


########## scRNA_v2
library(ggstatsplot)
setwd('/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data')
mtcars<-read.csv('TCGA-PAAD-merge.csv')
ggscatterstats(mtcars, 
               y =DDIT4, 
               x =EI.Score,
               type = "pearson",
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#009E73",
               yfill = "#D55E00", 
               marginal.type = "histogram",  
               title = "PAAD Relationship between DDIT4 and EI_Score")



###############循环
setwd('/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data')

list=c('TCGA-COAD','TCGA-READ','TCGA-KIRC','TCGA-UCEC')
for (j in list){

rda_file=paste("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/",j,".rda",sep='')
load(rda_file)
ls()
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
 samples=TCGAquery_SampleTypes(b, c("TP"))

#R=TCGAquery_MatchedCoupledSampleTypes(b,c("NT","TP"))

#### 基因名转换
EI_gene_list<-vector()
EI_gene=c('ANGPT2','CDH5','ESAM','ESM1','ERG','ICAM2','TIE1')
for (i in EI_gene){
  print(i)
  gene_id<-rowData(data)[rowData(data)$gene_name==i,'gene_id']
  print(gene_id)
  EI_gene_list<-append(EI_gene_list,gene_id)
}
EI_gene_list
####
gd<-as.data.frame(geneexp)
samples_express <- gd[EI_gene_list,samples,drop=FALSE]
##转置并加一列type
samples_express<-t(samples_express)
samples_express
colnames(samples_express) <- EI_gene

#csv_file='TCGA-PAAD.csv'
csv_file=paste(j,'.csv',sep='')
write.csv( samples_express, csv_file,row.names=T)

##########USP33 USP39 FBXO22 DDIT4
cor_gene_list<-vector()
cor_gene=c('USP33','USP39','FBXO22','DDIT4','RGS5')
for (i in cor_gene){
  print(i)
  gene_id<-rowData(data)[rowData(data)$gene_name==i,'gene_id']
  print(gene_id)
  cor_gene_list<-append(cor_gene_list,gene_id)
}
cor_gene_list
####
gd<-as.data.frame(geneexp)
samples_express <- gd[cor_gene_list,samples,drop=FALSE]
##转置并加一列type
samples_express<-t(samples_express)
samples_express[1:4,1:4]
colnames(samples_express) <- cor_gene

#csv_file='TCGA-PAAD-cor_gene.csv'
csv_file=paste(j,'-cor_gene.csv',sep='')
write.csv( samples_express, csv_file,row.names=T)
}


###########

./calculate_ei -o /home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/output_ei/TCGA-KIRC-ei.csv  /home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data/TCGA-KIRC.csv


####

list=c('TCGA-COAD','TCGA-READ','TCGA-KIRC')
for (j in list){

a_file=paste('/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/output_ei/',j,'-ei.csv.csv',sep='')
b_file=paste('/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data/',j,'-cor_gene.csv',sep='')


a<-read.csv(a_file)
b<-read.csv(b_file)
a[1:4,1:4]
b[1:4,1:4]

m=merge(a,b,by.x='Unnamed..0',by.y='X')
m
csv=paste('/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data/',j,'-merge.csv',sep='')
write.csv(m,csv)
}
################


########## scRNA_v2
library(ggstatsplot)
library(ggplot2)
setwd('/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data')
mtcars<-read.csv('/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data/ TCGA-KIRC -merge.csv')


ggscatterstats(mtcars, 
               y =USP39,
               x =EI.Score,
               type = "pearson",
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#009E73",
               yfill = "#D55E00", 
               marginal.type = "histogram",  
               title = "PAAD Relationship between USP39 and EI_Score")

ggsave('PAAD_USP39.png')

ggscatterstats(mtcars, 
               y =USP33,
               x =EI.Score,
               type = "pearson",
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#009E73",
               yfill = "#D55E00", 
               marginal.type = "histogram",  
               title = "PAAD Relationship between USP33 and EI_Score")

ggsave('PAAD_USP33.png')

ggscatterstats(mtcars, 
               y =FBXO22,
               x =EI.Score,
               type = "pearson",
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#009E73",
               yfill = "#D55E00", 
               marginal.type = "histogram",  
               title = "PAAD Relationship between FBXO22 and EI_Score")

ggsave('PAAD_FBXO22.png')

ggscatterstats(mtcars, 
               y =DDIT4,
               x =EI.Score,
               type = "pearson",
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#009E73",
               yfill = "#D55E00", 
               marginal.type = "histogram",  
               title = "PAAD Relationship between DDIT4 and EI_Score")

ggsave('PAAD_DDIT4.png')



####################
library(ggplot2)

# 基本箱线图

p = ggplot(mtcars,  aes(y =EI.Score,x=x)) + 
  geom_boxplot()
p


mtcars<-read.csv('/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data/KIRC.csv')
p = ggplot(mtcars,  aes(y =Y,x=X)) + 
  geom_boxplot()
p


############


mtcars<-read.csv('/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data/KIRC-TEST-v2.csv')



ggscatterstats(mtcars, 
               y =USP39,
               x =EI,
               type = "pearson",
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#009E73",
               yfill = "#D55E00", 
               marginal.type = "histogram",  
               title = "KIRC Relationship between USP39 and EI_Score")


mtcars<-read.csv('/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data/KIRC-TEST.csv')
log2(mtcars$EI)

ks.test(mtcars$EI,"pnorm",mean=mean(mtcars$EI),sd=sqrt(var(mtcars$EI)))

shapiro.test(mtcars$EI)



##############
##  complete
library("ggplot2")
data = read.csv('/home/wus/2023-5-23_014_TCGA-Endothelial_Index_score/calculate_ei_linux/TCGA_sample_data/KIRC-TEST-v2.csv')
## wilcox
#bgi = data[data$group=="BGISEQ",]$Completeness
# illu = data[data$group=="Illumina",]$Completeness
# mean(bgi)
# mean(illu)
# wilcox.test(bgi, illu, paired = TRUE)


library(ggpubr)
#before <-c(200.1, 190.9, 192.7, 213, 241.4, 196.9, 172.2, 185.5, 205.2, 193.7)
#after <-c(392.9, 393.2, 345.1, 393, 434, 427.9, 422, 383.9, 392.3, 352.2)
#d <- data.frame(before = before, after = after)
ggpaired(data, cond1 = "USP39", cond2 = "EI",
fill = "condition", palette = "jco")



#########
#step1 read packages
install.packages("MASS")  
library(MASS)   #基于此包进行box-cox转换
install.packages("moments")
library(moments) #基于此包进行峰度偏度计算
#step2 构建线性模型并检查是否满足正态分布
lm.sales1<-lm(yi~xi,data=data)
res2<-lm.sales1$residual 
hist(res2)  #绘制残差频率分布图
#head(res2)
skewness(res2)  #计算偏度，结果：0.13
#通过残差频率分布图和偏度值可以发现残差并不满足正态分布
#step3 进行box-cox转换，确定λ值
b2=boxcox(yi~xi, data=mydata2) # 定义函数类型和数据
I=which(b2$y==max(b2$y))
b2$x[I]#lambda=0.55,即为图中最高点
