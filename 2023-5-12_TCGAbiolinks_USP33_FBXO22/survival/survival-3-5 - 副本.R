##### TCGA
###最终使用时间：3-5 19：43
library(dplyr)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(ggplot2)


########### 循环测试
library(survival)
library(survminer)

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
fig_file_T=paste('/home/wus/2023-2-23_008_TCGAbiolinks_download/survival/',projects$project_id[i],"_TCGAbiolinks.png")
fig_file_S=paste('/home/wus/2023-2-23_008_TCGAbiolinks_download/survival/',projects$project_id[i],"_survival.png")
fig_file_B=paste('/home/wus/2023-2-23_008_TCGAbiolinks_download/survival/',projects$project_id[i],"_bestcut.png")
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

samplesTP <- TCGAquery_SampleTypes(colnames(geneexp), typesample = c("TP"))
print(paste0('   length(samplesTP):',length(samplesTP)))
#####
if(length(samplesTP)==0){
  print(paste0(projects$project_id[i],'have no matched sample'))
  next;}
else{

#### 基因名转换
gene_id<-rowData(data)[rowData(data)$gene_name=="USP39",'gene_id']

####
gd<-as.data.frame(geneexp)
Pair_sample_express <- gd[gene_id,samplesTP,drop=FALSE]
##转置并加一列type
Pair_sample_express<-t(Pair_sample_express)
##样本信息
sample_type_information<-colData(data)[ samplesTP  ,c('vital_status','days_to_death','days_to_last_follow_up'),drop=FALSE]

###合并
df3<-merge(Pair_sample_express, sample_type_information, by = "row.names", all = F)

print(df3)

#############以上构建出需要的数据
df<-df3

# 8、根据ABCB1基因的表达情况进行分组，取基因的平均值，大于平均值的为H，小于为L
df$ENSG00000168883.20 <- as.numeric(df$ENSG00000168883.20)
df$exp <- ''
df[df$ENSG00000168883.20 >= median(df$ENSG00000168883.20),]$exp <- "H"
df[df$ENSG00000168883.20 < median(df$ENSG00000168883.20),]$exp <- "L"

# 9、用 TCGAanalyze_survival 函数进行生存分析
TCGAanalyze_survival(df,
                     clusterCol="exp",
                     risk.table = FALSE,
                     conf.int = FALSE,
                     color = c('#da6207','#149c74'),
                     legend = projects$project_id[i],
                     height = 8, width = 8,
                     filename = fig_file_T)

#############
df2 <- df
df2 <- df2[!is.na(df2$vital_status),]
# 10.2 将status表示患者结局，1表示删失，2表示死亡
df2[df2$vital_status=='Dead',]$vital_status <- 2
df2[df2$vital_status=='Alive',]$vital_status <- 1
df2$vital_status <- as.numeric(df2$vital_status)

df2$time <- df2$days_to_death
df2$time[which(is.na(df2$time))] <- df2$days_to_last_follow_up[which(is.na(df2$time))]


res.cut <- surv_cutpoint(df2, #数据集
                         time = "time", #生存状态
                         event = "vital_status", #生存时间
                         variables = c("ENSG00000168883.20") #需要计算的数据列名
                         )

print(summary(res.cut)) #查看数据最佳截断点及统计量
res.cat <- surv_categorize(res.cut)

fit <- survfit(Surv(time, vital_status) ~ENSG00000168883.20, data = res.cat)#拟合生存分析
#绘制生存曲线并显示P值
gg<-ggsurvplot(fit,
           data = res.cat,
           risk.table = TRUE,
           pval = T,title=projects$project_id[i])
              
png(fig_file_B,width = 2500, 
    height = 2200,res=300)
print(gg, newpage = FALSE,title=projects$project_id[i])
dev.off()
                    
}}
dev.new()