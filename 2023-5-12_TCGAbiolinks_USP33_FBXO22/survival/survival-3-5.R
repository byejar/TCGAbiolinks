##### R4.0



setwd('/home/wus/2023-2-23_008_TCGAbiolinks_download/survival/')

d<-read.table('forlist.txt',sep='\t',header=F,row.name=1,check.names=FALSE)
for_list<-rownames(d)
for_list

for (i in for_list){
try({   
gene_file=paste('TCGA_',i,'.FPKM.txt.subset.txt',sep='',collapse = "")
sample_file=paste('/home/wus/TCGA_analysis/all/survival_',i,'_survival.txt',sep='',collapse = "")
write_csv=paste(i,'_for_fig.csv',sep='',collapse = "")

df<-read.table(gene_file,sep='\t',header=1,row.name=1,check.names=FALSE)
exp<-df

####

exp_tum = exp[,str_sub(colnames(exp),14,15)<=10]

exp_tum = exp_tum[,str_sub(colnames(exp_tum),16,16) =='A']


df<-exp_tum


data<-as.data.frame(t(df))
data

summary(data)

#有可能替换的位置：筛选条件
d_up=subset(data,USP39>mean(data$USP39) )
d_down=subset(data,USP39<mean(data$USP39) )
#nrow(data_filter)
name_u_list=rownames(d_up)
name_d_list=rownames(d_down)

#字符串截取和拼接
library(stringr)
    name<-function(l){
    l2=paste(unlist(strsplit(l, "-"))[1:4],collapse='-')
    l2
    l3=substr(l2,1,nchar(l2)-1)
    return (l3)
    }

name_u_data<-lapply(name_u_list,name)
name_d_data<-lapply(name_d_list,name)


#得到转换好的名称列表，比对并取出病人信息
sd<-read.table(sample_file,sep='\t',header=1,row.name=1,check.names=FALSE)

sample_uc=as.vector(unlist(name_u_data))
sasu=sort(sample_uc)
up_sd=sd[sasu,]
nrow(up_sd)
up_sd$type='up'

sample_dc=as.vector(unlist(name_d_data))
sample_dc
sasd=sort(sample_dc)
down_sd=sd[sasd,]
nrow(down_sd)
down_sd$type='down'


bind<-rbind(up_sd,down_sd)
bind
bind$OS=bind$OS+1
bind
write.csv(bind,write_csv)
    })
        }





##可用----10-9-----OS------
rm(list=ls())
library("survival")
library("survminer")
setwd('/home/wus/TCGA_analysis/2-24_all')

d<-read.table('fig_list.txt',sep='\t',header=FALSE,row.name=1,check.names=FALSE)
fig_list<-rownames(d)
fig_list

for (i in fig_list){
    fig_name=paste('./OS/',i,'_OS.png',sep='',collapse = "")


    colon=read.csv(i)
    colon
    


    fit3 <- survfit( Surv( OS.time,OS) ~ type,
                 data = colon )
       
gg<-ggsurvplot(fit3, pval = TRUE,
           xlim=c(0,2000),
           risk.table = TRUE,
           risk.table.height = 0.5
              )

png(fig_name,width = 800, height = 1000)
  print(gg, newpage = FALSE)
  dev.off()



               }
dev.new()


######################## 

############## 使用下载好的数据


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

# 建模
#fit <- survfit(Surv(time, vital_status)~exp, data=df2) # 根据表达建模

#gg<-ggsurvplot(fit, pval = TRUE,
#           #xlim=c(0,2000),
#           risk.table = TRUE,
#           risk.table.height = 0.5
#              )

#png(fig_file_S,width = 800, height = 1000)
#  print(gg, newpage = FALSE)
#  dev.off()

############################# best cut
# 1. Determine the optimal cutpoint of variables

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
              
png(fig_file_B,width = 6, height = 8,res=300)
print(gg, newpage = FALSE,title=projects$project_id[i])
dev.off()
                    
}}
dev.new()