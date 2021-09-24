setwd('E:/JLZP_rnaseq_DE/')
#PCA TPM
pre_exp_data <- read.csv('all_sample_tpm_addSYMBOL.csv')
ids <- pre_exp_data[,1:2]
duplicate_exp_data <- aggregate(x=pre_exp_data[,1:28],by=list(pre_exp_data$SYMBOL),FUN = mean,na.rm = T)
duplicate_exp_data <- duplicate_exp_data[,-c(2:3)]
colnames(duplicate_exp_data)[1] <- "SYMBOL"
row.names(duplicate_exp_data) <- duplicate_exp_data[,1]
duplicate_exp_data <- duplicate_exp_data[,-1]
for_pca_exp_data <- log2(duplicate_exp_data+1)
for_pca_exp_data[1:4,1:4]
group_1 <- c("Oocyte","Oocyte",
             "GV","GV","MI","MI","MII",
             "Ctrl_MII","Ctrl_MII","Ctrl_MII",
             "Ctrl_MII","Ctrl_MII","Ctrl_MII",
             "Ctrl_MII","Ctrl_MII","Ctrl_MII",
             "Ctrl_GV","Ctrl_GV","Ctrl_GV","Ctrl_GV"
             ,'Ctrl_GV',"Ctrl_GV","Ctrl_GV","Ctrl_GV"
             ,"Ctrl_GV","Ctrl_GV")
group_2 <- c("JZLP","JZLP","JZLP","JZLP","JZLP","JZLP","JZLP","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl")
library("FactoMineR")
library("factoextra")
for_pca_exp_data_1 <- t(for_pca_exp_data)
pca_1 <- PCA(for_pca_exp_data_1,graph = F)
summary(pca_1)
fviz_pca_ind(pca_1,geom.ind = c('point',"text"),col.ind = group_1,
             addEllipses = T, mean.point=F ,legend.title ='JZLP vs Normal')
#heatmap 
library(pheatmap)
for_heat_map_data <- log2(duplicate_exp_data +1 )
cg =names(tail(sort(apply(for_heat_map_data,1,sd)),6000))
n = t(scale(t(for_heat_map_data[cg,])))
ac = data.frame(group=group_1)
n[n>1]=1
n[n<-3]=-3
rownames(ac) = colnames(n)
pheatmap(n,show_colnames = T  ,show_rownames = F,annotation_col = ac,fontsize = 6)
ac_2 = data.frame(group=group_2)
rownames(ac_2) = colnames(n)
pheatmap(n,show_colnames = T  ,show_rownames = F,annotation_col = ac_2,fontsize = 6)
ggsave("differ_stage_PCA3.png",dpi=1080)
