# jlzp de gene
# DEseq2
# according to PCA GV vs GV ; MII vs MII
setwd("E:/JLZP_rnaseq_DE/")
library(DESeq2)
library(clusterProfiler)
options(stringsAsFactors = F)
rm(list = ls())
JLZP_counts <- read_csv("gene_counts_all.csv")
ctrl_counts <- read_csv('ctrl_MII_GV.csv')
ids <- JLZP_counts[,c(1:2)]
pre_JZLP <- JLZP_counts[,-2]
duplicated_pre_JLZP <- aggregate(x=pre_JZLP,by=list(pre_JZLP$gene_id),FUN = mean,
                                 na.rm = T)
duplicated_pre_JLZP <- duplicated_pre_JLZP[,-2]
rownames(duplicated_pre_JLZP) <- duplicated_pre_JLZP[,1]
duplicated_pre_JLZP <- duplicated_pre_JLZP[,-1]
colnames(duplicated_pre_JLZP) <- c("Oocyte_1","Oocyte_2",
                              "GV_1","GV_2","MI_1","MI_2","MII_1")
keep_feature_JLZP <- rowSums(duplicated_pre_JLZP > 1) >=  1
table(keep_feature_JLZP)
keep_gene_jzlp <- duplicated_pre_JLZP[keep_feature_JLZP,]
write.csv(keep_gene_jzlp,"pre_clean_jlzp.csv")
for_use_jlzp <- read.csv('pre_clean_jlzp.csv')
colnames(for_use_jlzp)[1] <- c("gene_id")
merge_data <- merge(for_use_jlzp,ids,by = 'gene_id')
write.csv(merge_data,'1stclean_JLZP_counts.csv')
colnames(ctrl_counts) <- c('gene_id',"Ctrl_MII_1","Ctrl_MII_2","Ctrl_MII_3",
                           "Ctrl_MII_4","Ctrl_MII_5","Ctrl_MII_6",
                           "Ctrl_MII_7","Ctrl_MII_8","Ctrl_MII_9",
                           "Ctrl_GV_1","Ctrl_GV_2","Ctrl_GV_3","Ctrl_GV_4"
                           ,'Ctrl_GV_5',"Ctrl_GV_6","Ctrl_GV_7","Ctrl_GV_8"
                           ,"Ctrl_GV_9","Ctrl_GV_10")
write.csv(ctrl_counts,'1st_clean_ctrl_counts.csv')
rm(list=ls())

ctrl_counts_for <- read_csv('1st_clean_ctrl_counts.csv')
JLZP_counts_for <- read_csv('1stclean_JLZP_counts.csv')
ctrl_counts_for <- ctrl_counts_for[,-1]
JLZP_counts_for <- JLZP_counts_for[,-1]
GV_DE_data <- merge(JLZP_counts_for[,c(1,4:6,9)],ctrl_counts_for[,c(1,11:20)],by= "gene_id")
MII_DE_data <- merge(JLZP_counts_for[,c(1:3,9)],ctrl_counts_for[,c(1:10)],by = "gene_id")
GV_DE_data <- GV_group_DE_data[,c(1,5,2:4,6:15)]
MII_DE_DATA <- MII_DE_data[,c(1,4,2:3,5:ncol(MII_DE_data))]
gv_ids <-  GV_DE_data[,c(1:2)]
MII_ids <- MII_DE_DATA[,c(1:2)]

#pre-clean counts data gene counts > 1
DE_GV_list <- GV_DE_data
duplicate_GV_List <- aggregate(x = DE_GV_list,by = list(DE_GV_list$Symbol), FUN = mean, na.rm = T)
DE_GV <- duplicate_GV_List[,-c(2:3)]
rownames(DE_GV) <- DE_GV[,1]
keep_GV_ft <- rowSums(DE_GV > 1) > 1 

table(keep_GV_ft)

DE_GV <- DE_GV[keep_GV_ft,]
DE_GV_1 <- DE_GV[,-1]
Group_GV <- c(rep('case',3),rep('ctrl',10))
Group_GV <- factor(Group_GV,levels = c("case","ctrl"))
library(DESeq2)
condition <- Group_GV 
coldata <- data.frame(row.names = colnames(DE_GV_1),condition)
DE_GV_2 <- apply(DE_GV_1,2,as.integer)
row.names(DE_GV_2) <- row.names(DE_GV_1)
dds_GV <- DESeqDataSetFromMatrix(countData = DE_GV_2,colData = coldata,design = ~condition)
dds_GV$condition <- relevel(dds_GV$condition,ref='ctrl')
dds_GV_1 <- DESeq(dds_GV)
nrDEG_DESeq2_GV <- as.data.frame(results(dds_GV_1))
nrDEG_DESeq2_GV <- nrDEG_DESeq2_GV[order(nrDEG_DESeq2_GV$log2FoldChange),]
# sample clust
library(factoextra)
vsd <- vst(dds_GV_1,blind = T)
sampleDists <- dist(t(assay(vsd)))
res1 <- hcut(sampleDists, k = 2, stand = FALSE,hc_method ="average" ) 
fviz_dend(res1,
          rect_fill = T,
          # 字体大小
          cex = 1,
          # 字体颜色
          color_labels_by_k=T,
          # 平行放置
          horiz=T)
resultsNames(dds_GV_1)
res_GV <- results(dds_GV_1)
summary(res_GV)
table(res_GV$padj < 0.05)
res_GV <- res_GV[order(res_GV$padj < 0.05),]
GV_deg <- nrDEG_DESeq2_GV
GV_deg$Group ="N.S"
logFC_cutoff <- 1
GV_deg$Group[which((GV_deg$pvalue < 0.05) & (GV_deg$log2FoldChange > logFC_cutoff))] ="Up_regulated"
GV_deg$Group[which((GV_deg$pvalue < 0.05) & (GV_deg$log2FoldChange < -logFC_cutoff))] ="Down_regulated"
table(GV_deg$Group)
write.csv(GV_deg,'GV_DEG.csv')
# DE gene valconal demo('colors')
library(ggplot2)
GV_deg$color <- ifelse(GV_deg$pvalue<0.05 & abs(GV_deg$log2FoldChange)>= 1,ifelse(GV_deg$log2FoldChange > 1,'red1','blue2'),'gray37')
GV_deg$color[which((GV_deg$pvalue<0.05) & (GV_deg$log2FoldChange > 2))] ="darkorange"
GV_deg$color[which((GV_deg$pvalue <0.05) & (GV_deg$log2FoldChange < -2))] = "deepskyblue3"
color <- c(red1 = "red1",gray37 = "gray37",blue2 = "blue2",darkorange = "darkorange",deepskyblue3="deepskyblue3")
valcanol_GV <- ggplot(GV_deg,aes(log2FoldChange,-log10(pvalue),col = color)) +
                        geom_point() +
                        theme_bw() +
                        scale_color_manual(values = color,
                                           name = "Fold Change",
                                           breaks = c('red1',"gray37","blue2","darkorange","deepskyblue3"),
                                           labels = c('2X Fold Change','N.S','-2X Fold Change','4X Fold Change','-4X Fold Change')) +
                        labs(x = "log2_FC",y = "-log10(P_value)")+
                        geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6)+
                        geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6)+
theme(panel.grid = element_blank(),axis.title = element_text(size = 16),axis.text = element_text(size = 14))
                       


valcanol_GV 
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
BiocManager::install('pathview')
library(pathview)
#GO & KEGG
GV_for_GO <- GV_deg
GV_for_GO$Symbol <- rownames(GV_for_GO)
GV_for_GO_merge <- merge(GV_for_GO,gv_ids,by = 'Symbol')
GV_GO_use <- GV_for_GO_merge[order(GV_for_GO_merge$log2FoldChange),]
sort_GV_GO <- GV_GO_use %>%
  arrange(desc(abs(log2FoldChange)))
sort_GV_GO_Use <- sort_GV_GO[,c(1,3,10)]
write.csv(sort_GV_GO_Use,'DF_byFC_for_GO.csv')
rm(list = ls())
GO_GENE <- read.csv('DF_byFC_for_GO.csv')
GO_GENE <- GO_GENE[,-1]
#GESA KEGG
genelist <- GO_GENE[,c(1:2)]
genelist <- genelist %>%
  arrange(desc(log2FoldChange))
genelist_for_GO = genelist$log2FoldChange
names(genelist_for_GO) = genelist$Symbol
head(genelist_for_GO)
gese_go <- gseGO(genelist_for_GO,keyType = 'SYMBOL',ont = "BP",OrgDb = org.Hs.eg.db,pvalueCutoff = 0.5,                 pAdjustMethod = "BH")
head(gese_go,2)
gseaplot2(gese_go,1:nrow(gese_go@result))
gesa_go_res <- gese_go@result
write.csv(gesa_go_res,'GESA_GO_RES.csv')

genelist_for_KEGG <- GO_GENE[,c(2,3)]
gene=bitr(genelist_for_KEGG$gene_id,fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
colnames(gene)[1] = 'gene_id'
kegg_gene <- merge(genelist_for_KEGG,gene,by = 'gene_id')
genelist_for_KEGG_1 <- kegg_gene %>%
  arrange(desc(log2FoldChange))
genelist_use_KEGG <-genelist_for_KEGG_1$log2FoldChange
names(genelist_use_KEGG) <- genelist_for_KEGG_1$ENTREZID

gesa_kegg <- gseKEGG(genelist_use_KEGG,organism = 'hsa',keyType = 'kegg',pvalueCutoff = 0.5)
gesa_kegg_reslut <- gesa_kegg@result
write.csv(gesa_kegg_reslut,'gesa_kegg.csv')
gseaplot2(gesa_kegg,1:nrow(gesa_kegg_reslut))

#go & kegg enrich
library(cluster)
library(org.Hs.eg.db)
library(enrichplot)
library(pathview)
rm(list = ls())
GO_GENE <- read.csv('DF_byFC_for_GO.csv')
GO_GENE <- GO_GENE[,-1]
gene <- bitr(GO_GENE$gene_id,fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb="org.Hs.eg.db")
colnames(GO_GENE)[3] <- 'ENSEMBL'  
gene_list <- merge(gene,GO_GENE,by='ENSEMBL')
genelist_1 <- gene_list[which(abs(gene_list$log2FoldChange) > 1),]
ego_BP <- enrichGO(genelist_1$Symbol,OrgDb = 'org.Hs.eg.db',keyType = 'SYMBOL',ont = "BP",pvalueCutoff = 0.05)
ego_res_bp <- as.data.frame(ego_BP)
write.csv(ego_res_bp,'GOFC2.csv')

p1 <- barplot(ego_BP,showCategory = 20)
p1
p2 <- cnetplot(ego_BP, showCategory = 5)
p2
goplot(ego_BP)
ekg <- enrichKEGG(genelist_1$ENTREZID,organism = 'hsa')
ekg_res <- as.data.frame(ekg)
write.csv(ekg_res,'Kegg.csv')
p3 <- barplot(ekg,showCategory = 30)
p3
p4 <- dotplot(ekg,showCategory = 10)
p4

rownames(ekg_res) <- 1:nrow(ekg_res)
ekg_res$order=factor(rev(as.integer(rownames(ekg_res))),labels = rev(ekg_res$Description))
q <- ggplot(ekg_res,aes(y = order,x=count,fill=p.adjust)) +geom_bar(stat = 'identity',width = 0.7)+ scale_fill_gradient(low = 'red',high = 'blue') +labs(title='KEGG pathway Enrichment',x = 'Gene number', y = 'pathway')+
  theme(axis.title.x = element_text(face = "bold",size = 16),axis.title.y = element_text(face = "bold",size = 16),legend.title = element_text(face = "bold",size = 16)) +theme_bw()
                                                                              
                                                                             
  
