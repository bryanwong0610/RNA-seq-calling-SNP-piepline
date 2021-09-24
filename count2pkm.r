setwd('E:/JLZP_rnaseq_DE/')
BiocManager::install("GenomicFeatures")
#counts2TPM 
library("GenomicFeatures")
# withdrew gene lenth
txdb <- makeTxDbFromGFF("E:/secondpass/dbsnp/Homo_sapiens.GRCh38.99.gtf",format="auto")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
n=t(as.data.frame(exons_gene_lens))
write.table(n,file = 'gene_length.txt',col.names=F,row.names=T,quote=F,sep='\t')
gene_length <- read.table('gene_length.txt',header=F,sep='\t',check.names=F)
names(gene_length) <-c("gene_id","length")
# loading the counts file
JLZP_counts = read.csv(file = 'gene_counts_all.csv')
ctrl_counts = read.csv(file = 'ctrl_MII_GV.csv')
gene_counts_De_data = merge(JLZP_counts,ctrl_counts,by = "gene_id")
colnames(gene_counts_De_data) = c("gene_id","SYMBOL","Oocyte_1","Oocyte_2",
                                  "GV_1","GV_2","MI_1","MI_2","MII_1",
                                  "Ctrl_MII_1","Ctrl_MII_2","Ctrl_MII_3",
                                  "Ctrl_MII_4","Ctrl_MII_5","Ctrl_MII_6",
                                  "Ctrl_MII_7","Ctrl_MII_8","Ctrl_MII_9",
                                  "Ctrl_GV_1","Ctrl_GV_2","Ctrl_GV_3","Ctrl_GV_4"
                                  ,'Ctrl_GV_5',"Ctrl_GV_6","Ctrl_GV_7","Ctrl_GV_8"
                                  ,"Ctrl_GV_9","Ctrl_GV_10")
rownames(gene_counts_De_data) <- gene_counts_De_data[,1]
gene_counts_De_data_1 <- gene_counts_De_data[,-1]
gene_counts_De_data_1 <- gene_counts_De_data_1[,-1]
merge_lens_gene <- merge(gene_counts_De_data,gene_length,by='gene_id')
write.table(merge_lens_gene,'length_merge.txt')
dim(merge_lens_gene)
# calculate TPM&add symbol to it
a <- read.table('length_merge.txt')
rownames(a) <- a[,1]
kb <- a$length/1000
count <- a[,3:28]
rpk <- count/kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
write.table(tpm,'all_sample_tpm.txt')
b <- read.table('all_sample_tpm.txt',row.names = NULL)
colnames(b)[1] <- 'gene_id'
head(tpm)
symbols_id <- a[,1:2]
all_sample_tpm_addSYMBOL <- merge(symbols_id,b, by = 'gene_id')
write.csv(all_sample_tpm_addSYMBOL,'all_sample_tpm_addSYMBOL.csv',row.names = F,quote = F, sep ='\t')
write.table(all_sample_tpm_addSYMBOL,'all_sample_tpm_addSYMBOL.txt',row.names = F,quote = F, sep ='\t')
# calculate FPKM add symbol to it 
fpkm <- t(t(rpk)/colSums(count) * 10^6)
write.table(fpkm,'all_sample_fpkm.txt')
fkpm_before <- read.table('all_sample_fpkm.txt',row.names = NULL)
colnames(fkpm_before)[1] <- 'gene_id'
all_sample_fpkm_addSYMBOL <- merge(symbols_id,fkpm_before,by = 'gene_id')
write.table(all_sample_fpkm_addSYMBOL,'all_sample_fpkm_addSYMBOL.txt',row.names = F,quote = F,sep = '\t')
write.csv(all_sample_fpkm_addSYMBOL,'all_sample_fpkm_addSYMBOL.csv',row.names = F,quote = F)

#count2PCA 
res.pca <- PCA(gene_counts_De_data_1, scale.unit = TRUE, 
               graph = FALSE)


