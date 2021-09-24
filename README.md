# RNA-seq-calling-SNP-piepline

With star 2-pass mapping mode calling SNP by using the transcriptome sequencing file

# 1st mapping index generation

For 1st mapping ，the index was built by STAR 2.7.9a. 

The mapping command line for the index which built for 1-st mapping is shown as below:
       
       STAR --runThreadN 6 --runMode genomeGenerate \ --genomeDir /home/bryan0610/reference 
        --genomeFastaFiles \ /home/bryan0610/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa 
        --sjdbGTFfile \ /home/bryan0610/reference/Homo_sapiens.GRCh38.99.gtf --sjdbOverhang 99
# 2-pass mapping protocol
Due to the limitation of WSL ( Windows subsysterm for Linux), 
the "*.fq.gz" file could not be a suitable file for star mapping process. 

So pre-decompress the fq file by using the function gunzip for the subsequent analysis

Refer to the key commands in the demo file mapping.sh


    Ref=~/reference

    source activate rnaseq

    cd /mnt/e/WAL-mutTrans/

    mkdir 1st_mapping

    mkdir 2nd_index

    mkdir 2nd_mapping

    echo "----- folder creation is complete ------"

    cat sample_list.txt | while read line
      do
        gunzip /mnt/e/WAL-mutTrans/${line}/${line}_1.fq.gz
   
        gunzip /mnt/e/WAL-mutTrans/${line}/${line}_2.fq.gz
   
        STAR --runThreadN 8 --genomeDir $Ref --outSAMtype BAM Unsorted  \
                --readFilesIn /mnt/e/WAL-mutTrans/${line}/${line}_1.fq /mnt/e/WAL-mutTrans/${line}/${line}_2.fq \
                --outFileNamePrefix /mnt/e/WAL-mutTrans/1st_mapping/${line}.
 
        echo " ------ ${line} mapping completed ------ "
      done

     echo " ------ Preparing 2nd-pass index generate ------"

     cat sample_list.txt | while read line
       do
         
         mkdir ${line}_2nd_index
         
         STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /mnt/e/WAL-mutTrans/2nd_index/${line}_2nd_index \
                --genomeFastaFiles $Ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
              --sjdbGTFfile $Ref/Homo_sapiens.GRCh38.99.gtf \
             --sjdbFileChrStartEnd /mnt/e/WAL-mutTrans/1st_mapping/${line}.SJ.out.tab
       
     echo "------ ${line}_2nd_pass index generated ------"
      
     done
    
     echo "------ Starting 2nd-pass mapping ------"
    
     cat sample_list.txt | while read line
        
        do
            
            STAR --runThreadN 8 --genomeDir /mnt/e/WAL-mutTrans/2nd_index/${line}_2nd_index/ --outSAMtype BAM Unsorted \
                --readFilesIn /mnt/e/WAL-mutTrans/${line}/${line}_1.fq /mnt/e/WAL-mutTrans/${line}/${line}_2.fq \
                --outFileNamePrefix /mnt/e/WAL-mutTrans/2nd_mapping/${line}.
        
            echo " ------ ${line} 2nd-pass mapping completed ------ "

        done

      echo "  ----------   Hoooooray --------- " 

# 3 GATK mapping pipline

The demo of calling SNP&INDEL script by using GATK 

       # ref also need index so use samtools to generate the ref.genome index
       # parameter "samtools faidx ***.fa" to generate
       # and the dictionary of ref genome is also needded use gatk CreateSequenceDictionary to generate
       # the parameter " gatk CreateSequenceDictionary -I ****.fa -O ***.dict " (which same prefix name but not has the .fa)
       Bamdir=/mnt/e/WAL-mutTrans/2nd_mapping
       samplelist=/mnt/e/WAL-mutTrans/sample_list.txt
       Ref=/mnt/e/secondpass/dbsnp
       base=/mnt/e/WAL-mutTrans

       echo " ------ folder preparing ------- "
       cat $samplelist | while read line 
       do
	       mkdir $base/${line}.gatk_call
	       source activate rnaseq
	       echo "------ ${line} gatk starting ------"
	       gatk AddOrReplaceReadGroups -I $Bamdir/${line}.Aligned.out.bam \
		       -SO coordinate -RGID ${line} -RGLB ${line}_rna_paired  \
		       -RGPL BGI_smartseq -RGPU Smartseq2 -RGSM PolyPN_oocyte \
		       -O $base/${line}.gatk_call/${line}_added_sorted.bam 
	       echo "------ ${line}'s bam Sorted ------"
	       echo "------ preparing MarkDuplicates ------"
	       gatk MarkDuplicates -I $base/${line}.gatk_call/${line}_added_sorted.bam \
		       -O $base/${line}.gatk_call/${line}_dedupped_added_sorted.bam \
		       -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT \
		       -M ${line}_output.metrics
	       echo "------ ${line}_bam has been marked the Duplicates ------"
	       echo "------ preparing split the Intron AND Exon ------"
	       gatk SplitNCigarReads -R $Ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
		       -I $base/${line}.gatk_call/${line}_dedupped_added_sorted.bam \
		       -O $base/${line}.gatk_call/${line}_splicted_sorted.bam 
	       echo "------ ${line} split completed ------"
	       echo "------ preparing BQSR ------"
	       gatk BaseRecalibrator -R $Ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
		       -I $base/${line}.gatk_call/${line}_splicted_sorted.bam \
		       -known-sites $Ref/no_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
		       -known-sites $Ref/no_1000G_omni2.5.hg38.vcf.gz \
		       -known-sites $Ref/no_dbsnp_146.hg38.vcf.gz \
		       -known-sites $Ref/no_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
		       -O $base/${line}.gatk_call/${line}_recall_data.table
	       echo "------ ${line} BQSR complete ------"
	       echo "------ preparing apply BQSR result ------"
	       gatk ApplyBQSR -R $Ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
		       -I $base/${line}.gatk_call/${line}_splicted_sorted.bam \
		       -bqsr-recal-file $base/${line}.gatk_call/${line}_recall_data.table \
		       -O $base/${line}.gatk_call/${line}_recal_splicted_sorted.bam
	       echo "------ ${line}_bam has been BQSRed ------"
	       echo "------ preparing calling the SNP & INDEL ------"
	       gatk HaplotypeCaller -R $Ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
		       -I $base/${line}.gatk_call/${line}_recal_splicted_sorted.bam \
		       -dont-use-soft-clipped-bases \
		       -stand-call-conf 20 \
		       -O  $base/${line}.gatk_call/${line}_calling_SNP_INDEL.vcf
              done
              echo "------ SNP&INDEL calling complete ------"

# 4 Annotated SNP 
The annotation of VCF file which calling by GATK in RNA-seq situation.
      
       dbsnp=/mnt/e/secondpass/annovar/annovar_annotated_db/
       base=/mnt/e/WAL-mutTrans

       source activate rnaseq
       echo "----- preparing annotation process -------"
       cat $base/sample_list.txt | while read line
       do
	       table_annovar.pl $base/${line}.gatk_call/${line}_calling_SNP_INDEL.vcf $dbsnp  \
		       -buildver hg38 -out $base/${line}.gatk_call/${line}.annovar. \
		       -remove -protocol                  refGene,cytoBand,dbnsfp41a,exac03,clinvar_20210123,avsnp150,esp6500siv2_all,gnomad211_exome,
                     gnomad30_genome,ALL.sites.2015_08,EAS.sites.2015_08,AFR.sites.2015_08,AMR.sites.2015_08,EUR.sites.2015_08,SAS.sites.2015_08,dbscsnv11 \
		       -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish -thread 6 -dot2underline
	       echo "------ ${line} annotation complete ------"
       done
       echo "------- all process have been done ------"
 # 5 RNA-seq BAM to counts was using STAR-Stringtie get the Counts matrix
 
 # 6 Deferent Exprission Gene Analysis with R 
       #1 Counts to TPM
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
              colnames(gene_counts_De_data) = c("gene_id","SYMBOL",gruop
                                   )
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
       # PCA 
       
       
 
