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
