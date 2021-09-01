Ref=~/reference
source activate rnaseq
cd /mnt/e/WAL-mutTrans/
#mkdir 1st_mapping
#mkdir 2nd_index
#mkdir 2nd_mapping
#echo "----- folder creation is complete ------"
#cat sample_list.txt | while read line
#do
#	gunzip /mnt/e/WAL-mutTrans/${line}/${line}_1.fq.gz
#	gunzip /mnt/e/WAL-mutTrans/${line}/${line}_2.fq.gz
#	STAR --runThreadN 8 --genomeDir $Ref --outSAMtype BAM Unsorted  \
#		--readFilesIn /mnt/e/WAL-mutTrans/${line}/${line}_1.fq /mnt/e/WAL-mutTrans/${line}/${line}_2.fq \
#		--outFileNamePrefix /mnt/e/WAL-mutTrans/1st_mapping/${line}.
#	echo " ------ ${line} mapping completed ------ "
#done
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
