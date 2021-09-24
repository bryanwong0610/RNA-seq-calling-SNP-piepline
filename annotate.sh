dbsnp=/mnt/e/secondpass/annovar/annovar_annotated_db/
base=/mnt/e/WAL-mutTrans

source activate rnaseq
echo "----- preparing annotation process -------"
cat $base/sample_list.txt | while read line
do
	table_annovar.pl $base/${line}.gatk_call/${line}_calling_SNP_INDEL.vcf $dbsnp  \
		-buildver hg38 -out $base/${line}.gatk_call/${line}.annovar. \
		-remove -protocol refGene,cytoBand,dbnsfp41a,exac03,clinvar_20210123,avsnp150,esp6500siv2_all,gnomad211_exome,gnomad30_genome,ALL.sites.2015_08,EAS.sites.2015_08,AFR.sites.2015_08,AMR.sites.2015_08,EUR.sites.2015_08,SAS.sites.2015_08,dbscsnv11 \
		-operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish -thread 6 -dot2underline
	echo "------ ${line} annotation complete ------"
done
echo "------- all process have been done ------"
