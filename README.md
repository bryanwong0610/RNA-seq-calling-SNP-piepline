# RNA-seq-calling-SNP-piepline
With star 2-pass mapping mode calling SNP by using the transcriptome sequencing file
# 1st mapping index generation
For 1st mapping ，the index was built by STAR 2.7.9a. 
The mapping command line for the index which build for 1-st mapping is shown as below:
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir /home/bryan0610/reference 
--genomeFastaFiles /home/bryan0610/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa 
--sjdbGTFfile /home/bryan0610/reference/Homo_sapiens.GRCh38.99.gtf --sjdbOverhang 99
# 2nd mapping protocol
see as mapping.sh
