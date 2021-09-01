# RNA-seq-calling-SNP-piepline
With star 2-pass mapping mode calling SNP by using the transcriptome sequencing file
# 1st mapping index generation
For 1st mapping ，the index was built by STAR 2.7.9a. 
The mapping command line for the index which build for 1-st mapping is shown as below:
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


