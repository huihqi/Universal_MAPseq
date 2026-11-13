#!/bin/bash
# Load Zebra finch carrier protein screening fastq files
# After this step, you will get barcodes and UMI sequence file
Prefix=Bird_CP_Screening
for i in {1..48};do

	echo $i

	echo HH48-${i}_S[0-9]*_L007_R1_001.fastq.gz
        echo HH48-${i}_S[0-9]*_L007_R2_001.fastq.gz
	zcat HH48-${i}_S[0-9]*_L007_R1_001.fastq.gz | awk "NR%4==2" | cut -b 1-30 > ${Prefix}_${i}_R1_stripped.txt
 	zcat HH48-${i}_S[0-9]*_L007_R2_001.fastq.gz | awk "NR%4==2" | cut -b 1-20 | paste -d '' ${Prefix}_${i}_R1_stripped.txt - > ${Prefix}_BC${i}_PE.txt
        awk "NR==${i}" CSI.txt > CSI
	
         C=$(cut -b 1-8 CSI)
         grep "$C"$ ${Prefix}_BC${i}_PE.txt |  cut -b 1-26,31-42 | grep -v N | sort | uniq -c | sort -nr > ${Prefix}_BC${i}_BU.txt
         awk '{print $1}' ${Prefix}_BC${i}_BU.txt > ${Prefix}_BC${i}_BU_counts.txt
         awk '{print $2}' ${Prefix}_BC${i}_BU.txt > ${Prefix}_BC${i}_BU_seq.txt
done
