#!/bin/bash
# Load Mouse Carrier protein screening Fastq files and combine two batches of sequencing files
Prefix=Mouse_CP_Screening
for i in {1..65};do

	echo $i

	echo HH65-${i}_S[0-9]*_L008_R1_001.fastq.gz
        echo HH65-${i}_S[0-9]*_L008_R2_001.fastq.gz
	zcat HH65-${i}_S[0-9]*_L008_R1_001.fastq.gz | awk "NR%4==2" | cut -b 1-30 > ${Prefix}_${i}_R1_stripped.txt
 	zcat HH65-${i}_S[0-9]*_L008_R2_001.fastq.gz | awk "NR%4==2" | cut -b 1-20 | paste -d '' ${Prefix}_${i}_R1_stripped.txt - > ${Prefix}_BC${i}_PE.txt

done
