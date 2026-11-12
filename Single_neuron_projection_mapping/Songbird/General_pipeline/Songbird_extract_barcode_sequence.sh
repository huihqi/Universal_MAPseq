#!/bin/bash

Prefix=Songbird_RA
for i in {1..57};do

	echo $i

	echo *[^0-9]${i}_S[0-9]*_L007_R1_001.fastq.gz
        echo *[^0-9]${i}_S[0-9]*_L007_R2_001.fastq.gz
	zcat *[^0-9]${i}_S[0-9]*_L007_R1_001.fastq.gz | awk "NR%4==2" | cut -b 1-28 > ${Prefix}_${i}_R1_stripped.txt
 	zcat *[^0-9]${i}_S[0-9]*_L007_R2_001.fastq.gz | awk "NR%4==2" | cut -b 1-20 | paste -d '' ${Prefix}_${i}_R1_stripped.txt - > ${Prefix}_BC${i}_PE.txt
        awk "NR==${i}" CSI_songbird.txt > CSI
	
         C=$(cut -b 1-8 CSI)
         grep "$C"$ ${Prefix}_BC${i}_PE.txt | grep -v N | cut -b 1-40 | sort | uniq -c | sort -nr > ${Prefix}_BC${i}_BU.txt
         awk '{print $1}' ${Prefix}_BC${i}_BU.txt > ${Prefix}_BC${i}_BU_counts.txt
         awk '{print $2}' ${Prefix}_BC${i}_BU.txt > ${Prefix}_BC${i}_BU_seq.txt
done
