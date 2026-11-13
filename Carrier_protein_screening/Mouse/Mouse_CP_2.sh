#!/bin/bash
# Filter sequences based on CSI
Prefix=Mouse_CP_Screening
for i in {1..65};do



        awk "NR==${i}" CSI.txt > CSI
	
         C=$(cut -b 1-8 CSI)
         grep "$C"$ ${Prefix}_BC${i}_PE.txt |  cut -b 1-26,31-42 | grep -v N | sort | uniq -c | sort -nr > ${Prefix}_BC${i}_BU.txt
         awk '{print $1}' ${Prefix}_BC${i}_BU.txt > ${Prefix}_BC${i}_BU_counts.txt
         awk '{print $2}' ${Prefix}_BC${i}_BU.txt > ${Prefix}_BC${i}_BU_seq.txt
done
