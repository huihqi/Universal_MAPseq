#!/bin/bash

Prefix=Mouse_CP_Screening
for i in {1..65};do

	echo $i

	echo ${Prefix}_BC${i}_PE.txt
        
        echo >> ${Prefix}_BC${i}_PE.txt && cat ./Last_Batch/${Prefix}_BC${i}_PE.txt >> ${Prefix}_BC${i}_PE.txt

        awk "NR==${i}" CSI.txt > CSI
	
         C=$(cut -b 1-8 CSI)
         grep "$C"$ ${Prefix}_BC${i}_PE.txt | grep -v N | cut -b 1-42 | sort | uniq -c | sort -nr > ${Prefix}_BC${i}_BU.txt
         awk '{print $1}' ${Prefix}_BC${i}_BU.txt > ${Prefix}_BC${i}_BU_counts.txt
         awk '{print $2}' ${Prefix}_BC${i}_BU.txt > ${Prefix}_BC${i}_BU_seq.txt
done
