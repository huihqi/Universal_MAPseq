

Prefix=Parrot_AAC
#Set the first position as 0 and give your thresholds from the second position, for example, 
#threshold=(0 3 3 3 3 3 3 3 3 3 3 3 0 2 3 3 3 3 3 3 3 3 2 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 3 2 2 3)

threshold=(0	2	5	2	2	5	20	20	80	80	80	80	3	4	2	3	20	10	20	80	80	20	80	80	20	1)

#do a quick collpase around the UMI tags
for i in {1..22}; do
    
    j=${threshold[$i]} 
    echo ${j}
    
    if [ "$j" != "1" ];then
	awk '$1 < '$j' {print NR}' ../${Prefix}_BC${i}_BU_counts.txt | head -1 >t
	thresh=$(cat t)
	echo $thresh
	head ../${Prefix}_BC${i}_BU_seq.txt -n $thresh | cut -b 1-30 | sort | uniq -c | sort -nr > ${Prefix}_BC${i}_quickout.txt
   else
	cut -b 1-30 ../${Prefix}_BC${i}_BU_seq.txt | sort | uniq -c | sort -nr > ${Prefix}_BC${i}_quickout.txt
   fi
done
