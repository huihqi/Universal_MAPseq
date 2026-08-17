

Prefix=Parrot_AACd
#Set the first position as 0 and give your thresholds from the second position, for example, 
#threshold=(0 3 3 3 3 3 3 3 3 3 3 3 0 2 3 3 3 3 3 3 3 3 2 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 3 2 2 3)

threshold=(0	3	20	2	10	10	10	80	80	80	50	8	10	80	80	80	80	80	80	80	80	80	80	8	8	80	100	80	80	80	80	80	10	100	80	100	10	20	8	100	3	100	3	10	10	80	80	80	10	20	80)

#do a quick collpase around the UMI tags
for i in {1..50}; do
    
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
