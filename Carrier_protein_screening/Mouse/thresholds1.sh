
cd thresholds

Prefix=Mouse_CP_Screening
#Set the first position as 0 and give your thresholds from the second position, for example, 
#threshold=(0 3 3 3 3 3 3 3 3 3 3 3 0 2 3 3 3 3 3 3 3 3 2 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 3 2 2 3)

threshold=(0	50	50	50	50	50	50	50	50	50	50	50	50	50	50	50	50	50	50	50	50	10	10	10	1	10	50	50	50	20	50	50	2	3	8	20	50	10	50	80	80	7	50	10	50	7	50	50	50	8	50	20	50	2	2	2	2	2	2	2	2	2	2	2	2	2)
#do a quick collpase around the UMI tags
for i in {1..65}; do
    j=${threshold[$i]} 
    echo ${j}
    
    if [ "$j" != "1" ];then
	awk '$1 < '$j' {print NR}' ../${Prefix}_BC${i}_BU_counts.txt | head -1 >t
	thresh=$(cat t)
	echo $thresh
	head ../${Prefix}_BC${i}_BU_seq.txt -n $thresh | cut -b 1-26 | sort | uniq -c | sort -nr > ${Prefix}_BC${i}_quickout.txt
   else
	cut -b 1-26 ../${Prefix}_BC${i}_BU_seq.txt | sort | uniq -c | sort -nr > ${Prefix}_BC${i}_quickout.txt
   fi
done
