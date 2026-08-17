

Prefix=Songbird_RA
#Set the first position as 0 and give your thresholds from the second position, for example, 
#threshold=(0 3 3 3 3 3 3 3 3 3 3 3 0 2 3 3 3 3 3 3 3 3 2 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 3 2 2 3)

threshold=(0	2	2	50	50	10	3	8	4	50	50	50	50	50	50	50	50	50	2	2	50	50	10	10	10	10	50	50	50	50	50	50	50	50		100	2	10	50	50	50	50	50	50	50	50	50	50	50	50	50	50	50
)

#do a quick collpase around the UMI tags
for i in {1..51}; do
    
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
