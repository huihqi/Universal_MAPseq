Prefix='Caleb_Bird'
#Set the first position as 0 and give your thresholds from the second position, for example, 
threshold=(0 2 3 5 3 5 5 4 0 10 6 20 20 20 10 20 20 10 10 10 10 9 10 10 10 9 9 5 5 4 4 4 4 10 10 10 6 1 0 0 0)


#do a quick collpase around the UMI tags
for i in {1..40}; do
    j=${threshold[$i]} 
    echo ${j}
   # k=$[$i-1]
    if [ "$j" != "1" ];then
	awk '$1 < '$j' {print NR}' ../${Prefix}_BC${i}_BU_counts.txt | head -1 >t
	thresh=$(cat t)
	echo $thresh
	head ../${Prefix}_BC${i}_BU_seq.txt -n $thresh | cut -b 1-28 | sort | uniq -c | sort -nr > ${Prefix}_BC${i}_quickout.txt
   else
	cut -b 1-28 ../${Prefix}_BC${i}_BU_seq.txt | sort | uniq -c | sort -nr > ${Prefix}_BC${i}_quickout.txt
   fi
done
