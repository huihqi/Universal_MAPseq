Prefix=Songbird_RA
for i in {1..51};do
echo ${i}
in=${Prefix}_BC${i}_quickout.txt
grep -v 'GGTTCT[A-Z]\{4\}$' ${in} > ${Prefix}_BC${i}_quickprocessed.txt
awk '{print $1}' ${Prefix}_BC${i}_quickprocessed.txt > ${Prefix}_${i}_counts.txt
awk '{print $2}' ${Prefix}_BC${i}_quickprocessed.txt > ${Prefix}_${i}_seq.txt
nl ${Prefix}_${i}_seq.txt | awk '{print ">" $1 "\n" $2}' > ${Prefix}_${i}fasta2u.txt
bowtie-build -q ${Prefix}_${i}fasta2u.txt indexes/${i}fasta2u
bowtie -v 3 -p 10 -f --best -a indexes/${i}fasta2u ${Prefix}_${i}fasta2u.txt bowtiealignment${i}_2u.txt
awk '{print $1}' bowtiealignment${i}_2u.txt > bowtie${i}_2u_1.txt
awk '{print $3}' bowtiealignment${i}_2u.txt > bowtie${i}_2u_3.txt
done
