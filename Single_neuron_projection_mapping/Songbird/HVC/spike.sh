Prefix=Caleb_Bird
for i in {1..40};do
echo ${i}
in=${Prefix}_BC${i}_quickout.txt
grep 'CGTC$' ${in} > ${Prefix}_spikesBC${i}_quickprocessed.txt
awk '{print $1}' ${Prefix}_spikesBC${i}_quickprocessed.txt > ${Prefix}_spikes${i}_counts.txt
awk '{print $2}' ${Prefix}_spikesBC${i}_quickprocessed.txt > ${Prefix}_spikes${i}_seq.txt
nl ${Prefix}_spikes${i}_seq.txt | awk '{print ">" $1 "\n" $2}' > ${Prefix}_spikes${i}fasta2u.txt
bowtie-build -q ${Prefix}_spikes${i}fasta2u.txt indexes/spikes${i}fasta2u
bowtie -v 3 -p 10 -f --best -a indexes/spikes${i}fasta2u ${Prefix}_spikes${i}fasta2u.txt bowtiealignmentspikes${i}_2u.txt
awk '{print $1}' bowtiealignmentspikes${i}_2u.txt > bowtiespikes${i}_2u_1.txt
awk '{print $3}' bowtiealignmentspikes${i}_2u.txt > bowtiespikes${i}_2u_3.txt
done
