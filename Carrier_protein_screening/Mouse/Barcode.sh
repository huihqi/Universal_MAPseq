# filter sequence by carrier protein
Prefix=Mouse_CP_Screening
brain=sample9
for i in 33 34 35 36 61;do
echo ${i}
        for j in {1..31}; do
        awk "NR==${j}" CP_index.txt > CP_index	
        C=$(cut -b 1-6 CP_index)
        echo ${C}
        Carrier_Protein_pool=(0 AP2B1	AP2M1	CACNG4	CADM1	DNM1	ERC2	GAP43	PTPRA	RAB3A	RAB6B	SNAP25	STX1B	STX6	SYNAPSINI	VAMP2_SYNAPTOPHYSIN	SYNAPTOPHYSIN	SYNGR3	SYT1	VAMP1	VAMP2	VAMP4	STX1B_4nl	STX1B_PP7	SYNAPTOPHYSIN_4nl	SYNAPTOPHYSIN_PP7	VAMP2_4nl	VAMP2_PP7	DNM1_4nl	DNM1_PP7	mVAMP2_HK16	MAPP_JK100)
        Carrier_Protein=${Carrier_Protein_pool[$j]} 
        echo ${Carrier_Protein}
        in=${Prefix}_BC${i}_quickout.txt
        grep "$C"$ ${in} > ${brain}_${Prefix}_BC${i}_${Carrier_Protein}_quickprocessed.txt
        awk '{print $1}' ${brain}_${Prefix}_BC${i}_${Carrier_Protein}_quickprocessed.txt > ${brain}_${Prefix}_${i}_${Carrier_Protein}_counts.txt
        awk '{print $2}' ${brain}_${Prefix}_BC${i}_${Carrier_Protein}_quickprocessed.txt > ${brain}_${Prefix}_${i}_${Carrier_Protein}_seq.txt
        nl ${brain}_${Prefix}_${i}_${Carrier_Protein}_seq.txt | awk '{print ">" $1 "\n" $2}' > ${brain}_${Prefix}_${i}${Carrier_Protein}_fasta2u.txt
        bowtie-build -q ${brain}_${Prefix}_${i}${Carrier_Protein}_fasta2u.txt indexes/${brain}_${i}${Carrier_Protein}_fasta2u
        bowtie -v 3 -p 10 -f --best -a indexes/${brain}_${i}${Carrier_Protein}_fasta2u ${brain}_${Prefix}_${i}${Carrier_Protein}_fasta2u.txt ${brain}_bowtiealignment${i}_${Carrier_Protein}_2u.txt
        awk '{print $1}' ${brain}_bowtiealignment${i}_${Carrier_Protein}_2u.txt > ${brain}_bowtie${i}_${Carrier_Protein}_2u_1.txt
        awk '{print $3}' ${brain}_bowtiealignment${i}_${Carrier_Protein}_2u.txt > ${brain}_bowtie${i}_${Carrier_Protein}_2u_3.txt
        done
done



 
