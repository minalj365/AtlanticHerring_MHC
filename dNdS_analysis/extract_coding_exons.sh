#!/bin/sh

genome="/data1/CluHar/PacBio/CCS_data/denovo/hifiasm/RagTag/assemblies/"
gtf="/data1/CluHar/PacBio/CCS_data/MHC/analysis_for_paper/RagTag_annotations/Without_Null/"
out="/data1/CluHar/PacBio/CCS_data/MHC/analysis_for_paper/Raw_sequences/complete_exons/"

samples="Reference CS2_hap1 CS2_hap2 CS4_hap1 CS4_hap2 CS5_hap1 CS5_hap2 CS7_hap1 CS7_hap2 CS8_hap1 CS8_hap2 CS10_hap1 CS10_hap2 BS1_hap1 BS1_hap2 BS2_hap1 BS2_hap2 BS3_hap1 BS3_hap2 BS4_hap1 BS4_hap2 BS5_hap1 BS5_hap2 BS6_hap1 BS6_hap2 NSSH2_hap1 NSSH2_hap2 NSSH10_hap1 NSSH10_hap2"
genes="DAA1.1 DAA1.2 DAA1.3 DAA1.4 DAA1.5 DAA2.1 DAA2.2 DAA2.3 DAA4.1 DAA4.2 DAA5.1 DAA7.1 DAA8.1 DAA8.2 DAA9.1 DBA3.1 DBA4.1 DBA4.2 DBA6.1 DAB1.1 DAB1.2 DAB1.3 DAB1.4 DAB1.5 DAB2.1 DAB2.2 DAB2.3 DAB4.1 DAB4.2 DAB5.1a DAB5.1b DAB7.1 DAB8.1 DAB8.2 DAB9.1 DBB2.1 DBB3.1 DBB4.1 DBB4.2 DBB6.1"
exons="1 2 3 4"

mkdir -p ${out}concatenate # to store combined sequences

# get all exons with gffread
for i in ${samples}; do
    for j in ${genes}; do
        for k in ${exons}; do
                grep "${j}" ${gtf}$i.gtf | grep "exon_number \"${k}\"" | ~/bin/gffread/gffread -w - --w-nocds -g ${genome}$i.fa | sed "s/>/>${i}_/g" - > ${out}${j}_${i}_E${k}.fa
        done;
    done;
done

wait;

## Now change DA1.4 and 1.5 gene pairs from CS7_hap1 to DA1.1 and DA1.2 in CS7_hap2

from_genes=("DAA1.4" "DAA1.5" "DAB1.4" "DAB1.5")
to_genes=("DAA1.1" "DAA1.2" "DAB1.1" "DAB1.2")

for exon in $exons; do
  for i in "${!from_genes[@]}"; do
    from_gene=${from_genes[i]}
    to_gene=${to_genes[i]}
    old_file="${from_gene}_CS7_hap1_E${exon}.fa"
    new_file="${to_gene}_CS7_hap2_E${exon}.fa"
    
    # Renaming the file
    mv "$old_file" "$new_file"
    
    # Using sed to replace content in the new file
    sed -i -e "s/CS7_hap1_${from_gene}/CS7_hap2_${to_gene}/g" "$new_file"
  done
done


# index all exon fasta files
for i in ${out}*fa; do samtools faidx ${i}; done

########## make exon 2 files ##########

# extract last nt from exon1 AND extract all but last nt from exon2
for i in ${samples}; do
    for j in ${genes}; do
                echo "samtools faidx ${out}${j}_${i}_E1.fa ${i}_${j}T:`cut -f2 ${out}${j}_${i}_E1.fa.fai`-`cut -f2 ${out}${j}_${i}_E1.fa.fai` | sed "s/T:`cut -f2 ${out}${j}_${i}_E1.fa.fai`-`cut -f2 ${out}${j}_${i}_E1.fa.fai`//" > ${out}${j}_${i}_E1_for_E2.fa" >> to_bash1.sh
                echo "samtools faidx ${out}${j}_${i}_E2.fa ${i}_${j}T:1-$(expr `cut -f2 ${out}${j}_${i}_E2.fa.fai` - 1) | sed "s/T:1-$(expr `cut -f2 ${out}${j}_${i}_E2.fa.fai` - 1)//" > ${out}${j}_${i}_E2_for_E2.fa" >> to_bash2.sh
        done;
done;

########## make exon 3 files ##########

# extract last nt from exon2 AND extract all but last nt from exon3
for i in ${samples}; do
    for j in ${genes}; do
                echo "samtools faidx ${out}${j}_${i}_E2.fa ${i}_${j}T:`cut -f2 ${out}${j}_${i}_E2.fa.fai`-`cut -f2 ${out}${j}_${i}_E2.fa.fai` | sed "s/T:`cut -f2 ${out}${j}_${i}_E2.fa.fai`-`cut -f2 ${out}${j}_${i}_E2.fa.fai`//" > ${out}${j}_${i}_E2_for_E3.fa" >> to_bash3.sh
                echo "samtools faidx ${out}${j}_${i}_E3.fa ${i}_${j}T:1-$(expr `cut -f2 ${out}${j}_${i}_E3.fa.fai` - 1) | sed "s/T:1-$(expr `cut -f2 ${out}${j}_${i}_E3.fa.fai` - 1)//" > ${out}${j}_${i}_E3_for_E3.fa" >> to_bash4.sh
        done;
done;


bash to_bash1.sh
bash to_bash2.sh
bash to_bash3.sh
bash to_bash4.sh

wait;


for i in ${samples}; do
    for j in ${genes}; do
                # combine exon1 and exon2
                seqkit concat ${out}${j}_${i}_E1_for_E2.fa ${out}${j}_${i}_E2_for_E2.fa > ${out}concatenate/${j}_${i}_E2.fa
                # combine exon2 and exon3
                seqkit concat ${out}${j}_${i}_E2_for_E3.fa ${out}${j}_${i}_E3_for_E3.fa > ${out}concatenate/${j}_${i}_E3.fa
        done;
done;

############## clean the data ###############

mkdir ${out}sequences
for j in ${genes}; do
        cat ${out}concatenate/${j}*E2.fa > ${out}sequences/${j}_E2_CDS.fa
        cat ${out}concatenate/${j}*E3.fa > ${out}sequences/${j}_E3_CDS.fa
        seqkit translate ${out}sequences/${j}_E2_CDS.fa > ${out}sequences/${j}_E2_AA.fa
        seqkit translate ${out}sequences/${j}_E3_CDS.fa > ${out}sequences/${j}_E3_AA.fa
done;

wait;

rm -r ${out}concatenate ${out}*fa*  ${out}to_bash*
mv ${out}sequences/*fa ${out}
rm -r ${out}sequences
rm *1.4* *1.5*
