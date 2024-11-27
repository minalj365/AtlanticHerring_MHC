#!/bin/sh

genome="/data1/CluHar/PacBio/CCS_data/denovo/hifiasm/RagTag/assemblies/"
gtf="/data1/CluHar/PacBio/CCS_data/MHC/analysis_for_paper/RagTag_annotations/Without_Null/"
out="/data1/CluHar/PacBio/CCS_data/MHC/analysis_for_paper/Raw_sequences/pi/"

samples="Reference CS2_hap1 CS2_hap2 CS4_hap1 CS4_hap2 CS5_hap1 CS5_hap2 CS7_hap1 CS7_hap2 CS8_hap1 CS8_hap2 CS10_hap1 CS10_hap2 BS1_hap1 BS1_hap2 BS2_hap1 BS2_hap2 BS3_hap1 BS3_hap2 BS4_hap1 BS4_hap2 BS5_hap1 BS5_hap2 BS6_hap1 BS6_hap2 NSSH2_hap1 NSSH2_hap2 NSSH10_hap1 NSSH10_hap2"
genes="DAA8.1 DAA1.1 DAA1.2 DAA1.3 DAA1.4 DAA1.5 DAA2.1 DAA2.2 DAA2.3 DAA4.1 DAA4.2 DAA5.1 DAA7.1 DAA8.1 DAA8.2 DAA9.1 DBA3.1 DBA4.1 DBA4.2 DBA6.1 DAB8.1 DAB1.1 DAB1.2 DAB1.3 DAB1.4 DAB1.5 DAB2.1 DAB2.2 DAB2.3 DAB4.1 DAB4.2 DAB5.1a DAB5.1b DAB7.1 DAB8.1 DAB8.2 DAB9.1 DBB2.1 DBB3.1 DBB4.1 DBB4.2 DBB6.1"
exons="1 2 3 4"

########## extract entire sequence
for i in ${samples}; do
    for j in ${genes}; do
        grep "${j}" ${gtf}$i.gtf | ~/bin/gffread/gffread -w - --w-nocds -g ${genome}$i.fa | sed "s/>${j}T/>${i}_${j}/g" - | cat - >> ${out}${j}_entire_CDS.fa
    done;
done;

## Now change DA1.4 and 1.5 gene pairs from CS7_hap1 to DA1.1 and DA1.2 in CS7_hap2 AND add them to DA1.1 and DA1.2 files
# alpha
sed 's/CS7_hap1_DAA1.4/CS7_hap2_DAA1.1/' ${out}DAA1.4_entire_CDS.fa >> ${out}DAA1.1_entire_CDS.fa
sed 's/CS7_hap1_DAA1.5/CS7_hap2_DAA1.2/' ${out}DAA1.5_entire_CDS.fa >> ${out}DAA1.2_entire_CDS.fa

# beta
sed 's/CS7_hap1_DAB1.4/CS7_hap2_DAB1.1/' ${out}DAB1.4_entire_CDS.fa >> ${out}DAB1.1_entire_CDS.fa
sed 's/CS7_hap1_DAB1.5/CS7_hap2_DAB1.2/' ${out}DAB1.5_entire_CDS.fa >> ${out}DAB1.2_entire_CDS.fa

seqkit translate ${out}${j}_entire_CDS.fa > ${out}${j}_entire_AA.fa

########## extract exon 1,2,3,4 sequences
for i in ${samples}; do
    for j in ${genes}; do
        for k in ${exons}; do
                grep "${j}" ${gtf}$i.gtf | grep "exon_number \"${k}\"" | ~/bin/gffread/gffread -w - --w-nocds -g ${genome}$i.fa | sed "s/>${j}T/>${i}_${j}/g" - | cat - >> ${out}${j}_E${k}_CDS.fa
        done;
    done;
done

## Now change DA1.4 and 1.5 gene pairs from CS7_hap1 to DA1.1 and DA1.2 in CS7_hap2 AND add them to DA1.1 and DA1.2 files
# alpha
sed 's/CS7_hap1_DAA1.4/CS7_hap2_DAA1.1/' ${out}DAA1.4_E1_CDS.fa >> ${out}DAA1.1_E1_CDS.fa
sed 's/CS7_hap1_DAA1.4/CS7_hap2_DAA1.1/' ${out}DAA1.4_E2_CDS.fa >> ${out}DAA1.1_E2_CDS.fa
sed 's/CS7_hap1_DAA1.4/CS7_hap2_DAA1.1/' ${out}DAA1.4_E3_CDS.fa >> ${out}DAA1.1_E3_CDS.fa
sed 's/CS7_hap1_DAA1.4/CS7_hap2_DAA1.1/' ${out}DAA1.4_E4_CDS.fa >> ${out}DAA1.1_E4_CDS.fa

sed 's/CS7_hap1_DAA1.5/CS7_hap2_DAA1.2/' ${out}DAA1.5_E1_CDS.fa >> ${out}DAA1.2_E1_CDS.fa
sed 's/CS7_hap1_DAA1.5/CS7_hap2_DAA1.2/' ${out}DAA1.5_E2_CDS.fa >> ${out}DAA1.2_E2_CDS.fa
sed 's/CS7_hap1_DAA1.5/CS7_hap2_DAA1.2/' ${out}DAA1.5_E3_CDS.fa >> ${out}DAA1.2_E3_CDS.fa
sed 's/CS7_hap1_DAA1.5/CS7_hap2_DAA1.2/' ${out}DAA1.5_E4_CDS.fa >> ${out}DAA1.2_E4_CDS.fa

# beta
sed 's/CS7_hap1_DAB1.4/CS7_hap2_DAB1.1/' ${out}DAB1.4_E1_CDS.fa >> ${out}DAB1.1_E1_CDS.fa
sed 's/CS7_hap1_DAB1.4/CS7_hap2_DAB1.1/' ${out}DAB1.4_E2_CDS.fa >> ${out}DAB1.1_E2_CDS.fa
sed 's/CS7_hap1_DAB1.4/CS7_hap2_DAB1.1/' ${out}DAB1.4_E3_CDS.fa >> ${out}DAB1.1_E3_CDS.fa
sed 's/CS7_hap1_DAB1.4/CS7_hap2_DAB1.1/' ${out}DAB1.4_E4_CDS.fa >> ${out}DAB1.1_E4_CDS.fa

sed 's/CS7_hap1_DAB1.5/CS7_hap2_DAB1.2/' ${out}DAB1.5_E1_CDS.fa >> ${out}DAB1.2_E1_CDS.fa
sed 's/CS7_hap1_DAB1.5/CS7_hap2_DAB1.2/' ${out}DAB1.5_E2_CDS.fa >> ${out}DAB1.2_E2_CDS.fa
sed 's/CS7_hap1_DAB1.5/CS7_hap2_DAB1.2/' ${out}DAB1.5_E3_CDS.fa >> ${out}DAB1.2_E3_CDS.fa
sed 's/CS7_hap1_DAB1.5/CS7_hap2_DAB1.2/' ${out}DAB1.5_E4_CDS.fa >> ${out}DAB1.2_E4_CDS.fa


for k in ${exons}; do seqkit translate ${out}${j}_E${k}_CDS.fa > ${out}${j}_E${k}_AA.fa; done


# remove DA1.4 and DA1.5 files
rm *1.4* *1.5*
