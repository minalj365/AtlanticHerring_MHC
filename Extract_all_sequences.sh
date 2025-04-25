#!/bin/bash
### extracting sequences
genome="/data1/CluHar/PacBio/CCS_data/denovo/hifiasm/RagTag/assemblies/"
gtf="/data1/CluHar/PacBio/CCS_data/MHC/analysis_for_paper/RagTag_annotations/Without_Null/"
out="/data1/CluHar/PacBio/CCS_data/MHC/analysis_for_paper/Raw_sequences/all/"


alpha="DAA1.1 DAA1.2 DAA1.3 DAA1.4 DAA1.5 DAA2.1 DAA2.2 DAA2.3 DAA4.1 DAA4.2 DAA4.3 DAA5.1 DAA7.1 DAA7.2 DAA8.1 DAA9.1 DBA3.1 DBA4.1 DBA4.2 DBA7.1"
beta="DAB1.1 DAB1.2 DAB1.3 DAB1.4 DAB1.5 DAB2.1 DAB2.2 DAB2.3 DAB4.1 DAB4.2 DAB4.3 DAB5.1a DAB5.1b DAB7.1 DAB7.2 DAB8.1 DAB9.1 DBB2.1 DBB3.1 DBB4.1 DBB4.2 DBB6.1 DBB7.1"
samples="Reference CS2_hap1 CS2_hap2 CS4_hap1 CS4_hap2 CS5_hap1 CS5_hap2 CS7_hap1 CS7_hap2 CS8_hap1 CS8_hap2 CS10_hap1 CS10_hap2 BS1_hap1 BS1_hap2 BS2_hap1 BS2_hap2 BS3_hap1 BS3_hap2 BS4_hap1 BS4_hap2 BS5_hap1 BS5_hap2 BS6_hap1 BS6_hap2 NSSH2_hap1 NSSH2_hap2 NSSH10_hap1 NSSH10_hap2"
exons="1 2 3 4"

##################################################################################### alpha ##############################################################################################

########## extract entire sequence
for i in ${samples}; do
    for j in ${alpha}; do
        grep "${j}" ${gtf}$i.gtf | ~/bin/gffread/gffread -w - --w-nocds -g ${genome}$i.fa | sed "s/>${j}T/>${i}_${j}/g" - | cat - >> ${out}alpha_entire_CDS.fa
    done;
done;


## Now change DA1.4 and 1.5 gene pairs from CS7_hap1 to DA1.1 and DA1.2 in CS7_hap2
sed -i -e 's/CS7_hap1_DAA1.4/CS7_hap2_DAA1.1/'  -e 's/CS7_hap1_DAA1.5/CS7_hap2_DAA1.2/' ${out}alpha_entire_CDS.fa

seqkit translate ${out}alpha_entire_CDS.fa > ${out}alpha_entire_AA.fa

########## extract exon 1,2,3,4 sequences
for i in ${samples}; do
    for j in ${alpha}; do
        for k in ${exons}; do
                grep "${j}" ${gtf}$i.gtf | grep "exon_number \"${k}\"" | ~/bin/gffread/gffread -w - --w-nocds -g ${genome}$i.fa | sed "s/>${j}T/>${i}_${j}/g" - | cat - >> ${out}alpha_E${k}_CDS.fa
        done;
    done;
done

for k in ${exons}; do seqkit translate ${out}alpha_E${k}_CDS.fa > ${out}alpha_E${k}_AA.fa; done


##################################################################################### beta ##############################################################################################
########## extract entire sequence
for i in ${samples}; do
    for j in ${beta}; do
        grep "${j}" ${gtf}$i.gtf | ~/bin/gffread/gffread -w - --w-nocds -g ${genome}$i.fa | sed "s/>${j}T/>${i}_${j}/g" - | cat - >> ${out}beta_entire_CDS.fa
    done;
done;

## Now change DA1.4 and 1.5 gene pairs from CS7_hap1 to DA1.1 and DA1.2 in CS7_hap2
sed -i -e 's/CS7_hap1_DAB1.4/CS7_hap2_DAB1.1/'  -e 's/CS7_hap1_DAB1.5/CS7_hap2_DAB1.2/' ${out}beta_entire_CDS.fa

seqkit translate ${out}beta_entire_CDS.fa > ${out}beta_entire_AA.fa

########## extract exon 1,2,3,4 sequences
for i in ${samples}; do
    for j in ${beta}; do
        for k in ${exons}; do
                grep "${j}" ${gtf}$i.gtf | grep "exon_number \"${k}\"" | ~/bin/gffread/gffread -w - --w-nocds -g ${genome}$i.fa | sed "s/>${j}T/>${i}_${j}/g" - | cat - >> ${out}beta_E${k}_CDS.fa
        done;
    done;
done

for k in ${exons}; do seqkit translate ${out}beta_E${k}_CDS.fa > ${out}beta_E${k}_AA.fa; done
