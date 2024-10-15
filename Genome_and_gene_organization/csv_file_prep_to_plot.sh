#!/bin/bash
gtf="/data1/CluHar/PacBio/CCS_data/MHC/analysis_for_paper/RagTag_annotations/gtf/"
out="/data1/CluHar/PacBio/CCS_data/MHC/genomic_organization/regions/RagTag/"

mkdir -p ${out}gtf
### prepare gtf files ###
samples="Reference CS2_hap1 CS2_hap2 CS4_hap1 CS4_hap2 CS5_hap1 CS5_hap2 CS7_hap1 CS7_hap2 CS8_hap1 CS8_hap2 CS10_hap1 CS10_hap2 BS1_hap1 BS2_hap1 BS3_hap1 BS1_hap2 BS2_hap2 BS3_hap2 BS4_hap1 BS5_hap1 BS6_hap1 BS4_hap2 BS5_hap2 BS6_hap2 NSSH2_hap1 NSSH2_hap2 NSSH10_hap1 NSSH10_hap2"

for i in ${samples}; do
    grep -e "A1." -e "B1." ${gtf}$i.gtf >> ${out}gtf/Locus1.gtf
    grep -e "A2." -e "B2." ${gtf}$i.gtf >> ${out}gtf/Locus2.gtf
    grep -e "A3." -e "B3." ${gtf}$i.gtf >> ${out}gtf/Locus3.gtf
    grep -e "A4." -e "B4." ${gtf}$i.gtf >> ${out}gtf/Locus4.gtf
    grep -e "A5." -e "B5." ${gtf}$i.gtf >> ${out}gtf/Locus5.gtf
    grep -e "A6." -e "B6." ${gtf}$i.gtf >> ${out}gtf/Locus6.gtf
    grep -e "A7." -e "B7." ${gtf}$i.gtf >> ${out}gtf/Locus7.gtf
    grep -e "A8." -e "B8." ${gtf}$i.gtf >> ${out}gtf/Locus8.gtf
    grep -e "A9." -e "B9." ${gtf}$i.gtf >> ${out}gtf/Locus9.gtf
done;

wait;

###### prepare files for plotting in gggenome ###
mkdir -p ${out}GeneTrack

regions="Locus1 Locus2 Locus3 Locus4 Locus5 Locus6 Locus7 Locus8 Locus9"

for i in ${regions}; do
        ## get gene positions for GeneTrack
        cat ${out}gtf/${i}.gtf | awk 'BEGIN{OFS=","}{if($3=="gene") print $1, $4, $5, $7, $10;}' | sed 's/"//g;s/;//g' | sed '1 i\seq_id,start,end,strand,gene\' > ${out}GeneTrack/${i}_positions.csv
done
