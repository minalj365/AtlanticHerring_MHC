# Convert hifiasm assemblies to RagTag scaffolded assemblies
#!/bin/sh
PacBio="/data1/CluHar/PacBio/CCS_data/denovo/hifiasm/assemblies/Renamed/"
out="/data1/CluHar/PacBio/CCS_data/denovo/hifiasm/RagTag/"
ref="/data1/CluHar/Genomes/Assembly_v2.0.2/Ch_v2.0.2.fasta"
#samples="CS2_hap1 CS4_hap1 CS5_hap1 CS7_hap1 CS8_hap1 CS10_hap1 F1_hap1 F2_hap1 F3_hap1 F4_hap1 F5_hap1 F6_hap1 CS2_hap2 CS4_hap2 CS5_hap2 CS7_hap2 CS8_hap2 CS10_hap2 F1_hap2 F2_hap2 F3_hap2 F4_hap2 F5_hap2 F6_hap2"
samples="CS7_hap1 CS8_hap1 CS10_hap1 CS7_hap2 CS8_hap2 CS10_hap2 F4_hap1 F5_hap1 F6_hap1 F4_hap2 F5_hap2 F6_hap2"
for i in ${samples}; do
ragtag.py scaffold -r -t 2 -u -o ${out}$i ${ref} ${PacBio}$i.fa &
done;
