#!/bin/sh
## This script first get all files (alignment and tree files) together at one place for .codeml to execute
## It edits the codeml.ctl file for each run
## It then runs codeml.
## Prerequisite: codeml.ctl general file that can be edited and codeml.exe file that can be copied

pal2nal="/data1/CluHar/PacBio/CCS_data/MHC/dNdS/PAML/pal2nal/"
trees="/data1/CluHar/PacBio/CCS_data/MHC/dNdS/PAML/trees/"
out="/data1/CluHar/PacBio/CCS_data/MHC/dNdS/PAML/codeml/"

genes="DAA1.1 DAA1.2 DAA1.3 DAA2.1 DAA2.2 DAA2.3 DAA4.1 DAA4.2 DAA5.1 DAA7.1 DAB1.1 DAB1.2 DAB1.3 DAB2.1 DAB2.2 DAB2.3 DAB4.1 DAB4.2 DAB5.1a DAB5.1b DAB7.1 DBA3.1 DBA4.1 DBB2.1 DBB3.1 DBB4.1 DBB6.1"

exons="E2 E3"
for i in ${genes}; do
        for j in ${exons}; do
                mkdir ${out}${i}_${j}
                cp ${pal2nal}${i}_${j}.nuc ${out}${i}_${j}
                cp ${trees}${i}_${j}_CDS.tre ${out}${i}_${j}
                cp ${out}codeml.exe ${out}${i}_${j}
                cp /home/minal03/bin/paml4.9j/bin/codeml ${out}${i}_${j}
                cp ${out}codeml.ctl ${out}${i}_${j}
                sed -i -e "s/pal2nal_file/${i}_${j}.nuc/"  -e "s/tree_file/${i}_${j}_CDS.tre/" -e "s/output_file/${i}_${j}.out/" ${out}${i}_${j}/codeml.ctl
        done;
done;

wait;

for i in ${genes}; do
        for j in ${exons}; do
                cd ${out}${i}_${j} && ./codeml &
        done
done
