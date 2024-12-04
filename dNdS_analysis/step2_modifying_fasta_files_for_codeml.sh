#!/bin/sh
alignments="/data1/CluHar/PacBio/CCS_data/MHC/dNdS/PAML/alignments/"
out="/data1/CluHar/PacBio/CCS_data/MHC/dNdS/PAML/pal2nal/"
pal2nal="/home/minal03/bin/pal2nal.v14/pal2nal.pl"

genes="DAA1.1 DAA1.2 DAA1.3 DAA2.1 DAA2.2 DAA2.3 DAA4.1 DAA4.2 DAA5.1 DAA7.1 DBA3.1 DBA4.1 DAB1.1 DAB1.2 DAB1.3 DAB2.1 DAB2.2 DAB2.3 DAB4.1 DAB4.2 DAB5.1a DAB5.1b DAB7.1 DBB2.1 DBB3.1 DBB4.1 DBB6.1"
exons="E2 E3"

for i in ${genes}; do
        for j in ${exons}; do
                ${pal2nal} ${alignments}${i}_${j}_AA.fa ${alignments}${i}_${j}_CDS.fa -output paml -nogap > ${out}${i}_${j}.nuc
        done;
done;
