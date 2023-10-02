#!/bin/bash

# hisat2Index.sh
#  
#
#module load hisat

HomeDir=$(dirname `readlink -f $0`)
cd  ${HomeDir}

############## make hisat index #############

# hg38
mkdir hg38 # directory to store the index
cd hg38
# download the reference file for hg38 from gencode
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz
gunzip -c gencode.v38.transcripts.fa.gz > gencode.v38.transcripts.fa
# make hisat index
hisat2-build -f gencode.v38.transcripts.fa genome # this requires significant RAM
# finishing
rm gencode.v38.transcripts.fa.gz
rm gencode.v38.transcripts.fa
cd ..

# mm10
mkdir mm10 # directory to store the index
cd mm10
# download the reference file for hg38 from gencode
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/GRCm38.p4.genome.fa.gz
gunzip -c GRCm38.p4.genome.fa.gz > GRCm38.p4.genome.fa
# make hisat index
hisat2-build -f GRCm38.p4.genome.fa genome # this requires significant RAM
# finishing
rm GRCm38.p4.genome.fa.gz
rm GRCm38.p4.genome.fa


