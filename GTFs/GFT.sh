#!/bin/bash

# hisat2Index.sh
#
#
#module load hisat

HomeDir=$(dirname `readlink -f $0`)
cd  ${HomeDir}

############## download GTF #############
# hg38
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
gunzip -c gencode.v38.annotation.gtf.gz >gencode.v38.annotation.gtf
rm gencode.v38.annotation.gtf.gz

# mm10
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz
gunzip -c gencode.vM10.annotation.gtf.gz >gencode.vM10.annotation.gtf
rm gencode.vM10.annotation.gtf.gz
