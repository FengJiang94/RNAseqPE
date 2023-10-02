#!/bin/bash

# exit and report the failure if any command fails
#exit_trap () {
  #local lc="$BASH_COMMAND" rc=$?
  #echo "Command [$lc] exited with code [$rc]"
#}
#trap exit_trap EXIT
#set -e

################################ RNAseqPE.sh ###########################
#Get pipeline directory
HomeDir=$(dirname `readlink -f $0`)
#set default parameters
array=( "$@" ) #read all command line variable to "array"
counter=0 #counter to keep track of the variable indexes
Run=$(pwd) # current working directory
input="" # names of raw fastq files as input
adapters="adapters.fa" # A fasta file containing the adapter sequences 
trim="PE -threads 16 -phred33" #running parameters for trimmomatic
CLIP=":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:25 MINLEN:32 HEADCROP:0" # clipping parameters for trimmomatic
hisat="-p 24 -x" # running parameters for hisat
genome="genome" # the hisat index used for mapping
featureCounts="-T 10  -g gene_name -B -C -p --ignoreDup --fracOverlap 0.1" # running parameters for featureCounts
GTF="GTF.gtf" # gtf annotation for featureCounts
output=$(pwd) # output directory
trimmomatic="/software/trimmomatic/0.36/trimmomatic-0.36.jar"

if [ $# == 0 ]
    then
        echo "********************************************************************************************"
        echo "*                    RNAseqPE: pipeline for paired-end RNAseq analysis.                    *"
        echo "*                              Version 6, 2021-06-05, F.J                                  *"
        echo "* Usage: RNAseqPE.sh -i name1 name2 name3 [options]                                        *"
        echo "*        Required (paired end data) :                                                      *"
        echo "*                 -i --input names of raw fastq/fastq.gz                                   *"
        echo "*                    For multiple files use -i name1 name2 name3                           *"
        echo "*                    Paired fastq files should be named as [name1]_1.fq.gz [name1]_2.fq.gz.*"
        echo "*        Optional flags:                                                                   *"
        echo "*                 -a --adapters a fasta file containing the adapter sequences.             *"
        echo "*                 -t path to trimmomatic                                                   *"
        echo "*                    Defualt:  /software/trimmomatic/0.36/trimmomatic-0.36.jar             *"
        echo "*                 --trim running parameters passed to trimmomatic. must be in quote        *"
        echo "*                    Default: \"PE -threads 16 -phred33\"                                  *"
        echo "*                 -C --CLIP clipping parameters passed to trimmomatic. must be in quote    *"
        echo "*                    Default: \":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:25 MINLEN:32 *"
        echo "*                    HEADCROP:0\"                                                          *"
        echo "*                 -h --hisat running parameters passed to hisat. must be in quote          *"
        echo "*                    Default: \"-p 24 -x\"                                                 *"
        echo "*                 -g --genome prefix or path to the hisat index used for mapping.          *"
        echo "*                    By defualt, RNAseqPE will look for index in RNAseqPE/hisat2Index/     *"
        echo "*                    For example: hg38                                                     *"
        echo "*                    For example: home/reference/hisat2/index/hg38                         *"
        echo "*                 -f --featureCounts running parameters passed to featureCounts. in quote  *"
        echo "*                    Default: \"-T 10 -g gene_name -B -C -p --ignoreDup --fracOverlap 0.1\"*"
        echo "*                 --GTF prefix or path to the GTF annotation file for featureCounts.       *"
        echo "*                    By defualt, RNAseqPE will look for GTF in RNAseqPE/GTFs/              *"
        echo "*                    For example: hg38.gtf                                                 *"
        echo "*                    For example: home/reference/GTF/hg38.gtf                              *"
        echo "*                 -o --output the ouput directory.                                         *"
        echo "*                    Default: pwd                                                          *"
        echo "* This analysis pipeline takes in raw fastq and generate a count table for each fastq pair.*"
        echo "* Analysis includes fastqc, trimminng by trimmomatic, mapping by hisat, and featureCounts. *"
        echo "* NOTE: fastqc, hisat, samtools, subread, multiqc are required and should be loaded in the *"
        echo "* the way your cluster required.                                                           *"
        echo "********************************************************************************************"
        exit 1
fi
echo "*****************************************************************************"
echo "*            RNAseqPE: pipeline for paired-end RNAseq analysis.             *"
echo "*                      Version 6, 2021-06-05,  F.J                          *"
echo "*****************************************************************************"
echo "0. Loading softwares:"
#echo "${array[@]}"
#echo "${#array[@]}"
#echo "${#Input[@]}"
#echo "${Input[@]}"
#echo "${Name[@]}"
#echo "${#Name[@]}"

#Get pipeline directory
echo " Home Directory for RNAseqPE: $HomeDir"
#echo " ${array[@]}"
#echo " ${#array[@]}"
for var in "$@" #for every variable
do 
#if the variable is "-i", set up the name of input files
    if [[ $var == "-i" || $var == "--input" ]] 
    then 
        Input=( ${array[$counter+1]} )
        #echo " Input filess ${array[$counter+1]} "
        #echo "${Input[@]}"
        # if there is more than one input file, set up it as well
        index=$((counter+2))
        while [ $index -lt ${#array[@]} ]
        do
          if [[ ${array[$index]} != "-a" && ${array[$index]} != "--adapters" && ${array[$index]} != "-i" && ${array[$index]} != "--input" && ${array[$index]} != "-t" && ${array[$index]} != "--trim" && ${array[$index]} != "-h" && ${array[$index]} != "--hist" && ${array[$index]} != "-C" && ${array[$index]} != "--CLIP" && ${array[$index]} != "--help" && ${array[$index]} != "-g" && ${array[$index]} != "--genome" && ${array[$index]} != "-f" && ${array[$index]} != "--featureCounts" && ${array[$index]} != "--GTF" && ${array[$index]} != "-o" && ${array[$index]} != "--output" ]]
          then
              Input+=( ${array[$index]} )
              #echo " Input files ${array[$index]} "
              ((index+=1))
              #echo "${index}"
          else
              index=${#array[@]}
              #echo "${index}"
          fi
        done
        #echo " All files loaded"
#if the variable is -a, set up the adapter file
    elif [[ $var == "-a" || $var == "--adapters" ]]
    then
        adapters=( ${array[$counter+1]} )
#if the variable is -t set up the running parameters for trimmomatic
    elif [[ $var == "--trim" ]]
    then
        trim=( ${array[$counter+1]} )
#if the variable is -t set up the running parameters for trimmomatic
    elif [[ $var == "-t" ]]
    then
        trimmomatic=( ${array[$counter+1]} )
#if the variable is -C, set up the clipping parameters for trimmomatic
    elif [[ $var == "-C" || $var == "--CLIP" ]]
    then
        CLIP=( ${array[$counter+1]} )
#if the variable is -h, set up the running parameters for hisat
    elif [[ $var == "-h" || $var == "--hisat" ]]
    then
        hisat=( ${array[$counter+1]} )
#if the variable is -g, set up the hisat index for mapping
    elif [[ $var == "-g" || $var == "--genome" ]]
    then
        genome=( ${array[$counter+1]} )
#if the variable is -f, set up the parameters for featureCounts
    elif [[ $var == "-f" || $var == "--featureCounts" ]]
    then
        featureCounts=( ${array[$counter+1]} )
#if the variable is --GTF, set up the GTF annotation file for featureCounts
    elif [[ $var == "--GTF" ]]
    then
        GTF=( ${array[$counter+1]} )
#if the variable is -o, set up the ouput directory
    elif [[ $var == "-o" || $var == "--output" ]]
    then
        output=( ${array[$counter+1]} )
#if the variable is -h, generate help information
    elif [[ $var == "--help"  ]]
    then
        echo "********************************************************************************************"
        echo "*                    RNAseqPE: pipeline for paired-end RNAseq analysis.                    *"
        echo "*                              Version 6, 2021-06-05, F.J                                  *"
        echo "* Usage: RNAseqPE.sh -i name1 name2 name3 [options]                                        *"
        echo "*        Required (paired end data) :                                                      *"
        echo "*                 -i --input names of raw fastq/fastq.gz                                   *"
        echo "*                    For multiple files use -i name1 name2 name3                           *"
        echo "*                    Paired fastq files should be named as [name1]_1.fq.gz [name1]_2.fq.gz.*"
        echo "*        Optional flags:                                                                   *"
        echo "*                 -a --adapters a fasta file containing the adapter sequences.             *"
        echo "*                 -t path to trimmomatic                                                   *"
        echo "*                    Defualt:  /software/trimmomatic/0.36/trimmomatic-0.36.jar             *"
        echo "*                 --trim running parameters passed to trimmomatic. must be in quote        *"
        echo "*                    Default: \"PE -threads 16 -phred33\"                                  *"
        echo "*                 -C --CLIP clipping parameters passed to trimmomatic. must be in quote    *"
        echo "*                    Default: \":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:25 MINLEN:32 *"
        echo "*                    HEADCROP:0\"                                                          *"
        echo "*                 -h --hisat running parameters passed to hisat. must be in quote          *"
        echo "*                    Default: \"-p 24 -x\"                                                 *"
        echo "*                 -g --genome prefix or path to the hisat index used for mapping.          *"
        echo "*                    By defualt, RNAseqPE will look for index in RNAseqPE/hisat2Index/     *"
        echo "*                    For example: hg38                                                     *"
        echo "*                    For example: home/reference/hisat2/index/hg38                         *"
        echo "*                 -f --featureCounts running parameters passed to featureCounts. in quote  *"
        echo "*                    Default: \"-T 10 -g gene_name -B -C -p --ignoreDup --fracOverlap 0.1\"*"
        echo "*                 --GTF prefix or path to the GTF annotation file for featureCounts.       *"
        echo "*                    By defualt, RNAseqPE will look for GTF in RNAseqPE/GTFs/              *"
        echo "*                    For example: hg38.gtf                                                 *"
        echo "*                    For example: home/reference/GTF/hg38.gtf                              *"
        echo "*                 -o --output the ouput directory.                                         *"
        echo "*                    Default: pwd                                                          *"
        echo "* This analysis pipeline takes in raw fastq and generate a count table for each fastq pair.*"
        echo "* Analysis includes fastqc, trimminng by trimmomatic, mapping by hisat, and featureCounts. *"
        echo "* NOTE: fastqc, hisat, samtools, subread, multiqc are required and should be loaded in the *"
        echo "* the way your cluster required.                                                           *"
        echo "********************************************************************************************"
        exit 1
    fi
    ((counter+=1)) #iterate counter
done

# check each argument of the function

# check whether input files are valid
# if no input provided
if [[ ${#Input[@]} == 0 ]]
then
    echo "[Error] No input provided."
    exit 1
fi

# check whether the paired sequenced input files exist
for i in ${Input[@]}
do
 if [ -s ${i}_1.fq.gz ]
 then
     if [ -s ${i}_2.fq.gz ]
     then
         echo " load sequence files: ${i}_1.fq.gz ${i}_2.fq.gz "
     else
         echo "[Error] ${i}_2.fq.gz dose not exist. Please make sure paired sequence files are named as [name1]_1.fq.gz [name1]_2.fq.gz."
         exit 1
     fi
 else
     echo "[Error] ${i}_1.fq.gz dose not exist. Please make sure paired sequence files are named as [name1]_1.fq.gz [name1]_2.fq.gz."
     exit 1
 fi
done

# check adapter file
if [ -s $adapters ]
then
    echo " adapter file: $adapters "
else
    echo "[Error] adapter file $adapters does not exist"
    exit 1
fi

# check trimmomatic
if [ -s ${trimmomatic} ]
then
    echo " trimmomatic: ${trimmomatic} "
else
    echo "[Error] trimmomatic ${trimmomatic} does not exist"
    exit 1
fi
echo " set trimmomatic parameter ${trim[@]} "
echo " set ILLUMINACLIP:${adapters}${CLIP[@]} "

# check hist2 index file
if [ -d ${HomeDir}/hisat2Index/$genome ]
then
    genome=${HomeDir}/hisat2Index/$genome
    echo " use hisat2 index: $genome "
elif [ -d $genome ]
then
    echo " use hisat2 index: $genome "
else
    echo "[Error] hisat2 index does not exist. Please provide the correct path or choose from ${HomeDir}/hisat2Index/ "
    exit 1
fi
echo " set hisat2 ${hisat[@]} $genome/genome "

# check GTF file
if [ -s ${HomeDir}/GTFs/$GTF ]
then
    GTF=${HomeDir}/GTFs/$GTF
    echo " use GTF annotation: $GTF "
elif [ -s $GTF ]
then
    echo " use GTF annotation: $GTF "
else
    echo "[Error] GTF annotation does not exist. Please provide the correct path or choose from ${HomeDir}/GTFs/ "
    exit 1
fi
echo " set featureCounts ${featureCounts[@]} -a $GTF "

# check output directory
# check hist2 index file
if [ -d $output ]
then
    echo " results stored in $output "
else
    mkdir $output
    echo " results stored in $output "
fi

cd $outpath
# directory storing the trimmed data
if [ -d "trimmed_data" ] 
then
    echo "[Warning] Directory $outpath/trimmed_data already exists."
else
    mkdir trimmed_data
fi
# directory storing the fastqc reports before trimming
if [ -d "fastqc_before" ] 
then
    echo "[Warning] Directory $outpath/fastqc_before already exists."
else
    mkdir fastqc_before
fi
# directory storing the fastqc reports after trimming
if [ -d "fastqc_after" ] 
then
    echo "[Warning] Directory $outpath/fastqc_after already exists."
else
    mkdir fastqc_after
fi
# directory storing the mapped data 
if [ -d "mapped_data" ] 
then
    echo "[Warning] Directory $outpath/mapped_data already exists."
else
    mkdir mapped_data
fi
# directory storing the raw counts
if [ -d "raw_counts" ] 
then
    echo "[Warning] Directory $outpath/raw_counts already exists."
else
    mkdir raw_counts
fi

# time to start the analysis
start=$(date +%s.%N)
echo "***** analysis began........"

## creat report files
echo "" > trimmomatic.log
echo "" > hisat2.log
# trimming and mapping
for i in ${Input[@]}
do
 echo "***** generating fastqc report for raw ${i}........"
 o=$(echo "$i" |rev|cut -d "/" -f 1|rev)
 # generate fastqc reports before trimming
 fastqc -o ${outpath}/fastqc_before -t 16 ${i}_1.fq.gz ${i}_2.fq.gz
 if [ ! -s ${outpath}/fastqc_before/${i}_2_fastqc.html ]
 then
     echo "[Error] failed to generate fastqc report for ${i}."
     exit 1
 fi
 echo "***** trimming ${i}........"
 # trimming
 java -jar ${trimmomatic} \
 ${trim[@]} \
 ${i}_1.fq.gz ${i}_2.fq.gz \
 ${outpath}/trimmed_data/${o}_1.tmp.fq.gz \
 ${outpath}/trimmed_data/${o}_1.tmu.fq.gz  \
 ${outpath}/trimmed_data/${o}_2.tmp.fq.gz \
 ${outpath}/trimmed_data/${o}_2.tmu.fq.gz \
 ILLUMINACLIP:${adapters}${CLIP} 2>>trimmomatic.log
 if [ ! -s ${outpath}/trimmed_data/${o}_1.tmp.fq.gz ]
 then
     echo "[Error] failed to trim ${i}. Please check trimmomatic.log for details"
     exit 1
 fi
 
 echo "***** generating fastqc report for trimmed ${i}........"
 # generate fastqc reports after trimming
 fastqc -o ${outpath}/fastqc_after -t 16 ${outpath}/trimmed_data/${o}_1.tmp.fq.gz ${outpath}/trimmed_data/${o}_2.tmp.fq.gz
 
 echo "***** mapping ${i}........"
 # mapping
 hisat2 ${hisat[@]} $genome/genome -1 ${outpath}/trimmed_data/${o}_1.tmp.fq.gz -2 ${outpath}/trimmed_data/${o}_2.tmp.fq.gz -S ${outpath}/mapped_data/${o}.sam 2>>hisat2.log
 if [ ! -s ${outpath}/mapped_data/${o}.sam ]
 then
     echo "[Error] failed to map ${i}. Please check hisat2.log for details"
     exit 1
 fi
 
 # convert sam to bam (only keep the paired the reads)
 samtools view -bS --threads 10 ${outpath}/mapped_data/${o}.sam > ${outpath}/mapped_data/${o}.bam
 samtools view -bF 12 --threads 10 ${outpath}/mapped_data/${o}.bam > ${outpath}/mapped_data/${o}_mapped.bam
 samtools sort --threads 10 -m 2G -o ${outpath}/mapped_data/${o}_mapped.sort.bam ${outpath}/mapped_data/${o}_mapped.bam
 samtools index ${outpath}/mapped_data/${o}_mapped.sort.bam
 # remove the orignal sam
 rm ${outpath}/mapped_data/${o}.sam
 echo "${i} mapped"
done
echo "ALL file mapped"

# counting
cd ${outpath}/mapped_data
echo "" > featureCounts.log
for i in ${Input[@]}
do
 echo "***** counting ${i}........"
 o=$(echo "$i" |rev|cut -d "/" -f 1|rev)
 featureCounts ${featureCounts[@]} -a ${GTF} -o ${outpath}/raw_counts/${o}.txt ${o}_mapped.sort.bam 2>>featureCounts.log
 if [ ! -s ${outpath}/raw_counts/${o}.txt ]
 then
     echo "[Error] failed to count ${i}. Please check ${outpath}/mapped_data/featureCounts.log for details"
     exit 1
 fi
 echo "${i} counted"
done
echo "ALL file counted"

# generate multiqc reports
cd ${outpath}
echo "***** generating multiqc reports........"
multiqc fastqc_before -n multi_fastqc_before
if [ ! -s multi_fastqc_before.html ]
then
    echo "[Error] failed to generate multiqc report for ${i}."
    exit 1
fi
multiqc fastqc_after -n multi_fastqc_after
if [ ! -s multi_fastqc_after.html ]
then
    echo "[Error] failed to generate multiqc report for ${i}."
    exit 1
fi
# generate summary of analysis
# generate tmp files cotaining the mapping information
cd ${Run}
echo "***** generating mapping report........"
printf "%s\n" "${Input[@]}" > ${outpath}/sample.csv
grep "Input Read Pairs:" trimmomatic.log |cut -d " " -f 4 >${outpath}/raw_reads.csv
grep "Both Surviving:" trimmomatic.log |cut -d " " -f 7 >${outpath}/raw_reads.csv
grep "aligned concordantly exactly 1 time" hisat2.log |cut -d " " -f 5 >${outpath}/unique_aligned_reads.csv
grep 'Successfully assigned alignments' ${outpath}/mapped_data/featureCounts.log |cut -d ":" -f 2 |cut -d " " -f 2 >${outpath}/assigned_reads.csv
# generate a final report with all above mapping information
Samples=($(less ${outpath}/sample.csv))
RawReads=($(less ${outpath}/raw_reads.csv))
CleanReads=($(less ${outpath}/raw_reads.csv))
UalignReads=($(less ${outpath}/unique_aligned_reads.csv))
AssignReads=($(less ${outpath}/assigned_reads.csv))
# creat the report
echo -e "Sample\tRawReads\tCleanReads\tUniuqeAlignedReads\tAlignRate\tAssignedReads\tAssignRate" >${outpath}/MappingReport.txt
for i in $(seq 1 ${#Input[@]})
do
 Sample=${Samples[i]}
 Raw=${RawReads[i]}
 Clean=${CleanReads[i]}
 Ualign=${UalignReads[i]}
 Assign=${AssignReads[i]}
 UalignRate=$(echo "scale=2 ; $Ualign / $Clean" |bc)
 AssignRate=$(echo "scale=2 ; $Assign / $Clean" |bc)
 echo -e "${Sample}\t${Raw}\t${Clean}\t${Ualign}\t${UalignRate}\t${Assign}\t${AssignRate}" >>${outpath}/MappingReport.txt
done
echo "report generated"
# remove tmp files
rm ${outpath}/raw_reads.csv
rm ${outpath}/unique_aligned_reads.csv
rm ${outpath}/unique_aligned_reads.csv
rm ${outpath}/assigned_reads.csv

# time to end the analysis
Finish=$(date +%s.%N)
# execution time
dur=$(echo "$Finish - $start" | bc)
echo "analysis started at $start, ended at $Finish"
echo "Excution time $dur"

echo "********************************************************************************************"
echo "*                                   analysis report                                        *"
echo "********************************************************************************************"
cat ${outpath}/MappingReport.txt
