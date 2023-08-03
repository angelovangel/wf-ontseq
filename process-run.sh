#! /usr/bin/env bash

# required: faster, parallel, nextflow, docker

# cat, compress, rename fastq files from a runfolder based on the samplesheet from the Shiny app
# run epi2me-labs/wf-clone-validation for every user in the samplesheet
# arg1 - csv file
# arg2 - path to fastq_pass

# output - everything goes in the results directory

# get col index as they are not very consistent
user_idx=$(head -1 ${1} | sed 's/,/\n/g' | nl | grep 'user' | cut -f 1)
size_idx=$(head -1 ${1} | sed 's/,/\n/g' | nl | grep 'dna_size' | cut -f 1)
samplename_idx=$(head -1 ${1} | sed 's/,/\n/g' | nl | grep 'sample' | cut -f 1)
barcode_idx=$(head -1 ${1} | sed 's/,/\n/g' | nl | grep 'barcode' | cut -f 1)

# make the samplesheet and size sheet file (headers are alias, barcode and alias, approx_size)
# echo "alias,barcode" > results/wf-clone-validation/samplesheet.csv

[ -d results/fastq ] && \
echo "Directory results/fastq exists, will be deleted ..." && \
rm -rf results/fastq
mkdir -p results/fastq
#exit 0
[ -d results/wf-clone-validation ] && \
echo "Directory results/wf-clone-validation exists, will be deleted ..." && \
rm -rf results/wf-clone-validation
mkdir -p results/wf-clone-validation

while IFS="," read line
do
# check if dir exists and do the work - merge/rename, make 1 samplesheet per user
userid=$(echo $line | cut -f $user_idx -d,)
barcode=$(echo $line | cut -f $barcode_idx -d,)
samplename=$(echo $line | cut -f $samplename_idx -d,)
currentdir=${2}/${barcode// /}

[ -d $currentdir ] && 
# generate 1 samplesheet per user
echo "${samplename},${barcode}"  >> results/wf-clone-validation/$userid-samplesheet.csv && \
# generate 1 sizesheet per user
echo "${samplename},${barcode}"  >> results/wf-clone-validation/$userid-samplesheet.csv && \

echo "merging ${samplename}-${barcode}" && \
cat $currentdir/*.fastq.gz > results/fastq/$samplename.fastq.gz || \
echo folder ${currentdir} not found!
done < "$1"

# get fastq stats for the merged files
[ "$(ls -A results/fastq/)" ] && \
echo "Running faster ..." && parallel faster -ts ::: results/fastq/* > results/fastq-stats.tsv || \
echo "no fastq files found"

echo "Merging fastq done, will start the epi2me-labs/wf-clone-validation pipeline..."



#awk -F "," -v si="$samplename_idx" -v bc="$barcode_idx" '{print $si "," $bc}' $1 >> results/wf-clone-validation/samplesheet.csv

exit 0

for i in *samplesheet.csv; \
do nextflow run epi2me-labs/wf-clone-validation --fastq fastq_pass \
--sample_sheet $i \
--approx_size_sheet $(basename $i -samplesheet.csv)-size.csv \
--out_dir $(basename $i -samplesheet.csv)-out; \
done
