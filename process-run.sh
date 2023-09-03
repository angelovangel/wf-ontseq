#! /usr/bin/env bash

# required: faster, parallel, nextflow, docker

# cat, compress, rename fastq files from a runfolder based on the samplesheet from the Shiny app
# run epi2me-labs/wf-clone-validation for every user in the samplesheet
# arg1 - csv file
# arg2 - path to fastq_pass

# output - everything goes in the results directory

# checks
if [[ $# -ne 2 ]]; then
    echo "Two parameters needed: path/to/csv and path/to/fastq_pass" >&2
    exit 2
fi

if [[ ! -f ${1} ]] || [[ ! -d ${2} ]]; then
    echo "File ${1} or ${2} does not exist" >&2
    exit 2
fi

[ -d results ] && \
echo "Directory results exists, will be deleted ..." && \
rm -rf results
mkdir -p results/fastq
mkdir -p results/wf-clone-validation
#exit 0


# get col index as they are not very consistent
user_idx=$(head -1 ${1} | sed 's/,/\n/g' | nl | grep 'user' | cut -f 1)
size_idx=$(head -1 ${1} | sed 's/,/\n/g' | nl | grep 'dna_size' | cut -f 1)
samplename_idx=$(head -1 ${1} | sed 's/,/\n/g' | nl | grep 'sample' | cut -f 1)
barcode_idx=$(head -1 ${1} | sed 's/,/\n/g' | nl | grep 'barcode' | cut -f 1)

# make the samplesheet (headers are alias, barcode and approx_size)
# 

while IFS="," read line; do
    # check if dir exists and do the work - merge/rename, make 1 samplesheet per user
    userid=$(echo $line | cut -f $user_idx -d,)
    barcode=$(echo $line | cut -f $barcode_idx -d,)
    samplename=$(echo $line | cut -f $samplename_idx -d,)
    dna_size=$(echo $line | cut -f $size_idx -d,)
    currentdir=${2}/${barcode// /}
    # skip if barcode is NA or is not a valid barcode name
    if [[ $barcode != barcode[0-9][0-9] ]]; then
        #echo "skipping ${barcode}"
        continue
    fi

    [ -d $currentdir ] && 
    [ "$(ls -A $currentdir)" ] &&
    # generate 1 samplesheet per user
    echo "${barcode},${samplename},${dna_size}"  >> results/wf-clone-validation/$userid-samplesheet.csv && \
    # generate 1 sizesheet per user
    # echo "${samplename},${dna_size}"  >> results/wf-clone-validation/$userid-sizesheet.csv && \
    echo "merging ${samplename}-${barcode}" && \
    cat $currentdir/*.fastq.gz > results/fastq/$samplename.fastq.gz || \
    echo folder ${currentdir} not found!
done < "$1"

# add headers
for f in results/wf-clone-validation/*-samplesheet.csv; do
    printf "%s\n" 1 i "barcode,alias,approx_size" . w | ed $f > /dev/null
done

# for f in results/wf-clone-validation/*-sizesheet.csv; do
#     printf "%s\n" 1 i "alias,approx_size" . w | ed $f > /dev/null
# done

# get fastq stats for the merged files
[ "$(ls -A results/fastq/)" ] && \
echo "Running faster ..." && parallel faster -ts ::: results/fastq/* > results/fastq-stats.tsv || \
echo "No fastq files found"

echo "Merging fastq done, starting the epi2me-labs/wf-clone-validation pipeline..."
exit 1
for i in results/wf-clone-validation/*samplesheet.csv; do 
    nextflow run epi2me-labs/wf-clone-validation \
    --fastq $2 \
    --sample_sheet $i \
    --approx_size_sheet results/wf-clone-validation/$(basename $i -samplesheet.csv)-sizesheet.csv \
    --out_dir results/wf-clone-validation/$(basename $i -samplesheet.csv)-assembly; \
done