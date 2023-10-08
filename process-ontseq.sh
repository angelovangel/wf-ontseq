#! /usr/bin/env bash

# required: faster, parallel, nextflow, docker, (faster-report.R - optional)

# cat, compress, rename fastq files from a runfolder based on the samplesheet from the ONT rapid Shiny app
# run epi2me-labs/wf-clone-validation or wf-bacterial-genomes (de novo assembly) for every user in the samplesheet
# output - everything goes in a results-ontseq/userid folder

set -e
usage="$(basename "$0") [-c SAMPLESHEET] [-p FASTQ_PASS] [-h] [-r]

Process ONT plasmid sequencing run - cat, compress, rename fastq files from a fastq_pass folder
based on the samplesheet from the Shiny app and run epi2me-labs/wf-clone-validation for every user in the samplesheet. 
    -h  show this help text
    -c  (required) samplesheet.csv, downloaded from the ONT rapid barcoding Shiny app. 
        Alternatively, a custom csv sample sheet with columns 'user', 'sample', 'dna_size' and 'barcode'
    -p  (required) path to ONT fastq_pass folder
    -w  (optional) ONT workflow to run, can be 'plasmid' or 'genome'. If unset 'plasmid' will be used.
    -r  (optional flag) generate faster-report html file"

REPORT=false;
#SUBSAMPLE=false;

options=':hrsc:p:w:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    c) SAMPLESHEET=$OPTARG;;
    p) FASTQ_PASS=$OPTARG;;
    w) WORKFLOW=$OPTARG;;
    r) REPORT=true;;
    #s) SUBSAMPLE=true;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

# set default to 'plasmid'
WORKFLOW=${WORKFLOW:-plasmid}

# mandatory arguments
if [ ! "$SAMPLESHEET" ] || [ ! "$FASTQ_PASS" ]; then
  echo "arguments -c and -p must be provided"
  echo "$usage" >&2; exit 1
fi
# get execution directory
EXECDIR=$(dirname $(readlink -f "$0"))

if [[ ! -f ${SAMPLESHEET} ]] || [[ ! -d ${FASTQ_PASS} ]]; then
    echo "File ${SAMPLESHEET} or directory ${FASTQ_PASS} does not exist" >&2
    exit 1
fi

[ -d results-ontseq ] && \
echo -e "results-ontseq folder exists, will be deleted ...\n====================" && \
rm -rf results-ontseq
mkdir -p results-ontseq
#exit 0


# get col index as they are not very consistent
user_idx=$(head -1 ${SAMPLESHEET} | sed 's/,/\n/g' | nl | grep 'user' | cut -f 1)
size_idx=$(head -1 ${SAMPLESHEET} | sed 's/,/\n/g' | nl | grep 'dna_size' | cut -f 1)
samplename_idx=$(head -1 ${SAMPLESHEET} | sed 's/,/\n/g' | nl | grep 'sample' | cut -f 1)
barcode_idx=$(head -1 ${SAMPLESHEET} | sed 's/,/\n/g' | nl | grep 'barcode' | cut -f 1)

# check samplesheet is valid
num='[0-9]+'
if [[ ! $user_idx =~ $num ]] || [[ ! $size_idx =~ $num ]] || [[ ! $samplename_idx =~ $num ]] || [[ ! $barcode_idx =~ $num ]]; then
    echo "Samplesheet is not valid, check that columns 'user','sample','dna_size','barcode' exist" >&2
    exit 1
fi

while IFS="," read line || [ -n "$line" ]; do
    # check if dir exists and do the work - merge/rename, make 1 samplesheet per user
    userid=$(echo $line | cut -f $user_idx -d, | awk '{$1=$1};1')
    barcode=$(echo $line | cut -f $barcode_idx -d,)
    samplename=$(echo $line | cut -f $samplename_idx -d, | awk '{$1=$1};1' | tr -s '[:blank:]' '[\-*]')
    dna_size=$(echo $line | cut -f $size_idx -d,)
    currentdir=${FASTQ_PASS}/${barcode// /}
    # skip if barcode is NA or is not a valid barcode name. Also skip if there is no user or sample name specified
    if [[ $barcode != barcode[0-9][0-9] ]] || [[ $userid == 'NA' ]] || [[ $samplename == 'NA' ]]; then
        #echo "skipping ${barcode}"
        continue
    fi

    [ -d $currentdir ] && 
    [ "$(ls -A $currentdir)" ] &&
    # generate 1 samplesheet per user
    mkdir -p results-ontseq/$userid
    echo "${barcode},${samplename},${dna_size}"  >> results-ontseq/$userid/samplesheet.csv && \
    echo "merging ${samplename}-${barcode}" && \
    mkdir -p results-ontseq/$userid/fastq && \
    cat $currentdir/*.fastq.gz > results-ontseq/$userid/fastq/$samplename.fastq.gz || \
    #cat $currentdir/*.fastq.gz > results-ontseq/fastq/$samplename.fastq.gz || \
    echo folder ${currentdir} not found!
done < "$SAMPLESHEET"
echo -e "===================="

# add headers
for f in results-ontseq/*/samplesheet.csv; do
    printf "%s\n" 1 i "barcode,alias,approx_size" . w | ed $f > /dev/null
done

# # SUBSAMPLE =================================================
# # downsample to given cov based on the dna_size in the samplesheet - experimentally verified that higher cov leads to wrong assemblies
# # this directory has the same structure as fastq_pass

# if [[ $SUBSAMPLE == 'true' ]] && [[ $(command -v faster2) ]]; then
#     echo -e "Subsampling reads down to 100x...\n===================="
#     #mkdir -p results-ontseq/filtered_fastq
#     # read sample sheet and downsample down to 100x based on dna_size
#     while IFS="," read line; do
#         userid=$(echo $line | cut -f $user_idx -d,)
#         barcode=$(echo $line | cut -f $barcode_idx -d,)
#         samplename=$(echo $line | cut -f $samplename_idx -d,)
#         dna_size=$(echo $line | cut -f $size_idx -d,)
#         currentdir=${FASTQ_PASS}/${barcode// /}
#         # skip if something is wrong
#         if [[ $barcode != barcode[0-9][0-9] ]] || [[ $userid == 'NA' ]] || [[ $samplename == 'NA' ]] || [[ ! $dna_size =~ $num ]]; then
#             echo "skipping ${barcode}"
#             continue
#         fi
#         total_bases=$(cat $currentdir/*.fastq.gz | faster2 -l - | paste -s -d+ - | bc) &&
#         subs_fraction=$(echo $dna_size*100/$total_bases | bc -l) ||
#         echo "Could not calculate subsample fraction for ${barcode}"
#         if [[ $subs_fraction > 1 ]]; then subs_fraction=1; fi
#         echo "$barcode -- subsample fraction: $subs_fraction" &&
#         mkdir -p results-ontseq/_downsampled_fastq_pass/$barcode &&
#         faster --sample $subs_fraction results-ontseq/$userid/fastq/$samplename.fastq.gz > results-ontseq/_downsampled_fastq_pass/$barcode/$samplename.fastq
#     done < "$SAMPLESHEET"
#     #exit 1

# fi
#============================================================================

# optionally get faster-report for the merged files
if [[ $REPORT == 'true' ]] && [[ $(command -v faster-report.R) ]]; then
    for i in results-ontseq/*/fastq; do
        [ "$(ls -A $i)" ] &&
        echo -e "Running faster-report.R in $i\n====================" &&
        faster-report.R -p $i &&
        mv faster-report.html $(dirname $i)/faster-report.html ||
        echo "No fastq files found"
    done
fi

# run faster stats
for i in results-ontseq/*/fastq; do
    nsamples=$(ls -A $i | wc -l)
    [ "$(ls -A $i)" ] && \
    echo -e "Running faster on $nsamples samples in $i...\n====================" && 
    parallel -k faster -ts ::: $i/* > $(dirname $i)/fastq-stats.tsv || 
    echo "No fastq files found"
done

if [ $WORKFLOW == 'plasmid' ]; then
    pipeline='wf-clone-validation'
elif 
    [ $WORKFLOW == 'genome' ]; then
    pipeline='wf-bacterial-genomes'
else
    echo "Use either '-w plasmid' or '-w genome'"; exit 1; 
fi

echo -e "Starting the epi2me-labs/${pipeline} pipeline...\n===================="
#exit 1

# set the CPUs and memory settings depending on where this is executed
if [ $(uname) == 'Linux' ];then 
    myconfig='prod.config'
    threads=12
else 
    myconfig='dev.config'
    threads=2
fi

# if [[ $SUBSAMPLE == 'true' ]]; then
#     FASTQ_PASS="results-ontseq/_downsampled_fastq_pass"
#     echo -e "Will use ${FASTQ_PASS} for assembly\n===================="
# fi

for i in results-ontseq/*/samplesheet.csv; do 
    echo -e "Starting $WORKFLOW assembly for $(dirname $i)\n===================="
    nextflow run epi2me-labs/${pipeline} \
    --fastq $FASTQ_PASS \
    --sample_sheet $i \
    --out_dir $(dirname $i)/assembly \
    --threads $threads \
    -c $EXECDIR/$myconfig
done


rm -rf work
echo -e "====================\nwf-ontseq finished successfully!"