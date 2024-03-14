#! /usr/bin/env bash

# required: faster, parallel, nextflow, docker, (faster-report.R, minimap2, samtools, perbase - optional)

# cat, compress, rename fastq files from a runfolder based on the samplesheet from the ONT rapid Shiny app
# run epi2me-labs/wf-clone-validation or wf-bacterial-genomes (de novo assembly) for every user in the samplesheet
# output - everything goes in a results-ontseq/userid folder

set -e
usage="$(basename "$0") -c SAMPLESHEET -p FASTQ_PASS [OPTIONS]

Process ONT plasmid sequencing run - cat, compress, rename fastq files from a fastq_pass folder
based on the samplesheet from the Shiny app and run epi2me-labs/wf-clone-validation for every user in the samplesheet. 
    -h  show this help text
    -c  (required) samplesheet.csv (or xlsx), downloaded from the ONT rapid barcoding Shiny app. 
        Alternatively, a custom csv/xlsx sample sheet with columns 'user', 'sample', 'dna_size' and 'barcode'
    -p  (required) path to ONT fastq_pass folder
    -w  (optional) ONT workflow to run, can be 'plasmid', 'genome' or 'amplicon'. If unset 'plasmid' will be used.
    -r  (optional flag) generate faster-report html file
    -s  (optional flag) use singularity profile (docker by default)
    -m  (optional flag) do mapping of reads to assembly after the wf-clone-validation pipeline (for plasmids only)
    -t  (optional flag) zip results (per user) and transfer to a webdav endpoint. For this to work, env variables have to be specified"

REPORT=false;
SINGULARITY=false;
MAPPING=false;
TRANSFER=false;

function logmessage () {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] - $1 \n================================================"
}

function timestamp () {
    date +'%Y%m%d-%H%M%S'
}

options=':hrsmtc:p:w:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    c) SAMPLESHEET=$OPTARG;;
    p) FASTQ_PASS=$OPTARG;;
    w) WORKFLOW=$OPTARG;;
    r) REPORT=true;;
    s) SINGULARITY=true;;
    m) MAPPING=true;;
    t) TRANSFER=true;;
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

if [[ ! -f ${SAMPLESHEET} ]] || [[ ! -d ${FASTQ_PASS} ]]; then
    echo "File ${SAMPLESHEET} or directory ${FASTQ_PASS} does not exist" >&2
    exit 1
fi
# get execution directory
EXECDIR=$(dirname $(readlink -f "$0"))

# which workflow and check if correct
if [ $WORKFLOW == 'plasmid' ]; then
    pipeline='wf-clone-validation'
elif 
    [ $WORKFLOW == 'genome' ]; then
    pipeline='wf-bacterial-genomes'
elif
    [ $WORKFLOW == 'amplicon' ]; then
    pipeline='wf-amplicon'
else
    echo "Use either '-w plasmid' or '-w genome' or '-w amplicon'"; echo "You used $WORKFLOW"; exit 1; 
fi
# setup results directory
RESULTS=$(dirname $FASTQ_PASS)/results-ontseq-$WORKFLOW

[ -d $RESULTS ] && \
logmessage "results-ontseq folder exists, will be deleted ..." && \
rm -rf $RESULTS
mkdir -p $RESULTS

# convert to csv if excel is provided, from here on $csvfile is used
infile_ext=${SAMPLESHEET##*.}
if [ $infile_ext == 'xlsx' ]; then
    logmessage 'Excel file provided, will be converted to csv ...'
    excel2csv.R $SAMPLESHEET && # writes the csv to the same location as the excel file
    csvfile=$(dirname $SAMPLESHEET)/$(basename $SAMPLESHEET .$infile_ext).csv && 
    logmessage "CSV file generated ==> ${csvfile}" ||
    logmessage 'Converting Excel to csv failed...!'
else
    logmessage 'CSV file provided...'
    csvfile=$SAMPLESHEET
fi


# get col index as they are not very consistent
user_idx=$(head -1 ${csvfile} | sed 's/,/\n/g' | nl | grep 'user' | cut -f 1)
size_idx=$(head -1 ${csvfile} | sed 's/,/\n/g' | nl | grep 'dna_size' | cut -f 1)
samplename_idx=$(head -1 ${csvfile} | sed 's/,/\n/g' | nl | grep 'sample' | cut -f 1)
barcode_idx=$(head -1 ${csvfile} | sed 's/,/\n/g' | nl | grep 'barcode' | cut -f 1)

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
    mkdir -p $RESULTS/$userid
    echo "${barcode},${samplename},${dna_size}"  >> $RESULTS/$userid/samplesheet.csv && \
    echo "merging ${samplename}-${barcode}" && \
    mkdir -p $RESULTS/$userid/01-fastq && \
    cat $currentdir/*.fastq.gz > $RESULTS/$userid/01-fastq/$samplename.fastq.gz || \
    echo folder ${currentdir} not found!
done < "$csvfile"
echo -e "===================="

# add headers
for f in $RESULTS/*/samplesheet.csv; do
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
    for i in $RESULTS/*/01-fastq; do
        [ "$(ls -A $i)" ] &&
        logmessage "Running faster-report.R for $i" &&
        faster-report.R -p $i &&
        mv faster-report.html $(dirname $i)/faster-report.html ||
        echo "No fastq files found"
    done
fi

# run faster stats
for i in $RESULTS/*/01-fastq; do
    nsamples=$(ls -A $i | wc -l)
    [ "$(ls -A $i)" ] && \
    logmessage "Running faster on $nsamples samples in $i..." &&
    parallel -k faster -ts ::: $i/* > $(dirname $i)/fastq-stats.tsv || 
    echo "No fastq files found"
done


logmessage "Starting the epi2me-labs/${pipeline} pipeline..."
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

# enable singularity
if [ $SINGULARITY == 'true' ]; then
    profile='singularity'
else
    profile='standard'
fi

# run once for every user
for i in $RESULTS/*/samplesheet.csv; do
    logmessage "Starting $WORKFLOW assembly for $(dirname $i)"
    nextflow run epi2me-labs/${pipeline} \
    --fastq $FASTQ_PASS \
    --sample_sheet $i \
    --out_dir $(dirname $i)/02-assembly \
    --threads $threads \
    -c $EXECDIR/$myconfig \
    -profile $profile
done

rm -rf work
logmessage "$WORKFLOW assembly finished successfully!"

function mapper() {
    queryname=$(basename $2 | cut -d. -f1)
    minimap2 -t $threads -ax map-ont $1 $2 > $3/$queryname.sam
    samtools view -S -b -@ $threads -T $1 $3/$queryname.sam | \
    samtools sort -@ $threads -o $3/$queryname.bam -
    samtools index -@ $threads $3/$queryname.bam
    rm $3/$queryname.sam
    perbase base-depth $3/$queryname.bam -F 260 > $3/$queryname.perbase.tsv
    # take only primary alignments, the flag is htslib thing
    #https://github.com/sstadick/perbase/issues/68
    # perbase only-depth $3/$queryname.bam > $3/$queryname.depth.tsv
    # get problem positions
    $EXECDIR/bin/perpos_freq.sh $3/$queryname.perbase.tsv > $3/$queryname.problems.tsv
}

# do mapping of reads to assembly 
# do once for every user and sample
if [ $MAPPING == 'true' ]; then
    # outer loop - per user
    for i in $RESULTS/*; do 
        user=$(basename $i)
        logmessage "Starting mapping reads to assembly for user $user..."
        mkdir -p $RESULTS/$user/03-mapping
        mapping_output=$RESULTS/$user/03-mapping
        # inner loop - per sample
        [ $(ls -A $RESULTS/$user/02-assembly/*.final.fasta) ] && # only go here if assembly produced something
        for j in $RESULTS/$user/02-assembly/*.final.fasta; do
            k=$(basename $j .final.fasta)
            gbk=$RESULTS/$user/02-assembly/$k.annotations.gbk
            query=$RESULTS/$user/01-fastq/$k.fastq.gz
            cov=$mapping_output/$k.perbase.tsv
            problems=$mapping_output/$k.problems.tsv

            mapper $j $query $mapping_output
            logmessage "Generating coverage plot for $k"
            #echo -e "Generating coverage plot for $k"
            $EXECDIR/bin/plot_plasmid.py $gbk $cov $problems
        done
    done
fi

if [ $TRANSFER == 'true' ]; then
    # load sensitive env variables USERNAME PASS and URL
    #eval "$(direnv export bash)" && direnv allow $EXECDIR
    logmessage "Compress and transfer to wahadrive..."
    #for i in $RESULTS/*; do zip -r $i $i; done
    for i in $RESULTS/*; do tar -cvf $i.tar -C $i .; done
    # make a folder on the endpoint for this analysis run
    thisrun=$(timestamp)-$WORKFLOW
    curl -u $USERNAME:$PASS -X MKCOL $URL/$thisrun && \
    for i in $RESULTS/*.tar; do curl -T $i -u $USERNAME:$PASS $URL/$thisrun/; done &&
    logmessage "Transfer finished ..." || \
    logmessage "Transfer failed ..."
fi

logmessage "wf-ontseq finished successfully!"
