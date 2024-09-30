#! /usr/bin/env bash

# required: faster, parallel, nextflow, docker, gzip, (faster-report.R, minimap2, samtools, bedtools, perbase, IGV-report - optional)

# cat, compress, rename fastq files from a runfolder based on the samplesheet from the ONT rapid Shiny app
# run epi2me-labs/wf-clone-validation or wf-bacterial-genomes (de novo assembly) for every user in the samplesheet
# output - everything goes in a results-ontseq/userid folder
# | tee -a $output_directory/logfile.log

set -e
usage="$(basename "$0") -c SAMPLESHEET -p FASTQ_PASS [OPTIONS]

Process ONT plasmid sequencing run - cat, compress, rename fastq files from a fastq_pass folder
based on the samplesheet from the Shiny app and run epi2me-labs/wf-clone-validation for every user in the samplesheet. 
    -h  show this help text
    -c  (required) samplesheet.csv (or xlsx), downloaded from the ONT rapid barcoding Shiny app. 
        Alternatively, a custom csv/xlsx sample sheet with columns 'user', 'sample', 'dna_size' and 'barcode'
    -p  (required) path to ONT fastq_pass folder
    -w  (optional) ONT workflow to run, can be 'plasmid', 'genome' or 'amplicon'. If unset 'plasmid' will be used.
    -x  (optional) config profile to use for nextflow run, can be 'dev' or 'prod' 
    -l  (optional flag) Do not filter reads by length based on approx_size parameter (plasmid workflow only).
    -r  (optional flag) generate faster-report html file
    -s  (optional flag) use singularity profile (docker by default)
    -m  (optional flag) do mapping of reads to assembly after the wf-clone-validation pipeline (for plasmids only)
    -t  (optional flag) zip results (per user) and transfer to a webdav endpoint. For this to work, env variables have to be specified
    -n  (optional) run name, default is {timestamp-workflow name}"

REPORT=false;
SINGULARITY=false;
MAPPING=false;
TRANSFER=false;
LARGE_CONSTRUCT=false;

function logmessage () {
    echo -e "================================================\n[$(date +'%Y-%m-%d %H:%M:%S')] - $1 \n================================================" \
    2>&1 | tee -a $2
}

function timestamp () {
    date +'%Y%m%d-%H%M%S'
}

options=':hrsmtln:c:p:w:x:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    c) SAMPLESHEET=$OPTARG;;
    p) FASTQ_PASS=$OPTARG;;
    w) WORKFLOW=$OPTARG;;
    x) CONFIG=$OPTARG;;
    l) LARGE_CONSTRUCT=true;;
    n) RUNNAME=$OPTARG;;
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

# set default nxf config to 'prod'
CONFIG=${CONFIG:-prod}

thisrun=$(timestamp)-$WORKFLOW

# set default run name to timestamp-workflow if not provided as parameter
RUNNAME=${RUNNAME:-$thisrun}
# setup results directory
RESULTS=$(dirname $FASTQ_PASS)/$RUNNAME
mkdir -p $RESULTS
#touch $RESULTS/$RUNNAME.log

logmessage "Starting run $RUNNAME" $RESULTS/$RUNNAME.log

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


# convert to csv if excel is provided, from here on $csvfile is used
infile_ext=${SAMPLESHEET##*.}
if [ $infile_ext == 'xlsx' ]; then
    #logmessage 'Excel file provided, will be converted to csv ...'
    ./bin/excel2csv.R $SAMPLESHEET && # writes the csv to the same location as the excel file
    csvfile_orig=$(dirname $SAMPLESHEET)/$(basename $SAMPLESHEET .$infile_ext).csv && 
    logmessage "CSV file generated from Excel --> ${csvfile_orig}" $RESULTS/$RUNNAME.log ||
    logmessage 'Converting Excel to csv failed...!' $RESULTS/$RUNNAME.log
else
    #logmessage 'CSV file provided...'
    csvfile_orig=$SAMPLESHEET
fi

# validate samplesheet 
logmessage "Validating samplesheet --> $csvfile_orig" $RESULTS/$RUNNAME.log
cat $csvfile_orig > $RESULTS/samplesheet-original.csv
./bin/validate_samplesheet.R $csvfile_orig $FASTQ_PASS $RESULTS 2>&1 | tee -a $RESULTS/$RUNNAME.log

csvfile=$RESULTS/samplesheet-validated.csv
logmessage "Validated samplesheet --> $RESULTS/samplesheet-validated.csv" $RESULTS/$RUNNAME.log
logmessage "Problematic samples in samplesheet (OK if empty):" $RESULTS/$RUNNAME.log
awk -F, '$11 !~ /OK/' $RESULTS/samplesheet-validated.csv >> $RESULTS/$RUNNAME.log

# get col index as they are not very consistent
user_idx=$(head -1 ${csvfile} | sed 's/,/\n/g' | nl | grep 'user' | cut -f 1)
size_idx=$(head -1 ${csvfile} | sed 's/,/\n/g' | nl | grep 'dna_size' | cut -f 1)
samplename_idx=$(head -1 ${csvfile} | sed 's/,/\n/g' | nl | grep 'sample' | cut -f 1)
barcode_idx=$(head -1 ${csvfile} | sed 's/,/\n/g' | nl | grep 'barcode' | cut -f 1)

# check samplesheet indexes are valid
num='[0-9]+'
if [[ ! $user_idx =~ $num ]] || [[ ! $size_idx =~ $num ]] || [[ ! $samplename_idx =~ $num ]] || [[ ! $barcode_idx =~ $num ]]; then
    echo "Samplesheet is not valid, check that columns 'user','sample','dna_size','barcode' exist" >&2
    exit 1
fi

# check samplenames are unique (for the whole run)
# samplesheet validation done by validate_samplesheet.R
#
logmessage "Running merge-rename for $RUNNAME" $RESULTS/$RUNNAME.log
while IFS="," read line || [ -n "$line" ]; do
    # check if dir exists and do the work - merge/rename, make 1 samplesheet per user
    userid=$(echo $line | cut -f $user_idx -d, | awk '{$1=$1};1')
    barcode=$(echo $line | cut -f $barcode_idx -d,)
    samplename=$(echo $line | cut -f $samplename_idx -d, | awk '{$1=$1};1' | tr -s '[:blank:]' '[\-*]')
    dna_size=$(echo $line | cut -f $size_idx -d,)
    currentdir=${FASTQ_PASS}/${barcode// /}
    # skip if barcode is NA or is not a valid barcode name. Also skip if there is no user or sample name specified
    if [[ $barcode != barcode[0-9][0-9] ]] || \
    [[ $userid == 'NA' ]] || \
    [[ $samplename == 'NA' ]] || \
    #[[ ! -d $currentdir ]] || \ 
    [[ ! "$(ls -A $currentdir/*.fastq.gz)" ]]; then
        echo "== skipping ${barcode} ==" 2>&1 | tee -a $RESULTS/$RUNNAME.log
        continue
    fi

    # generate 1 samplesheet per user
    mkdir -p $RESULTS/$userid
    echo "${barcode},${samplename},${dna_size}"  >> $RESULTS/$userid/samplesheet.csv && \
    echo "merging ${samplename}-${barcode} for $userid" 2>&1 | tee -a $RESULTS/$RUNNAME.log && \
    mkdir -p $RESULTS/$userid/01-fastq && \
    cat $currentdir/*.fastq.gz > $RESULTS/$userid/01-fastq/$samplename.fastq.gz || \
    echo folder ${currentdir} not found!
done < "$csvfile"
echo -e "===================="

# get run info from fastq header to use in faster-report
# it is the same for all users
fastqfiles=($FASTQ_PASS/barcode*/*.fastq.gz)
FLOWCELL=$(gzip -cd ${fastqfiles[1]} | head -n 1 | grep -oE "flow_cell_id=.*" | cut -d" " -f1 | cut -d= -f2)
RUNDATE=$(gzip -cd ${fastqfiles[1]} | head -n 1 | grep -oE "start_time=.*" | cut -d" " -f1 | cut -d= -f2 | cut -dT -f1)
BC_MODEL=$(gzip -cd ${fastqfiles[1]} | head -n 1 | grep -oE "model_version_id=.*" | cut -d" " -f1 | cut -d= -f2 | cut -dT -f1)

if [[ -z "$FLOWCELL" ]]; then
    FLOWCELL="NA"
fi

if [[ -z "$RUNDATE" ]]; then
    RUNDATE="NA"
fi

if [[ -z "$BC_MODEL" ]]; then
    BC_MODEL="NA"
fi

# add headers to peruser samplesheet
for f in $RESULTS/*/samplesheet.csv; do
    printf "%s\n" 1 i "barcode,alias,approx_size" . w | ed $f > /dev/null
done

# optionally get faster-report for the merged files
if [[ $REPORT == 'true' ]] && [[ $(command -v faster-report.R) ]]; then
    for i in $RESULTS/*/01-fastq; do
        currentuser=$(basename $(dirname $i))
        nsamples=$(ls -A $i | wc -l | tr -d ' ')
        [ "$(ls -A $i)" ] &&
        logmessage "Running faster-report.R on $nsamples samples for $currentuser" $RESULTS/$RUNNAME.log &&
        #faster-report-docker.sh -p $(realpath $i) -d $RUNSTART -f $FLOWCELL &&
        faster-report.R -p $i --rundate $RUNDATE --flowcell $FLOWCELL --user $currentuser --basecall $BC_MODEL &&
        mv faster-report.html $(dirname $i)/$currentuser-faster-report.html \
        || echo "faster-report.html generation error for $currentuser" 2>&1 | tee -a $RESULTS/$RUNNAME.log
    done
fi

# run faster stats
for i in $RESULTS/*/01-fastq; do
    currentuser=$(basename $(dirname $i))
    nsamples=$(ls -A $i | wc -l | tr -d ' ')
    echo -e "file\treads\tbases\tn_bases\tmin_len\tmax_len\tN50\tGC_percent\tQ20_percent" > $(dirname $i)/$currentuser-fastq-stats.tsv
    [ "$(ls -A $i)" ] && \
    logmessage "Running faster on $nsamples samples for $currentuser" $RESULTS/$RUNNAME.log &&
    parallel -k faster2 -ts ::: $i/* >> $(dirname $i)/$currentuser-fastq-stats.tsv \
    || echo "No fastq files found"
done

# set the CPUs and memory settings depending on where this is executed
if [ $CONFIG == 'prod' ];then 
    myconfig='prod.config'
    threads=12
else 
    myconfig='dev.config'
    threads=2
fi


# enable singularity
if [ $SINGULARITY == 'true' ]; then
    profile='singularity'
else
    profile='standard'
fi

# do not filter reads by len for plasmid wf
if [ $LARGE_CONSTRUCT == 'true' ]; then
    largeconstruct='--large_construct'
else
    largeconstruct=""
fi

# exit 0
# run once for every user
for i in $RESULTS/*/samplesheet.csv; do
    currentuser=$(basename $(dirname $i))
    nsamples=$(tail -n +2 $i | wc -l | tr -d ' ')
    logmessage "Starting $WORKFLOW assembly on $nsamples samples for $currentuser " $RESULTS/$RUNNAME.log
    nextflow run epi2me-labs/${pipeline} \
        --fastq $FASTQ_PASS \
        --sample_sheet $i \
        --out_dir $(dirname $i)/02-assembly \
        --threads $threads \
        $largeconstruct \
        -c $EXECDIR/$myconfig \
        -profile $profile
    # move-rename report 
    mv $(dirname $i)/02-assembly/*-report.html $(dirname $i)/$currentuser-$WORKFLOW-assembly-report.html \
    || logmessage "No assembly report found for $currentuser" $RESULTS/$RUNNAME.log
    logmessage "$WORKFLOW assembly finished for $currentuser" $RESULTS/$RUNNAME.log
done

rm -rf work
logmessage "$WORKFLOW assembly finished for all users" $RESULTS/$RUNNAME.log

# inputs are final.fasta, fastq, output folder
function mapper() {
    queryname=$(basename $2 | cut -d. -f1)
    minimap2 -t $threads -ax map-ont $1 $2 > $3/$queryname.sam
    samtools view -S -b -@ $threads -T $1 $3/$queryname.sam | \
    samtools sort -@ $threads -o $3/$queryname.bam -
    samtools index -@ $threads $3/$queryname.bam
    rm $3/$queryname.sam
    perbase base-depth --threads 4 $3/$queryname.bam -F 260 > $3/$queryname.perbase.tsv
    # take only primary alignments, the flag is htslib thing
    #https://github.com/sstadick/perbase/issues/68
    # perbase only-depth $3/$queryname.bam > $3/$queryname.depth.tsv
    # get problem positions
    $EXECDIR/bin/perpos_freq.sh $3/$queryname.perbase.tsv > $3/$queryname.problems.tsv
}

# do mapping of reads to assembly 
# do once for every user and sample
if [ $MAPPING == 'true' ] && [ $WORKFLOW != 'amplicon' ]; then
    # outer loop - per user
    for i in $RESULTS/*/; do #directories (user) only
        user=$(basename $i)
        logmessage "Starting mapping reads to assembly for user $user..." $RESULTS/$RUNNAME.log
        [ $(ls $RESULTS/$user/02-assembly/*.final.fasta | wc -l) -gt 0 ] && mkdir -p $RESULTS/$user/03-mapping
        [ $(ls $RESULTS/$user/02-assembly/*.final.fasta | wc -l) -gt 0 ] && mkdir -p $RESULTS/$user/04-igv-reports
        mapping_output=$RESULTS/$user/03-mapping
        igvoutput=$RESULTS/$user/04-igv-reports
        # inner loop - per sample
        [ $(ls $RESULTS/$user/02-assembly/*.final.fasta | wc -l) -gt 0 ] && # only go here if assembly produced something
        for j in $RESULTS/$user/02-assembly/*.final.fasta; do
            k=$(basename $j .final.fasta)
            bed=$RESULTS/$user/02-assembly/$k.annotations.bed
            query=$RESULTS/$user/01-fastq/$k.fastq.gz
            cov=$mapping_output/$k.perbase.tsv
            problems=$mapping_output/$k.problems.tsv
            bam=$mapping_output/$k.bam
            len=$(faster2 -l $j)
            header=$(grep ">" $j | cut -c 2-)

            mapper $j $query $mapping_output
            logmessage "Generating IGV report for $user --- $k" $RESULTS/$RUNNAME.log
            # dynamic calculation for subsampling, subsample for > 500 alignments
            count=$(samtools view -c $bam)
            subsample=$(echo $count | awk '{if ($1 <200) {print 1} else {print 200/$1}}')
            #echo -e "Generating coverage plot for $k"
            #$EXECDIR/bin/plot_plasmid.py $gbk $cov $problems
            
            $EXECDIR/bin/fix_bed.R $bed $len > $igvoutput/annotations.bed || logmessage "Could not fix bed file!" $RESULTS/$RUNNAME.log
            echo -e "$header\t0\t$len\tsubsampled alignments (fraction $subsample)" > $igvoutput/bedfile.bed
            awk -v OFS='\t' -v chr=$k 'NR>1 {print chr, $1-1, $1, "HET"}' $problems >> $igvoutput/bedfile.bed
            
            create_report \
                $igvoutput/bedfile.bed \
                --fasta $j \
                --tracks $igvoutput/annotations.bed $bam \
                --subsample $subsample \
                --output $igvoutput/$k-igvreport.html \
            || logmessage "Create IGV report failed for $k!" $RESULTS/$RUNNAME.log
            #rm $igvoutput/bedfile.bed
            #rm $igvoutput/annotations.bed
        done || logmessage "No assemblies found for $user" $RESULTS/$RUNNAME.log
    done
fi

if [ $TRANSFER == 'true' ]; then
    # load sensitive env variables USERNAME PASS and URL
    #eval "$(direnv export bash)" && direnv allow $EXECDIR
    logmessage "Compress and transfer to wahadrive..." $RESULTS/$RUNNAME.log
    #for i in $RESULTS/*; do tar -cvf $i.tar -C $i .; done
    for i in $RESULTS/*/01-fastq; do folder=$(dirname $i); tar -cvf $folder.tar -C $folder .; done #only tar user folders
    # make a folder on the endpoint for this analysis run
    curl --fail -u $USERNAME:$PASS -X MKCOL $URL/$RUNNAME && \
    for i in $RESULTS/*.tar; do curl --fail -T $i -u $USERNAME:$PASS $URL/$RUNNAME/; done &&
    logmessage "Transfer finished OK!" $RESULTS/$RUNNAME.log || \
    logmessage "Transfer failed!" $RESULTS/$RUNNAME.log
fi

logmessage "$RUNNAME finished successfully!" $RESULTS/$RUNNAME.log
