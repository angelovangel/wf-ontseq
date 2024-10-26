#! /usr/bin/env bash

# requires fasterplot to be installed and in env
# arg1 is path to fastq file
# output is maxbin [integer]

seqkit seq -M 49999 -g $1 | fasterplot -l - |  grep "# maxbin:" | cut -f2