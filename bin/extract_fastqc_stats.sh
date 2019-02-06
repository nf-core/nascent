#!/bin/sh

# Extract read statistics from fastqc output, providing a stopping point if we
# don't have enough.

# Set options
set -e
set -o errexit
#set -o pipefail
#	Assume only	utf-8
export LC_ALL=C

usage()
{
		echo "extract_nascent_stats.sh - extract fastqc stats from output zip files"
		echo "Example:"
		echo "    ./extract_nascent_stats.sh --srr=SRR2084556"
		echo "Usage:"
		echo "    -h/--help -- Display this help message."
		echo "    --srr     -- SRR / FASTA to parse"
		exit 0
}

while [ "$1" != "" ]; do
		PARAM=$(echo "$1" | awk -F= '{print $1}')
		VALUE=$(echo "$1" | awk -F= '{print $2}')
		case $PARAM in
				-h | --help)
						usage
						exit
						;;
				--srr)
						SRR=$VALUE
						;;
				*)
						echo "ERROR: unknown parameter \"$PARAM\""
						usage
						exit 1
						;;
		esac
		shift
done
echo Extracting fastqc statistics for "$SRR"

GC=$(unzip -c "$(find . -name *_fastqc.zip)" "$SRR"_fastqc/fastqc_data.txt \
         | grep "%GC" | grep -o "[0-9]*")
SEQ=$(unzip -c "$(find . -name *_fastqc.zip)" "$SRR"_fastqc/fastqc_data.txt | \
          grep "Total Sequences" | \
          grep -o "[0-9]*")
DEDUP=$(unzip -c "$(find . -name *_fastqc.zip)" "$SRR"_fastqc/fastqc_data.txt | \
            grep "#Total Deduplicated Percentage" | \
            grep -o "[0-9,.]*")

echo -e "SRR\t%GC\tTotal_Sequences\t%Total_Deduplicated"
echo -e "$SRR""$(printf "\\t")""$GC""$(printf "\\t")""$SEQ""$(printf "\\t")""$DEDUP"
