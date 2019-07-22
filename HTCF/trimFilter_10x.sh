#!/bin/bash
# Maps reads derived from a 10x single cell calling 
# cards library and identifies insertions. 

#SBATCH --job-name=trimFilter_10x
#SBATCH --output=logs/trimFilter_10x_%a.out # standard out goes here
#SBATCH --error=logs/trimFilter_10x_%a.err # standard error goes here
#SBATCH --array=3-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=16000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amoudgil@wustl.edu

# Initialize variables
PROJECT_DIR=/scratch/rmlab/misc/2802_K562-HCT-HyPB-noTSO-2xBiotin
PROJECT_OUT=$PROJECT_DIR"/output_and_analysis"
PROJECT_RAW=$PROJECT_DIR"/raw"
SAMPLE_LINE=$( sed -n ${SLURM_ARRAY_TASK_ID}p 10x_trimSheet.csv )

# Initialize important global variables
ORIGINAL_STEM=
NEW_STEM=

# Read local variables from $SAMPLE_LINE in subshell,
# then pass them along to the global versions
while IFS=, read -r original_stem new_stem
do
    { ORIGINAL_STEM=$original_stem; }
    { NEW_STEM=$new_stem; }
done <<< "$SAMPLE_LINE"

# Load the necessary modules
module load cutadapt

# Set variables for old and new read1 and read2 files
OLD_READ1=$PROJECT_RAW/$ORIGINAL_STEM"_R1_001.fastq.gz"
OLD_READ1_TRIM1=$PROJECT_RAW/$ORIGINAL_STEM"-trim1_R1_001.fastq.gz"
OLD_READ1_TRIM2=$PROJECT_RAW/$ORIGINAL_STEM"-trim2_R1_001.fastq.gz"

OLD_READ2=$PROJECT_RAW/$ORIGINAL_STEM"_R2_001.fastq.gz"
OLD_READ2_TRIM1=$PROJECT_RAW/$ORIGINAL_STEM"-trim1_R2_001.fastq.gz"
OLD_READ2_TRIM2=$PROJECT_RAW/$ORIGINAL_STEM"-trim2_R2_001.fastq.gz"

NEW_READ1=$PROJECT_OUT/$NEW_STEM"_S1_L001_R1_001.fastq.gz"
NEW_READ2=$PROJECT_OUT/$NEW_STEM"_S1_L001_R2_001.fastq.gz"

PB="^GGTTAA"
SB="^TGTA"

# Filter on, and trim, leading transposase sequence in read2's.
cutadapt \
    -g $PB \
    -o $OLD_READ2_TRIM1 \
    -p $OLD_READ1_TRIM1 \
    --minimum-length 1 \
    --discard-untrimmed \
    -e 0 \
    --no-indels \
    $OLD_READ2 \
    $OLD_READ1

# Trim trailing P7 adapter, if found.
cutadapt \
    -a "AGAGACTGGCAAGTACACGTCGCACTCACCATGANNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG" \
    -o $OLD_READ2_TRIM2 \
    -p $OLD_READ1_TRIM2 \
    --minimum-length 1 \
    $OLD_READ2_TRIM1 \
    $OLD_READ1_TRIM1

# Move files to new directory; remove intermediate files
cp $OLD_READ1_TRIM2 $NEW_READ1
cp $OLD_READ2_TRIM2 $NEW_READ2
rm $OLD_READ1_TRIM1
rm $OLD_READ1_TRIM2
rm $OLD_READ2_TRIM1
rm $OLD_READ2_TRIM2
