#!/bin/bash
#
#SBATCH --job-name=bulkRNA_CC
#SBATCH --output=logs/bulkRNA_CC_%a.out
#SBATCH --error=logs/bulkRNA_CC_%a.err
#SBATCH --array=2-13
#SBATCH --cpus-per-task=8
#SBATCH --mem=16000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amoudgil@wustl.edu

# Initialize settings
PROJECT_DIR=/scratch/rmlab/misc/2484_HCT116_HyPBase
PROJECT_OUT=$PROJECT_DIR"/output_and_analysis"
PROJECT_RAW=$PROJECT_DIR"/raw"
SAMPLE_LINE=$( sed -n ${SLURM_ARRAY_TASK_ID}p manifest.csv )

# Initialize important global variables
FILENAME=
BARCODE=
INDEX1=
INDEX2=

# Read local variables from $SAMPLE_LINE in subshell,
# then pass them along to the global versions
while IFS=, read -r filename barcode index1 index2
do
    { FILENAME=$filename; }
    { BARCODE=$barcode; }
    { INDEX1=$index1; }
    { INDEX2=$index2; }
done <<< "$SAMPLE_LINE"

# Demultiplexes files. Needs a file of input filenames
# and a corresponding file of barcodes
BASE=`basename $FILENAME`
STEM=$BARCODE"_""${BASE%%.*}"

# Prepare output variables
OUT_STEM=$PROJECT_OUT/$STEM
OUT_SAM=$OUT_STEM".sam"
OUT_BAM=$OUT_STEM".bam"
OUT_MAP_SORT_PREFIX=$OUT_STEM"_map_sort"
OUT_BAM_MAP=$OUT_MAP_SORT_PREFIX".bam"
OUT_BED=$OUT_MAP_SORT_PREFIX".bed"
OUT_BEDGRAPH=$OUT_MAP_SORT_PREFIX".bedgraph"

# Variables to store transposon TR sequences
SB_ITR="TAAGTGTATGTAAACTTCCGACTTCAACTGTA"
ALT_SB_ITR="AAGTGTATGTAAACTTCCGACTTCAACTGTA"
PB_LTR="TTTACGCAGACTATCTTTCTAGGGTTAA"
LONG_PB_LTR="GCGTCAAT"$PB_LTR

# Transposase should be one of PB (piggyBac) or SB (SleepingBeauty)
TRANSPOSASE=PB
# Genome can be one of hg38 or mm10
GENOME=hg38

# Load the necessary modules
module load cutadapt
module load novoalign
module load samtools
module load bedops

# Trim and demultiplex by TR and insertion site (exact match)
cutadapt \
    -g "^"$BARCODE$PB_LTR \
    -o $OUT_STEM".temp.fastq.gz" \
    --minimum-length 1 \
    --discard-untrimmed \
    -e 0 \
    --no-indels \
    $PROJECT_RAW/$FILENAME

# Trim any trailing Nextera adapters (allowing mismatches)
cutadapt \
    -a "CTGTCTCTTATACACATCTCCGAGCCCACGAGACTNNNNNNNNNNTCTCGTATGCCGTCTTCTGCTTG" \
    -o $OUT_STEM".fastq.gz" \
    --minimum-length 1 \
    $OUT_STEM".temp.fastq.gz"

# Align the reads
novoalign \
    -d /scratch/ref/rmlab/novoalign_indexes/$GENOME/$GENOME.nvx \
    -f $OUT_STEM".fastq.gz" \
    -n 40 \
    -o SAM \
    -o SoftClip > $OUT_SAM

# Convert to BAM
samtools view \
    -bS -h \
    $OUT_SAM \
    -o $OUT_BAM

# Filter only mapped reads, convert to BAM, and sort
samtools view \
    -bS -h -F 260 \
    $OUT_SAM | \
samtools sort - -o $OUT_BAM_MAP

# Tag reads with barcodes
python TagBam.py \
    --tag XP:Z:$BARCODE \
    $OUT_BAM_MAP \
    $OUT_MAP_SORT_PREFIX"_tagged.bam"

# Tag reads wih i7 indexes
python TagBam.py \
    --tag XJ:Z:$INDEX1 \
    $OUT_MAP_SORT_PREFIX"_tagged.bam" \
    $OUT_MAP_SORT_PREFIX"_tagged2.bam"

# Tag reads with i5 indexes
python TagBam.py \
    --tag XK:Z:$INDEX2 \
    $OUT_MAP_SORT_PREFIX"_tagged2.bam" \
    $OUT_MAP_SORT_PREFIX"_tagged3.bam"

# Tag reads with transposon insertions
python AnnotateInsertionSites.py \
    --transposase $TRANSPOSASE \
    -f \
    $OUT_MAP_SORT_PREFIX"_tagged3.bam" \
    /scratch/ref/rmlab/novoalign_indexes/$GENOME/$GENOME.2bit \
    $OUT_MAP_SORT_PREFIX"_final.bam"

# Index the final BAM file
samtools index $OUT_MAP_SORT_PREFIX"_final.bam"

# Get an (unsorted) list of unique insertions
python BamToCallingCard.py -b XP XJ -i $OUT_MAP_SORT_PREFIX"_final.bam" -o $OUT_MAP_SORT_PREFIX"_unsorted.ccf"

# Sort the CCF file
sort-bed $OUT_MAP_SORT_PREFIX"_unsorted.ccf" > $OUT_MAP_SORT_PREFIX"_final.ccf"

# Clean up
rm $OUT_SAM
rm $OUT_BAM
rm $OUT_STEM".temp.fastq.gz"
rm $OUT_BAM_MAP
rm $OUT_MAP_SORT_PREFIX"_tagged.bam"
rm $OUT_MAP_SORT_PREFIX"_tagged2.bam"
rm $OUT_MAP_SORT_PREFIX"_tagged3.bam"
rm $OUT_MAP_SORT_PREFIX"_unsorted.ccf"
