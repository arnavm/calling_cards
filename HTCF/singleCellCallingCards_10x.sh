#!/bin/bash
# Maps reads derived from a 10x single cell calling 
# cards library and identifies insertions. 

#SBATCH --job-name=scCC_10x
#SBATCH --output=logs/scCC_10x_%a.out # standard out goes here
#SBATCH --error=logs/scCC_10x_%a.err # standard error goes here
#SBATCH --array=3
#SBATCH --cpus-per-task=16
#SBATCH --mem=32000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amoudgil@wustl.edu

# Initialize variables
PROJECT_DIR=/scratch/rmlab/1/calling_card_mammalian/4302_snow_sccc_031820_pilot
PROJECT_OUT=$PROJECT_DIR"/output_and_analysis"
PROJECT_RAW=$PROJECT_DIR"/raw"
SAMPLE_LINE=$( sed -n ${SLURM_ARRAY_TASK_ID}p 10x_sampleSheet.csv )

# Initialize important global variables
LANE=
SAMPLE_ID=
INDEX=
SAMPLE_PROJECT=

# Read local variables from $SAMPLE_LINE in subshell,
# then pass them along to the global versions
while IFS=, read -r lane sample_id index sample_project
do
    { LANE=$lane; }
    { SAMPLE_ID=$sample_id; }
    { INDEX=$index; }
    { SAMPLE_PROJECT=$sample_project; }
done <<< "$SAMPLE_LINE"

# GENOME can be one of hg38, mm10, or hg38_and_mm10
GENOME="mm10"

# N_BARCODES should be set to the expected number of cells
# Note that this does not have any practical implications,
# as the count of true cells will be obtained from the scRNA-seq library.
N_BARCODES=30000

# CHEMISTRY should be set to the assay configuration. 
CHEMISTRY="SC3Pv3"

# TRANSPOSASE should be set to "PB" for piggyBac
TRANSPOSASE="PB"

# Load the necessary modules
module load samtools
module load cellranger
module load bedops

# Map reads to the genome
cellranger count \
    --id=$SAMPLE_PROJECT"_scCC_map_"$GENOME \
    --fastqs=$PROJECT_OUT \
    --transcriptome=/scratch/ref/rmlab/cellranger/custom/$GENOME \
    --sample=$SAMPLE_ID \
    --expect-cells=$N_BARCODES \
    --nosecondary \
    --chemistry=$CHEMISTRY \
    --localcores=16 \
    --localmem=30

# Filter only primary alignments (unique reads)
samtools view \
    -b -h -F 260 \
    -o $PROJECT_OUT/$SAMPLE_ID"_uniq.bam" \
    $SAMPLE_PROJECT"_map_scCC/outs/possorted_genome_bam.bam"

# Annotate transposition sites
python AnnotateInsertionSites.py \
    --transposase $TRANSPOSASE \
    -f \
    $SAMPLE_PROJECT"_scCC_map_"$GENOME"/outs/possorted_genome_bam.bam" \
    /scratch/ref/rmlab/novoalign_indexes/$GENOME/$GENOME.2bit \
    $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unfiltered.bam"

## Convert to calling card format
#python BamToCallingCard.py -b CB -i $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unfiltered.bam" -o $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unfiltered.ccf"

# Filter for transpositions supported by ≥2 different UMIs
python UMIFilter.py \
    -p 10x \
    -i $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unfiltered.bam" \
    -o $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_final.bam"

# Convert to calling card format
python BamToCallingCard.py -b CB -i $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_final.bam" -o $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unsorted.ccf"

# Sort output file
grep -v '_' $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unsorted.ccf" | sort-bed - > $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_final.ccf"

#rm $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unfiltered.ccf"
rm $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unsorted.ccf"
