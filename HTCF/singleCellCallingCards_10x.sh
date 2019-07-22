#!/bin/bash
# Maps reads derived from a 10x single cell calling 
# cards library and identifies insertions. 

#SBATCH --job-name=scCC_10x
#SBATCH --output=logs/scCC_10x_%a.out # standard out goes here
#SBATCH --error=logs/scCC_10x_%a.err # standard error goes here
#SBATCH --array=3-15
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amoudgil@wustl.edu

# Initialize variables
PROJECT_DIR=/scratch/rmlab/misc/3193_InVitro2+Ctx_dualIndex
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
GENOME="hg38_mm10"

# N_BARCODES should be set to the expected number of cells
N_BARCODES=6000

# CHEMISTRY should be set to the assay configuration. Default is
# SC3pV2 (Single Cell 3' v1/v2)
CHEMISTRY="SC3Pv2"

# TRANSPOSASE should be one of PB (piggyBac), SB (SleepingBeauty), HelR (Helraiser)
TRANSPOSASE="PB"

# Load the necessary modules
module load samtools
module load bcl2fastq2
module load cellranger
module load cutadapt
module load bedops

for GENOME in hg38_mm10 hg38 mm10;
do
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
    --localmem=60

python TagBamWithTransposonInserts.py \
    --transposase $TRANSPOSASE \
    -f \
    $SAMPLE_PROJECT"_scCC_map_"$GENOME"/outs/possorted_genome_bam.bam" \
    /scratch/ref/rmlab/novoalign_indexes/$GENOME/$GENOME.2bit \
    $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unfiltered.bam"

python BamToCallingCard.py -b CB -i $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unfiltered.bam" -o $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unfiltered.ccf"

python scCC_dyadFilter.py \
    -p 10x \
    -c $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unfiltered.bam" \
    -o $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_final.bam"

python BamToCallingCard.py -b CB -i $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_final.bam" -o $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unsorted.ccf"

sort-bed $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unsorted.ccf" > $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_final.ccf"
rm $PROJECT_OUT/$SAMPLE_ID"_"$GENOME"_unsorted.ccf"
done
