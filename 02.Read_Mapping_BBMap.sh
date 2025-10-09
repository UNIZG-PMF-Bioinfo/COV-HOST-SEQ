#!/bin/bash
# ============================================================
# Script: 01.Read_Mapping_BBMap.sh
# Author: Doris Repusic
# Description: 
#   PBS job submission script for mapping paired-end reads 
#   using BBMap, generating BAM output and mapping statistics.
#
#   This script is intended for use on an HPC cluster with PBS.
#   It supports job arrays for processing multiple input files.
#
# Usage:
#   qsub -J 0-N 01.Read_Mapping_BBMap.sh
#   (where N is the number of files - 1)
#
# Requirements:
#   - BBMap installed and accessible
#   - Samtools installed (optional, if you want to post-process BAMs)
#   - Input files in INPUT_DIR ending with *_1.phixtrim.fq.gz
# ============================================================

# ---------------- PBS Directives ----------------
# Queue selection
#PBS -q q2
# Job name
#PBS -N BBmap_trimmed
# Number of CPUs, memory, and threads per node
#PBS -l select=ncpus=16:mem=256g:ompthreads=24
# Job array (adjust range based on number of samples)
#PBS -J 0-200
# Merge stdout and stderr
#PBS -j oe

# Move to the directory where the job was submitted
cd $PBS_O_WORKDIR

# ---------------- Configurable Parameters ----------------
MEMORY=256g          # Memory allocation for BBMap
THREADS=16           # Number of threads to use
INPUT_DIR=/Trimming  # Directory containing input reads

# Paths to software
SAMTOOLS=/samtools
BBDIR=/BBmap/bbmap
BBMAP=bbmap.sh

# Extra parameters for BBMap
BBMAP_PAR="pigz=t local=t"

# ---------------- File Selection ----------------
# Collect input files (_1.fq.gz paired reads after trimming)
IN_SEQ=($INPUT_DIR/*_1.phixtrim.fq.gz)

# Select the file corresponding to the current array index
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}

# Derive the sample base name (removes path and suffix)
BASE=${FILE##*/}
BASE=${BASE%%_*}

echo $BASE

# Input file pair placeholder (# replaces 1/2)
INFILE="${FILE/_1/_#}"

# ---------------- Run BBMap ----------------
# Generate BAM output + statistics
$BBDIR/$BBMAP \
    -Xmx$MEMORY \
    threads=$THREADS \
    $BBMAP_PAR \
    in1=${INFILE} \
    out=${BASE}.bam \
    rgid="$BASE" \
    rgsm="$BASE" \
    scafstats=${BASE}_scafstats.txt \
    mhist=${BASE}_mhist.txt \
    bhist=${BASE}_bhist.txt

# ---------------- Notes ----------------
# - Output BAM file: <sample>.bam
# - Scaffold stats:  <sample>_scafstats.txt
# - Mismatch hist:   <sample>_mhist.txt
# - Base hist:       <sample>_bhist.txt
# ============================================================

# ---------------- Post-processing with Samtools ----------------
# Convert gzipped SAM to BAM
zcat ${BASE}.sam.gz | $SAMTOOLS view -@ $THREADS -Sb - > ${BASE}.bam

# Sort BAM file
$SAMTOOLS sort -@ $THREADS ${BASE}.bam -o ${BASE}_sorted.bam
$SAMTOOLS index ${BASE}_sorted.bam

# Remove unsorted BAM if sorted version exists
[ -f "${BASE}_sorted.bam" ] && rm -f ${BASE}.bam

# Remove duplicates
$SAMTOOLS rmdup ${BASE}_sorted.bam ${BASE}_sorted_dedupe.bam
$SAMTOOLS index ${BASE}_sorted_dedupe.bam

# ---------------- Coverage track generation ----------------
# Generate bedGraph from BAM
genomeCoverageBed -ibam ${BASE}_sorted.bam -bg -split -g $CHRLEN > ${BASE}.bedGraph

# Convert bedGraph to BigWig
wigToBigWig ${BASE}.bedGraph $CHRLEN ${BASE}.bw

# Remove bedGraph if BigWig file was created
[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
