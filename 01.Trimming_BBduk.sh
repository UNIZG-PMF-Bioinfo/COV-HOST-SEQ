#!/bin/bash
# ============================================================
# Script: 01.Trimming_BBduk.sh
# Author: Doris Repusic
# Description:
#   PBS job submission script for trimming adapters and spike-in
#   sequences from paired-end FASTQ reads using BBDuk (BBMap suite).
#
#   Supports job arrays for parallel processing of multiple samples.
#
# Usage:
#   qsub -J 0-N 01.Trimming_BBduk.sh
#   (where N is the number of unique samples - 1)
#
# Requirements:
#   - BBMap suite (BBDuk) installed and accessible
#   - Adapter FASTA file available
# ============================================================

# ---------------- PBS Directives ----------------
# Queue selection
#PBS -q q2
# Job name
#PBS -N trim_adapters.PE
# Resources: CPUs, memory per node
#PBS -l select=ncpus=12:mem=40g
# Job array (adjust index range to your dataset)
#PBS -J 0-200
# Merge stdout and stderr
#PBS -j oe

# Move to the directory where the job was submitted
cd $PBS_O_WORKDIR

# ---------------- Configurable Parameters ----------------
THREADS=12         # Number of threads
MEMORY=40g         # Memory allocation

# Input directory with FASTQ files
IN_DIR=/Data

# Adapter reference file (edit path if necessary)
ADAPTER=/adapters_wxs.fasta

# BBDuk parameters
BBDUK_PAR="overwrite=t \
ktrim=r \
k=23 \
rcomp=t \
mink=11 \
hdist=1 \
minoverlap=8 \
minlength=80 \
qtrim=lr \
trimq=25 \
tbo \
copyundefined=t \
threads=$THREADS \
-Xmx$MEMORY"

# ---------------- File Selection ----------------
# Collect input FASTQ files
IN_SEQ=($(find $IN_DIR -maxdepth 1 -name "*.fq.gz"))

# Remove _1/_2 suffixes to create unique sample list
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_*fq.gz}" | sort -u))

# Select current sample based on job array index
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${IN_DIR}/}

echo $BASE

# ---------------- Run BBDuk ----------------
# Remove spike-in reads and adapter contamination
bbduk.sh \
    in1=${FILE}_1.fq.gz \
    in2=${FILE}_2.fq.gz \
    out1=${BASE}_1.phixtrim.fq.gz \
    out2=${BASE}_2.phixtrim.fq.gz \
    ref="artifacts,phix" \
    ref=$ADAPTER \
    stats=${BASE}.stats \
    $BBDUK_PAR 2> ${BASE}.trim.log

# ---------------- Notes ----------------
# - Input:  <sample>_1.fq.gz / <sample>_2.fq.gz
# - Output: <sample>_1.phixtrim.fq.gz / <sample>_2.phixtrim.fq.gz
# - Stats:  <sample>.stats
# - Log:    <sample>.trim.log
# ============================================================
