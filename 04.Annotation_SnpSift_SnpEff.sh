#!/bin/bash
# ============================================================
# Script: 04.Annotation_SnpSift_SnpEff.sh
# Author: Doris Repusic
# Description:
#   PBS job submission script for variant annotation:
#     1. Rename mitochondrial chromosome (chrM → chrMT)
#     2. Annotate SNPs with dbSNP and ClinVar using SnpSift
#     3. Annotate effects using SnpEff
#     4. Compress and index final VCFs
#
# Usage:
#   qsub 04.Annotation_SnpSift_SnpEff.sh
#
# Requirements:
#   - bcftools
#   - SnpSift
#   - SnpEff (database: GRCh38.p14 or equivalent)
#   - bgzip + tabix (htslib)
#
# Input:
#   - snps_recal_all_xym.vcf.gz (VCF from previous pipeline step)
#   - chr_name_conv.txt (chromosome rename file, chrM → chrMT)
#   - dbSNP VCF
#   - ClinVar VCF
# ============================================================

# ---------------- PBS Directives ----------------
# Queue
#PBS -q q2
# Job name
#PBS -N SnpSift_SnpEff
# Resources
#PBS -l select=ncpus=16:mem=200g
# Merge stdout and stderr
#PBS -j oe

# Move to working directory
cd $PBS_O_WORKDIR

# ---------------- Configurable Parameters ----------------
# Reference annotation databases
CLINVAR=clinvar.vcf.gz
DBSNP=/dbSNP/00-All.vcf.gz

# Memory for SnpEff/SnpSift
JAVA_MEM="-Xmx20G"

# Input VCF (from variant calling pipeline)
INPUT_VCF=snps_recal_all_xym.vcf.gz

# Chromosome rename file (chrM → chrMT)
CHR_RENAME=chr_name_conv.txt

# ---------------- Chromosome Rename ----------------
bcftools annotate $INPUT_VCF \
    --rename-chrs $CHR_RENAME \
    -o snps_recal_all_xym_renamedM.vcf

# Compress and index renamed VCF
bgzip -c snps_recal_all_xym_renamedM.vcf > snps_recal_all_xym_renamedM.vcf.gz
tabix -p vcf snps_recal_all_xym_renamedM.vcf.gz

# ---------------- Annotation with SnpSift ----------------
SnpSift annotate -id $DBSNP $JAVA_MEM \
    snps_recal_all_xym_renamedM.vcf.gz > snps_ann_dbsnp.vcf

SnpSift annotate -id $CLINVAR $JAVA_MEM \
    snps_ann_dbsnp.vcf > snps_ann_both.vcf

# ---------------- Annotation with SnpEff ----------------
snpEff eff $JAVA_MEM -v GRCh38.p14 \
    snps_ann_both.vcf > snps_ann_all.vcf

# ---------------- Compress Outputs ----------------
bgzip snps_ann_dbsnp.vcf
bgzip snps_ann_both.vcf
bgzip snps_ann_all.vcf

# ---------------- Notes ----------------
# Input:
#   snps_recal_all_xym.vcf.gz
#
# Intermediate:
#   snps_recal_all_xym_renamedM.vcf.gz (chrM renamed → chrMT)
#
# Outputs:
#   snps_ann_dbsnp.vcf.gz  (annotated with dbSNP)
#   snps_ann_both.vcf.gz   (annotated with dbSNP + ClinVar)
#   snps_ann_all.vcf.gz    (final annotation with SnpEff)
# ============================================================
