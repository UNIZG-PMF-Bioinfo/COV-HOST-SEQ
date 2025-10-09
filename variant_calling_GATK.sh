#!/bin/bash
# ----------------QSUB Parameters-----------------
#PBS -q q2
#PBS -M dbrkic@bioinfo.hr
#PBS -m n
#PBS -N gatk_varcalling
#PBS -l select=ncpus=16:mem=400g
#PBS -J 0-200
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #

REF=/hg38.fa

INPUT_DIR=.

# input reads
IN_SEQ=($INPUT_DIR/*_sorted_dedupe_groups.bam)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*${INPUT_DIR}/}
BASE=${BASE%%_*}

# variant calling
gatk HaplotypeCaller \
        -R $REF \
        -I $FILE \
        -O ${BASE}_variants.vcf

# extracting SNPs
gatk SelectVariants \
       -R $REF \
       -V ${BASE}_variants.vcf \
       -select-type SNP \
       -O ${BASE}_raw_snps.vcf

# filtering SNPs
gatk VariantFiltration \
       -V ${BASE}_raw_snps.vcf \
       -filter "QD < 2.0" --filter-name "QD2" \
       -filter "QUAL < 30.0" --filter-name "QUAL30" \
       -filter "SOR > 3.0" --filter-name "SOR3" \
       -filter "FS > 60.0" --filter-name "FS60" \
       -filter "MQ < 40.0" --filter-name "MQ40" \
       -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
       -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
       -O ${BASE}_filtered_snps.vcf

# excluding filtered variants
gatk SelectVariants \
        --exclude-filtered \
        -V ${BASE}_filtered_snps.vcf \
        -O ${BASE}_bqsr_snps.vcf


#Base Quality Score Recalibration (BQSR)
gatk BaseRecalibrator \
       -I ${BASE}_sorted_dedupe_groups.bam \
       -R $REF \
       --known-sites ${BASE}_bqsr_snps.vcf \
       -O ${BASE}_recal_data.table

#Apply BSQR
gatk ApplyBQSR \
	-R $REF \
	-I ${BASE}_sorted_dedupe_groups.bam \
	--bqsr-recal-file ${BASE}_recal_data.table \
	-O ${BASE}_sorted_dedupe_recal.bam

# variant calling
gatk HaplotypeCaller \
        -R $REF \
        -I ${BASE}_sorted_dedupe_recal.bam \
        -O ${BASE}_recal_variants.vcf

# selecting SNPs
gatk SelectVariants \
        -R $REF  \
        -V ${BASE}_recal_variants.vcf \
        -select-type SNP \
        -O ${BASE}_recal_snps.vcf

# filtering SNPs
gatk VariantFiltration \
        -V ${BASE}_recal_snps.vcf \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O ${BASE}_recal_filtered_snps.vcf

# merging and indexing all files samples into one vcf file
bcftools merge /*_recal_filtered_snps.vcf.gz -m snps -0 -O z -o /snps_recal_all.vcf.gz
bcftools index snps_recal_all.vcf.gz

# removing SNPs from scaffolds, keeping only those from complete chromosomes
bcftools view /snps_recal_all.vcf.gz | grep -w '^#\|^#CHROM\|chr[1-9]\|chr[1-2][0-9]\|chr[X]\|chr[Y]\|chr[M]' > /snps_recal_all_xym.vcf.gz