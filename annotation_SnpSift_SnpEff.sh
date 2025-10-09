#!/bin/bash
# ----------------QSUB Parameters-----------------
#PBS -q q2
#PBS -M dbrkic@bioinfo.hr
#PBS -m n
#PBS -N SnpSift_SnpEff
#PBS -l select=ncpus=16:mem=200g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #

#renaming chrM to chrMT
bcftools annotate snps_recal_all_xym.vcf.gz --rename-chrs chr_name_conv.txt -o snps_recal_all_xym_renamedM.vcf

bgzip -c snps_recal_all_xym_renamedM.vcf > snps_recal_all_xym_renamedM.vcf.gz
tabix -p vcf snps_recal_all_xym_renamedM.vcf.gz

CLINVAR=clinvar.vcf.gz
DBSNP=dbSNP/00-All.vcf.gz

#snpSift
SnpSift annotate -id ${DBSNP} \
-Xmx20G \
snps_recal_all_xym_renamedM.vcf.gz > snps_ann_dbsnp.vcf

#snpSift
SnpSift annotate -id ${CLINVAR} \
-Xmx20G \
snps_ann_dbsnp.vcf > snps_ann_both.vcf

#snpEff
snpEff eff -Xmx20G -v GRCh38.p14 \
snps_ann_both.vcf > snps_ann_all.vcf

bgzip snps_ann_dbsnp.vcf
bgzip snps_ann_both.vcf
bgzip snps_ann_all.vcf



