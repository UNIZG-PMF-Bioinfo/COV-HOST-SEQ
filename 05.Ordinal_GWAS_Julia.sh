#!/bin/bash
###############################################################################
# Script Name: 05.Ordinal_GWAS.sh
# Author: Doris Repusic
# Description: End-to-end GWAS analysis pipeline using PLINK, R, Julia, and SnpSift
#
# Steps:
#   1. Convert VCF to PLINK format and perform QC filtering
#   2. Add/update sex and phenotype information
#   3. Perform pruning and PCA (to generate covariates)
#   4. Run Ordinal GWAS (Julia)
#   5. Visualize results (Manhattan, QQ plots in R)
#   6. Extract top SNPs and annotations from VCF
#   7. LD analysis of candidate regions
#
# HPC: Designed for PBS workload manager
###############################################################################

# ---------------- QSUB Parameters ---------------- #
#PBS -q q2
#PBS -N plink_pipeline
#PBS -l select=ncpus=15:mem=200g
#PBS -j oe

cd $PBS_O_WORKDIR

###############################################################################
# Step 1: Convert VCF to PLINK format
###############################################################################
plink --vcf snps_ann_all.vcf.gz \
    --recode \
    --keep keep_samples.txt \
    --out snps_ann

###############################################################################
# Step 2: Add and check sex information
###############################################################################
plink --file snps_ann \
    --update-sex sex_file.list \
    --make-bed \
    --out snps_ann_sex

# Impute and check sex consistency
plink --bfile snps_ann_sex \
    --impute-sex 0.964 0.969 \
    --make-bed \
    --out snps_ann_sex_imp

plink --bfile snps_ann_sex_imp \
    --check-sex \
    --out snps_ann_imp

###############################################################################
# Step 3: Add phenotype and perform QC filtering
###############################################################################
plink --bfile snps_ann_sex \
    --pheno pheno_file.list \
    --make-bed \
    --out snps_ann_sex_pheno

plink --bfile snps_ann_sex \
    --geno 0.05 \
    --hwe 0.001 \
    --keep pheno_july2025.list \
    --maf 0.05 \
    --mind 0.05 \
    --remove low_het_individuals.txt \
    --pheno pheno_july2025.list \
    --pheno-name severity \
    --make-bed \
    --out july2025_sev3_nohet_qcfiltered_maf005

# Notes:
# - Genotyping rate ~0.9999
# - 167,152 variants and 177 individuals pass QC

###############################################################################
# Step 4: LD pruning and PCA
###############################################################################
plink --bfile july2025_sev3_nohet_qcfiltered_maf005 \
    --indep 50 5 2 \
    --out july2025_sev3_nohet_qcfiltered_maf005_pruned

plink --bfile july2025_sev3_nohet_qcfiltered_maf005 \
    --extract july2025_sev3_nohet_qcfiltered_maf005_pruned.prune.in \
    --make-bed \
    --out july2025_sev3_nohet_qcfiltered_maf005_pruneddata

plink --bfile july2025_sev3_nohet_qcfiltered_maf005_pruneddata \
    --pca 20 \
    --out july2025_sev3_nohet_qcfiltered_maf005_pca

# Add header to PCA file
echo -e "FID\tIID\tPC1\tPC2\t...PC20" | \
cat - july2025_sev3_nohet_qcfiltered_maf005_pca.eigenvec \
> july2025_sev3_nohet_qcfiltered_maf005_pca_header.eigenvec

###############################################################################
# Step 5: Create covariate file (R)
###############################################################################
source activate r_env
Rscript - <<'EOF'
library(dplyr); library(magrittr); library(ggplot2)

pcs <- read.table("july2025_sev3_nohet_qcfiltered_maf005_pca_header.eigenvec",
                  header=TRUE) %>% as_tibble()
pheno <- read.table("pheno_july2025.list", header=TRUE) %>% as_tibble()
covariates <- right_join(pheno, pcs %>% dplyr::select(FID:PC10))

write.csv(covariates, "july2025_covariatesPC_177.csv", row.names=FALSE, quote=FALSE)
EOF

###############################################################################
# Step 6: Ordinal GWAS in Julia
###############################################################################
conda activate julia_env
julia -e '
using OrdinalGWAS, DataFrames
const datadir = "ind177_input"
const covfile = datadir * "/july2025_covariatesPC_177.csv"
const plkfile = datadir * "/july2025_sev3_nohet_qcfiltered_maf005"

ordinalgwas(@formula(severity ~ sex + age + PC1 + PC2 + PC3 + PC4 + PC5 +
                     PC6 + PC7 + PC8 + PC9 + PC10),
            covfile, plkfile, test=:lrt, analysistype="singlesnp")
'

###############################################################################
# Step 7: Visualization of GWAS results (R)
###############################################################################
source activate r_env
Rscript - <<'EOF'
library(qqman); library(data.table); library(dplyr); library(writexl); library(stringr)

results <- fread("ordinalgwas.pval.txt", sep=",", header=TRUE) %>%
           filter(hwepval>0.05, chr<23)

# Manhattan plot
png("Manhattan_plot.png", width=20, height=15, units="cm", res=600)
manhattan(results, chr="chr", bp="pos", snp="snpid", p="pval",
          main="Ordinal GWAS: Severity (177 individuals, covariates PC10)",
          col=c("blue","gray"),
          suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8))
dev.off()

# QQ plot
png("QQ_plot.png", res=600, width=15, height=15, units="cm")
qq(results$pval, main="QQ plot: Ordinal Severity GWAS")
dev.off()

# Export top SNPs
results_sorted <- arrange(results, pval)
write_xlsx(results_sorted[1:30,], "Top30_177_ordinal_variants.xlsx")

# Create SNP position lists for extraction
results %>%
  filter(pval<0.05) %>%
  arrange(pval) %>%
  select(chr, pos) %>%
  mutate(chr=case_when(chr==23~"X", chr==24~"Y", TRUE~as.character(chr)),
         chr=paste0("chr",chr)) %>%
  write.table("177sa_vcfextract_july2025.list", sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE)

# Full SNP list
results %>%
  arrange(pval) %>%
  select(chr, pos) %>%
  mutate(chr=case_when(chr==23~"X", chr==24~"Y", TRUE~as.character(chr)),
         chr=paste0("chr",chr)) %>%
  write.table("177sa_vcfextract_all_july2025.list", sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
EOF

###############################################################################
# Step 8: Extract SNPs from VCF and annotate (SnpSift)
###############################################################################
vcftools --gzvcf snps_ann_all.vcf.gz \
    --positions 177sa_vcfextract_july2025.list \
    --recode --recode-INFO-all \
    --out 177sa_pval

vcftools --gzvcf snps_ann_all.vcf.gz \
    --positions 177sa_vcfextract_all_july2025.list \
    --recode --recode-INFO-all \
    --out 177sa_all

SnpSift extractFields -s "," -e "." 177sa_pval.recode.vcf CHROM POS REF ALT ANN[*] > INFO_177sa_pval.txt
SnpSift extractFields -s "," -e "." 177sa_all.recode.vcf CHROM POS REF ALT ANN[*] > INFO_177sa_all.txt

###############################################################################
# Step 9: Parse VCF annotations (R)
###############################################################################
Rscript - <<'EOF'
library(data.table); library(dplyr); library(tidyr); library(vcfR)

# Load and parse
vcf <- read.vcfR("177sa_all.recode.vcf")
fix <- as_tibble(vcf@fix) %>% mutate(POS=as.numeric(POS)) %>% select(-INFO)

info <- read.table("INFO_177sa_all.txt", sep="\t", header=TRUE) %>%
        as_tibble() %>%
        separate_rows(ANN..., sep=",") %>%
        separate(ANN..., into=c("EFFECT_ALLELE","EFFECT","IMPACT","GENE","GENEID",
                                "FEATURE","FEATUREID","BIOTYPE","RANK","HGVS_C","HGVS_P",
                                "CDNA_POS","CDNA_LEN","CDS_POS","CDS_LEN","ERRORS"), sep="\\|") %>%
        select(-FEATUREID, -c(RANK:ERRORS)) %>%
        distinct()

vcf_info <- full_join(fix, info)
write.csv(vcf_info, "177sa_all_vcfinfo.csv", row.names=FALSE)
EOF

###############################################################################
# Step 10: LD analysis of candidate regions
###############################################################################
plink --bfile july2025_sev3_nohet_qcfiltered_maf005 \
    --chr 3 \
    --out chr3_ordinal177_highLD

plink --bfile july2025_sev3_nohet_qcfiltered_maf005 \
    --r2 \
    --ld-snp "rs7649631;1269087" \
    --ld-window-kb 1000 \
    --ld-window 99999 \
    --ld-window-r2 0
