#!/bin/bash
# ----------------QSUB Parameters-----------------
#PBS -q q2
#PBS -M dbrkic@bioinfo.hr
#PBS -m n
#PBS -N plink
#PBS -l select=ncpus=15:mem=200g
#PBS -j oe

cd $PBS_O_WORKDIR


# ----------------Commands------------------- #

#turning vcf to plink files
plink --vcf snps_ann_all.vcf.gz \
    --recode \
    --keep keep_samples.txt \
    --out snps_ann

#add sex
plink --file snps_ann \
    --update-sex sex_file.list \
    --make-bed \
    --out snps_ann_sex

#fix sex
plink --bfile snps_ann_sex \
    --impute-sex 0.964 0.969 \
    --make-bed \
    --out snps_ann_sex_imp

#check sex
plink --bfile snps_ann_sex_imp \
    --check-sex \
    --out snps_ann_imp

#add phenotype
plink --bfile snps_ann_sex \
    --pheno pheno_file.list \
    --make-bed \
    --out snps_ann_sex_pheno

# quality filtering 
plink --bfile snps_ann_sex \
    --geno 0.05 \
    --hwe 0.001 \
    --keep pheno_july2025.list \
    --maf 0.05 \
    --make-bed \
    --mind 0.05 \
    --out july2025_sev3_nohet_qcfiltered_maf005 \
    --pheno pheno_july2025.list \
    --pheno-name severity \
    --remove low_het_individuals.txt

#Total genotyping rate in remaining samples is 0.999948.
#1036 variants removed due to missing genotype data (--geno).
#--hwe: 1642211 variants removed due to Hardy-Weinberg exact test.
#3653128 variants removed due to minor allele threshold(s)
#(--maf/--max-maf/--mac/--max-mac).
#167152 variants and 177 people pass filters and QC.
#Phenotype data is quantitative.

plink --bfile july2025_sev3_nohet_qcfiltered_maf005 \
    --indep 50 5 2 \
    --out july2025_sev3_nohet_qcfiltered_maf005_pruned
#Pruning complete.  116899 of 167152 variants removed.

plink --bfile july2025_sev3_nohet_qcfiltered_maf005 \
    --extract july2025_sev3_nohet_qcfiltered_maf005_pruned.prune.in \
    --make-bed \
    --out july2025_sev3_nohet_qcfiltered_maf005_pruneddata
#Total genotyping rate is 0.999808.
#50377 variants and 177 people pass filters and QC.
#Phenotype data is quantitative.

plink --bfile july2025_sev3_nohet_qcfiltered_maf005_pruneddata \
    --pca 20 \
    --out july2025_sev3_nohet_qcfiltered_maf005_pca
#Excluding 1884 variants on non-autosomes from relationship matrix calc.
#Relationship matrix calculation complete.

echo -e "FID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10\tPC11\tPC12\tPC13\tPC14\tPC15\tPC16\tPC17\tPC18\tPC19\tPC20"| cat - july2025_sev3_nohet_qcfiltered_maf005_pca.eigenvec > july2025_sev3_nohet_qcfiltered_maf005_pca_header.eigenvec

################################

#Creating .csv file of PCs
source activate r_env

R

library(ggplot2)
library(magrittr)
library(dplyr)

ret_pca<-read.table("july2025_sev3_nohet_qcfiltered_maf005_pca_header.eigenvec", header=T) %>% as_tibble()

pheno <- read.table("pheno_july2025.list", header=T) %>% as_tibble() 

covariates <- right_join(pheno, ret_pca %>% select(FID:PC10))

write.csv(covariates, "july2025_covariatesPC_177.csv", row.names=F, quote=F)

q()

################################

conda activate julia_env

julia

using OrdinalGWAS
using DataFrames

# Define file paths
const datadir = "ind177_input"
const covfile = datadir * "/july2025_covariatesPC_177.csv"
const plkfile = datadir * "/july2025_sev3_nohet_qcfiltered_maf005"

# Run ordinal GWAS with covariates
ordinalgwas(@formula(severity ~ sex + age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10), covfile, plkfile, test = :lrt, analysistype = "singlesnp")

exit()

################################

#Plotting ordinal GWAS results
source activate r_env

R

library(qqman)
library(data.table)
library(writexl)
library(dplyr)
library(stringr)

results<-as.data.frame(fread("ordinalgwas.pval.txt", header=T, sep=",", stringsAsFactors=FALSE)) %>% filter(hwepval>0.05) %>% filter(chr<23)

png("Manhattan_plot.png", width=20,height= 15,units= "cm",res= 600)
manhattan(results, chr="chr", bp="pos", snp="snpid", p="pval", main="Manhattan Plot of ordinal severity\n(177 individuals, sex, age, PC10, MAF>5%) ", col = c("blue","gray"), chrlabs = NULL,
suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval=0.0001, ylim = c(0, 8.5))
dev.off()

png("qqlot_plot.png")
qq(results$pval, main="QQ plot of ordinal severity MAF>5%")
dev.off()

results_sorted <- results[order(results$pval), ]
top30<-results_sorted[1:30,]

write_xlsx(top30, "Top30_177_ordinal_variants_july2025.xlsx")

results %>% 
filter(pval<0.05) %>% 
arrange(pval) %>%
dplyr::select(chr, pos) %>%
mutate(chr=as.character(chr),
       chr=str_replace(chr, "23", "X"),
       chr=str_replace(chr, "24", "Y"),
       chr=str_replace(chr, "^", "chr")) %>% 
write.table("177sa_vcfextract_july2025.list", quote = F, sep = "\t", row.names = F, col.names = F)

results %>% 
arrange(pval) %>%
dplyr::select(chr, pos) %>%
mutate(chr=as.character(chr),
       chr=str_replace(chr, "23", "X"),
       chr=str_replace(chr, "24", "Y"),
       chr=str_replace(chr, "^", "chr")) %>% 
write.table("177sa_vcfextract_all_july2025.list", quote = F, sep = "\t", row.names = F, col.names = F)

################################

#Extracting SNPs from snps_ann_all.vcf.gz based on position 
vcftools --gzvcf snps_ann_all.vcf.gz \
    --positions 177sa_vcfextract_july2025.list \
    --recode \
    --recode-INFO-all \
    --out 177sa_pval

vcftools --gzvcf snps_ann_all.vcf.gz \
    --positions 177sa_vcfextract_all_july2025.list \
    --recode \
    --recode-INFO-all \
    --out 177sa_all

#Extracting INFO field form snps_ann_all.vcf.gz in terminal
SnpSift extractFields -s "," -e "." 177sa_pval.recode.vcf CHROM POS REF ALT ANN[*] > INFO_177sa_pval.recode.vcf.txt
SnpSift extractFields -s "," -e "." 177sa_all.recode.vcf CHROM POS REF ALT ANN[*] > INFO_177sa_all.recode.vcf.txt


################################

source activate r_env

R

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)
library(vcfR)

#MAF >5%
vcf <- read.vcfR("177sa_pval.recode.vcf")
vcf <- read.vcfR("177sa_all.recode.vcf")

vcf1 <- vcf@fix %>% 
  as_tibble() %>% 
  mutate(POS=as.numeric(POS)) %>% 
  dplyr::select(-INFO)
  
info <-
  read.table("INFO_177sa_pval.recode.vcf.txt", sep = "\t", header = T) %>% 
  as_tibble() %>%
  separate_rows(ANN..., sep = "\\,") %>%
  separate(ANN..., c("EFFECT_ALLELE","EFFECT","IMPACT","GENE","GENEID","FEATURE","FEATUREID","BIOTYPE","RANK","HGVS_C","HGVS_P","CDNA_POS","CDNA_LEN","CDS_POS","CDS_LEN","ERRORS"), sep = "\\|") %>% 
  dplyr::select(-FEATUREID, -c(RANK:ERRORS)) %>% 
  distinct()

  info <-
  read.table("INFO_177sa_all.recode.vcf.txt", sep = "\t", header = T) %>% 
  as_tibble() %>%
  separate_rows(ANN..., sep = "\\,") %>%
  separate(ANN..., c("EFFECT_ALLELE","EFFECT","IMPACT","GENE","GENEID","FEATURE","FEATUREID","BIOTYPE","RANK","HGVS_C","HGVS_P","CDNA_POS","CDNA_LEN","CDS_POS","CDS_LEN","ERRORS"), sep = "\\|") %>% 
  dplyr::select(-FEATUREID, -c(RANK:ERRORS)) %>% 
  distinct()

vcf_info <- full_join(vcf1, info)

write.csv(vcf_info, "177sa_pval_vcfinfo.csv", row.names = F)
write.csv(vcf_info, "177sa_all_vcfinfo.csv", row.names = F)

################################

#High LD region

plink --bfile july2025_sev3_nohet_qcfiltered_maf005 \
    --out chr3_ordinal177_highLD \
    --chr 3 \
    --from-bp 

plink --bfile july2025_sev3_nohet_qcfiltered_maf005 \
    --r2 \
    --ld-snp "rs7649631;1269087" \
    --ld-window-kb 1000 \
    --ld-window 99999 \
    --ld-window-r2 0

