library(data.table)
library(tidyverse)
library(DrugTargetMR)

###  数据预处理
### method1
dir <- list.files("F:/post-GWAS/data_gwas_list/data_raw/method1/")# GWAS所在文件夹
data_dir <- "F:/post-GWAS/data_gwas_list/data_raw/method1/"
out_data_dir <- "F:/post-GWAS/data_gwas_list/data_raw/"
n <- length(dir)

for (i in 1:n) {
  file_path <- paste0(data_dir,dir[i])
  a <- data.table:::fread(file_path,header = T)
  head(a)
  
  a_use <- a %>%
    select(SNP,effect_allele,other_allele,effect_allele_frequency,beta,se,pval,samplesize,chr,pos)
  head(a_use)
  
  prfx <- sub("\\.csv$", "", dir[i])
  out_file_path <- paste0(out_data_dir,prfx,".txt") 
  write.table(a_use, out_file_path, sep = "\t", row.names = FALSE, quote = FALSE)

}


dir <- list.files("F:/post-GWAS/data_gwas_list/data_raw/new_method1/")# GWAS所在文件夹
data_dir <- "F:/post-GWAS/data_gwas_list/data_raw/new_method1"
n <- length(dir)

for (i in 1:n) {
  file_path <- file.path(data_dir, dir[i])  # 构造文件路径
  prepare_others(
    file_path = file_path,
    out_path = "F:/post-GWAS/data_gwas_list/data_prepare/",
    generate_mr = T,
    generate_smr = T,
    sep = "\t",
    reg_sep = F,
    col_SNP = "SNP",
    col_effect_allele = "effect_allele",
    col_other_allele = "other_allele",
    col_eaf = "effect_allele_frequency",
    col_beta = "beta",
    col_se = "se",
    col_pval = "pval",
    col_samplesize = "samplesize",
    col_chr = "chr",
    col_pos = "pos"
  )
  # 输出处理状态
  cat("Processed:", i, "\n")
}

### method2

dir <- list.files("F:/post-GWAS/data_gwas_list/data_raw/method2/")# GWAS所在文件夹
data_dir <- "F:/post-GWAS/data_gwas_list/data_raw/method2/"
out_data_dir <- "F:/post-GWAS/data_gwas_list/data_raw/"
n <- length(dir)

for (i in 1:n) {
  file_path <- paste0(data_dir,dir[i])
  a <- data.table:::fread(file_path,header = T)
  head(a)
  
  a_use <- a %>%
    select(variant_id,effect_allele,other_allele,effect_allele_frequency,beta,standard_error,p_value,samplesize,chromosome,base_pair_location)
  head(a_use)
  
  prfx <- sub("\\.csv$", "", dir[i])
  out_file_path <- paste0(out_data_dir,prfx,".txt") 
  write.table(a_use, out_file_path, sep = "\t", row.names = FALSE, quote = FALSE)
  
}


dir <- list.files("F:/post-GWAS/data_gwas_list/data_raw/new_method2/")# GWAS所在文件夹
data_dir <- "F:/post-GWAS/data_gwas_list/data_raw/new_method2"
n <- length(dir)

for (i in 1:n) {
  file_path <- file.path(data_dir, dir[i])  # 构造文件路径
  prepare_others(
    file_path = file_path,
    out_path = "F:/post-GWAS/data_gwas_list/data_prepare/",
    generate_mr = T,
    generate_smr = T,
    sep = "\t",
    reg_sep = F,
    col_SNP = "variant_id",
    col_effect_allele = "effect_allele",
    col_other_allele = "other_allele",
    col_eaf = "effect_allele_frequency",
    col_beta = "beta",
    col_se = "standard_error",
    col_pval = "p_value",
    col_samplesize = "samplesize",
    col_chr = "chromosome",
    col_pos = "base_pair_location"
  )
  # 输出处理状态
  cat("Processed:", i, "\n")
}


### method3
a <- data.table:::fread("F:/post-GWAS/data_gwas_list/data_raw/method3/Renal_cell_cancer_FinnGenR10.csv",header = T)
head(a)

a_use <- a %>%
  select(variant_id,effect_allele,other_allele,effect_allele_frequency_cases,beta,standard_error,p_value,samplesize,chromosome,base_pair_location)
head(a_use)

prfx <- sub("\\.csv$", "", "Renal_cell_cancer_FinnGenR10.csv")
out_file_path <- paste0("F:/post-GWAS/data_gwas_list/data_raw/",prfx,".txt") 
write.table(a_use, out_file_path, sep = "\t", row.names = FALSE, quote = FALSE)


prepare_others(
  file_path = "F:/post-GWAS/data_gwas_list/data_raw/new_method3/Renal_cell_cancer_FinnGenR10.txt",
  out_path = "F:/post-GWAS/data_gwas_list/data_prepare/",
  generate_mr = T,
  generate_smr = T,
  sep = "\t",
  reg_sep = F,
  col_SNP = "variant_id",
  col_effect_allele = "effect_allele",
  col_other_allele = "other_allele",
  col_eaf = "effect_allele_frequency_cases",
  col_beta = "beta",
  col_se = "standard_error",
  col_pval = "p_value",
  col_samplesize = "samplesize",
  col_chr = "chromosome",
  col_pos = "base_pair_location"
)

### 添加pheno用于MR分析 ###
dir <- list.files("F:/post-GWAS/data_gwas_list/data_prepare/mr/")# GWAS所在文件夹
data_dir <- "F:/post-GWAS/data_gwas_list/data_prepare/mr/"
out_data_dir <- "F:/post-GWAS/data_gwas_list/data_prepare/mr_with_phenotype/"
n <- length(dir)

for (i in 1:n) {
  file_path <- paste0(data_dir,dir[i])
  a <- data.table:::fread(file_path,header = T)
  head(a)
  
  pheno <- sub("\\.txt$", "", dir[i])
  
  a_use <- a %>%
    mutate(phenotype_col = pheno)
  head(a_use)
  
  prfx <- sub("\\.txt$", "", dir[i])
  out_file_path <- paste0(out_data_dir,prfx,".txt") 
  write.table(a_use, out_file_path, sep = "\t", row.names = FALSE, quote = FALSE)
  
}


### coloc 数据准备 ####
dir <- list.files("F:/post-GWAS/data_gwas_list/data_prepare/coloc/hg19/")# GWAS所在文件夹
data_dir <- "F:/post-GWAS/data_gwas_list/data_prepare/coloc/hg19/"
n <- length(dir)

for (i in 1:n) {
  file_path <- paste0(data_dir,dir[i])
  transform_hg19ToHg38(file_path = file_path, 
                        out_path = "F:/post-GWAS/data_gwas_list/data_prepare/coloc/")

}

### 添加pheno用于MR分析 ###

out_data_dir <- "F:/post-GWAS/data_gwas_list/data_prepare/mr_with_phenotype/"



file_path <- "F:/post-GWAS/data_gwas_list/data_prepare/coloc/finngen_R10_C3_BREAST_EXALLC.txt"
a <- data.table:::fread(file_path,header = T)
head(a)

pheno <- "Breast_Cancer_FinnGen"

a_use <- a %>%
  mutate(phenotype_col = pheno)
head(a_use)

prfx <- "Breast_Cancer_FinnGen"
out_file_path <- paste0(out_data_dir,prfx,".txt") 
write.table(a_use, out_file_path, sep = "\t", row.names = FALSE, quote = FALSE)
