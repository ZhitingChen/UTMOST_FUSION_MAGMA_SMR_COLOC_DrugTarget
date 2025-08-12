### FUSION后SMR分析  ###
library(DrugTargetMR)
library(openxlsx)
library(data.table)
rm(list = ls())

######### 乳腺癌 ###########
## 数据预处理
# smr_prepare_ma(file_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/Melanoma_FinnGenR10.ma",
#                out_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/")

## 如果是多个组织 smr_qtl2gwas分析：
# ① 获取qtl文件前缀
qtl_besd <- dir("F:/post-GWAS/data_smr/GTEx_V8_cis_eqtl_summary_lite/",pattern = ".besd$",full.names = T)
qtl_besd <- gsub(".besd","",qtl_besd)

# ② 使用qtl文件前缀，进行批量smr分析
for (i in qtl_besd) {
  print(i)
  smr_qtl2gwas(smr_exe_path = "F:/post-GWAS/data_smr/smr软件/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
               bfile = "F:/post-GWAS/data_smr/1kg.v3/EUR",
               gwas_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/finngen_R10_C3_BREAST_EXALLC_treated.ma",
               qtls_path = i,
               probes_path = "E:/post_gwas_cancer/tmp1/BC_probe.txt",
               smr_multi_snp = T , # 基于多个SNP计算p值。
               smr_peqtl = 5e-08,
               smr_cis_wind = 1000,  # cis_wind值得是探针上游或者下游的碱基数量，即单侧碱基数量，单位为kb
               
               out_path = "E:/post_gwas_cancer/tmp/Breast_cancer_FinnGen_smr/result/",
               out_prefix = basename(i))
}

# 结果处理
dir <- list.files("E:/post_gwas_cancer/tmp/Breast_cancer_FinnGen_smr/result/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/tmp/Breast_cancer_FinnGen_smr/result/"
n <- length(dir)
for (i in 1:n) {
  file = paste0("E:/post_gwas_cancer/tmp/Breast_cancer_FinnGen_smr/result/",dir[i])
  tissues <- sub("\\.lite.*", "",basename(dir[i]))
  file_name = paste0("E:/post_gwas_cancer/tmp/Breast_cancer_FinnGen_smr/result_tissue/",tissues,"smr_results.csv")
  eqtl2gwas <- data.table::fread(file = file,
                                 data.table = F)
  eqtl2gwas$p_SMR_fdr <- p.adjust(eqtl2gwas$p_SMR,method = "fdr")
  # eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_HEIDI > 0.01,]
  # eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_SMR_fdr < 0.05,]
  eqtl2gwas <- na.omit(eqtl2gwas)
  if (nrow(eqtl2gwas) > 0) {
    eqtl2gwas$tissues <- tissues
    write.csv(eqtl2gwas, file = file_name, row.names = FALSE)
  } else {
    message(paste("No data for", tissues, "after filtering."))
  }
}

dir <- list.files("E:/post_gwas_cancer/tmp/Breast_cancer_FinnGen_smr/result_tissue/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/tmp/Breast_cancer_FinnGen_smr/result_tissue/"
n <- length(dir)
data1 <- read.csv(file = paste0(data_dir, dir[1]), sep = ",", header = TRUE)
for (i in 2:n) {
  input_path = paste0(data_dir, dir[i])
  data2 <- read.csv(file = input_path, sep = ",", header = TRUE)
  data1 <- rbind(data1, data2)
}
write.csv(data1, file = "E:/post_gwas_cancer/tmp/Breast_cancer_FinnGen_smr/Breast_cancer_FinnGen_smr汇总结果.csv", row.names = FALSE)

######### ER-乳腺癌 ###########
## 数据预处理
# smr_prepare_ma(file_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/Melanoma_FinnGenR10.ma",
#                out_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/")

## 如果是多个组织 smr_qtl2gwas分析：
# ① 获取qtl文件前缀
qtl_besd <- dir("F:/post-GWAS/data_smr/GTEx_V8_cis_eqtl_summary_lite/",pattern = ".besd$",full.names = T)
qtl_besd <- gsub(".besd","",qtl_besd)

# ② 使用qtl文件前缀，进行批量smr分析
for (i in qtl_besd) {
  print(i)
  smr_qtl2gwas(smr_exe_path = "F:/post-GWAS/data_smr/smr软件/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
               bfile = "F:/post-GWAS/data_smr/1kg.v3/EUR",
               gwas_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/finngen_R10_C3_BREAST_ERNEG_EXALLC_treated.ma",
               qtls_path = i,
               probes_path = "E:/post_gwas_cancer/tmp1/ER-BC_probe.txt",
               smr_multi_snp = T , # 基于多个SNP计算p值。
               smr_peqtl = 5e-08,
               smr_cis_wind = 1000,  # cis_wind值得是探针上游或者下游的碱基数量，即单侧碱基数量，单位为kb
               
               out_path = "E:/post_gwas_cancer/tmp/ER-_Breast_cancer_FinnGen_smr/result/",
               out_prefix = basename(i))
}

# 结果处理
dir <- list.files("E:/post_gwas_cancer/tmp/ER-_Breast_cancer_FinnGen_smr/result/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/tmp/ER-_Breast_cancer_FinnGen_smr/result/"
n <- length(dir)
for (i in 1:n) {
  file = paste0("E:/post_gwas_cancer/tmp/ER-_Breast_cancer_FinnGen_smr/result/",dir[i])
  tissues <- sub("\\.lite.*", "",basename(dir[i]))
  file_name = paste0("E:/post_gwas_cancer/tmp/ER-_Breast_cancer_FinnGen_smr/result_tissue/",tissues,"smr_results.csv")
  eqtl2gwas <- data.table::fread(file = file,
                                 data.table = F)
  eqtl2gwas$p_SMR_fdr <- p.adjust(eqtl2gwas$p_SMR,method = "fdr")
  # eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_HEIDI > 0.01,]
  # eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_SMR_fdr < 0.05,]
  eqtl2gwas <- na.omit(eqtl2gwas)
  if (nrow(eqtl2gwas) > 0) {
    eqtl2gwas$tissues <- tissues
    write.csv(eqtl2gwas, file = file_name, row.names = FALSE)
  } else {
    message(paste("No data for", tissues, "after filtering."))
  }
}

dir <- list.files("E:/post_gwas_cancer/tmp/ER-_Breast_cancer_FinnGen_smr/result_tissue/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/tmp/ER-_Breast_cancer_FinnGen_smr/result_tissue/"
n <- length(dir)
data1 <- read.csv(file = paste0(data_dir, dir[1]), sep = ",", header = TRUE)
for (i in 2:n) {
  input_path = paste0(data_dir, dir[i])
  data2 <- read.csv(file = input_path, sep = ",", header = TRUE)
  data1 <- rbind(data1, data2)
}
write.csv(data1, file = "E:/post_gwas_cancer/tmp/ER-_Breast_cancer_FinnGen_smr/ER-_Breast_cancer_FinnGen_smr汇总结果.csv", row.names = FALSE)

######### ER+乳腺癌 ###########
## 数据预处理
# smr_prepare_ma(file_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/Melanoma_FinnGenR10.ma",
#                out_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/")

## 如果是多个组织 smr_qtl2gwas分析：
# ① 获取qtl文件前缀
qtl_besd <- dir("F:/post-GWAS/data_smr/GTEx_V8_cis_eqtl_summary_lite/",pattern = ".besd$",full.names = T)
qtl_besd <- gsub(".besd","",qtl_besd)

# ② 使用qtl文件前缀，进行批量smr分析
for (i in qtl_besd) {
  print(i)
  smr_qtl2gwas(smr_exe_path = "F:/post-GWAS/data_smr/smr软件/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
               bfile = "F:/post-GWAS/data_smr/1kg.v3/EUR",
               gwas_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/finngen_R10_C3_BREAST_ERPLUS_EXALLC_treated.ma",
               qtls_path = i,
               probes_path = "E:/post_gwas_cancer/tmp1/ER+BC_probe.txt",
               smr_multi_snp = T , # 基于多个SNP计算p值。
               smr_peqtl = 5e-08,
               smr_cis_wind = 1000,  # cis_wind值得是探针上游或者下游的碱基数量，即单侧碱基数量，单位为kb
               
               out_path = "E:/post_gwas_cancer/tmp/ER+_Breast_cancer_FinnGen_smr/result/",
               out_prefix = basename(i))
}

# 结果处理
dir <- list.files("E:/post_gwas_cancer/tmp/ER+_Breast_cancer_FinnGen_smr/result/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/tmp/ER+_Breast_cancer_FinnGen_smr/result/"
n <- length(dir)
for (i in 1:n) {
  file = paste0("E:/post_gwas_cancer/tmp/ER+_Breast_cancer_FinnGen_smr/result/",dir[i])
  tissues <- sub("\\.lite.*", "",basename(dir[i]))
  file_name = paste0("E:/post_gwas_cancer/tmp/ER+_Breast_cancer_FinnGen_smr/result_tissue/",tissues,"smr_results.csv")
  eqtl2gwas <- data.table::fread(file = file,
                                 data.table = F)
  eqtl2gwas$p_SMR_fdr <- p.adjust(eqtl2gwas$p_SMR,method = "fdr")
  # eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_HEIDI > 0.01,]
  # eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_SMR_fdr < 0.05,]
  eqtl2gwas <- na.omit(eqtl2gwas)
  if (nrow(eqtl2gwas) > 0) {
    eqtl2gwas$tissues <- tissues
    write.csv(eqtl2gwas, file = file_name, row.names = FALSE)
  } else {
    message(paste("No data for", tissues, "after filtering."))
  }
}

dir <- list.files("E:/post_gwas_cancer/tmp/ER+_Breast_cancer_FinnGen_smr/result_tissue/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/tmp/ER+_Breast_cancer_FinnGen_smr/result_tissue/"
n <- length(dir)
data1 <- read.csv(file = paste0(data_dir, dir[1]), sep = ",", header = TRUE)
for (i in 2:n) {
  input_path = paste0(data_dir, dir[i])
  data2 <- read.csv(file = input_path, sep = ",", header = TRUE)
  data1 <- rbind(data1, data2)
}
write.csv(data1, file = "E:/post_gwas_cancer/tmp/ER+_Breast_cancer_FinnGen_smr/ER+_Breast_cancer_FinnGen_smr汇总结果.csv", row.names = FALSE)
