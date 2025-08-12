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
ts8 <- read.table(file = "E:/post_gwas_cancer/results/utmost,fusion,magma三者交集结果/乳腺癌/Breast_cancer_BCACIEU/overlap_gene49.txt",header = F)
probe <- Ensembl_GRCh37[Ensembl_GRCh37$gene_name %in% ts8$V1,"gene_id"]
write.table(probe,file = "E:/post_gwas_cancer/results/utmost,fusion,magma三者交集结果/乳腺癌/Breast_cancer_BCACIEU/overlap_gene49_probe.txt",quote = F)

# ② 使用qtl文件前缀，进行批量smr分析
for (i in qtl_besd) {
  print(i)
  smr_qtl2gwas(smr_exe_path = "F:/post-GWAS/data_smr/smr软件/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
               bfile = "F:/post-GWAS/data_smr/1kg.v3/EUR",
               gwas_path = "F:/post-GWAS/tmp8_smr_brain_axis/Liver_cancer_GWASCatalog.ma",
               qtls_path = i,
               probes_path = "E:/post_gwas_cancer/results/utmost,fusion,magma三者交集结果/乳腺癌/Breast_cancer_BCACIEU/overlap_gene49_probe.txt",
               smr_multi_snp = T , # 基于多个SNP计算p值。
               smr_peqtl = 5e-08,
               smr_cis_wind = 1000,  # cis_wind值得是探针上游或者下游的碱基数量，即单侧碱基数量，单位为kb
               
               out_path = "E:/post_gwas_cancer/results/Breast_cancer_BCACIEU_smr/",
               out_prefix = basename(i))
}

# 结果处理
dir <- list.files("E:/post_gwas_cancer/results/Breast_cancer_BCACIEU_smr/results/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/results/Breast_cancer_BCACIEU_smr/results/"
n <- length(dir)
for (i in 1:n) {
  file = paste0("E:/post_gwas_cancer/results/Breast_cancer_BCACIEU_smr/results/",dir[i])
  tissues <- sub("\\.lite.*", "",basename(dir[i]))
  file_name = paste0("E:/post_gwas_cancer/results/Breast_cancer_BCACIEU_smr/result_tissue/",tissues,"smr_results.csv")
  eqtl2gwas <- data.table::fread(file = file,
                                 data.table = F)
  eqtl2gwas$p_SMR_fdr <- p.adjust(eqtl2gwas$p_SMR,method = "fdr")
  eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_HEIDI > 0.01,]
  eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_SMR_fdr < 0.05,]
  eqtl2gwas <- na.omit(eqtl2gwas)
  if (nrow(eqtl2gwas) > 0) {
    eqtl2gwas$tissues <- tissues
    write.csv(eqtl2gwas, file = file_name, row.names = FALSE)
  } else {
    message(paste("No data for", tissues, "after filtering."))
  }
}

dir <- list.files("E:/post_gwas_cancer/results/Breast_cancer_BCACIEU_smr/result_tissue/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/results/Breast_cancer_BCACIEU_smr/result_tissue/"
n <- length(dir)
data1 <- read.csv(file = paste0(data_dir, dir[1]), sep = ",", header = TRUE)
for (i in 2:n) {
  input_path = paste0(data_dir, dir[i])
  data2 <- read.csv(file = input_path, sep = ",", header = TRUE)
  data1 <- rbind(data1, data2)
}
write.csv(data1, file = "E:/post_gwas_cancer/results/Breast_cancer_BCACIEU_smr/Breast_cancer_BCACIEU_smr汇总结果.csv", row.names = FALSE)


######### ER-乳腺癌 ###########
## 数据预处理
# smr_prepare_ma(file_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/Melanoma_FinnGenR10.ma",
#                out_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/")

## 如果是多个组织 smr_qtl2gwas分析：
# ① 获取qtl文件前缀
qtl_besd <- dir("F:/post-GWAS/data_smr/GTEx_V8_cis_eqtl_summary_lite/",pattern = ".besd$",full.names = T)
qtl_besd <- gsub(".besd","",qtl_besd)
ts8 <- read.table(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/ER-_Breast_cancer/overlap_gene3.txt",header = F)
probe <- Ensembl_GRCh37[Ensembl_GRCh37$gene_name %in% ts8$V1,"gene_id"]
write.table(probe,file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/ER-_Breast_cancer/overlap_gene3_probe.txt",quote = F)

# ② 使用qtl文件前缀，进行批量smr分析
for (i in qtl_besd) {
  print(i)
  smr_qtl2gwas(smr_exe_path = "F:/post-GWAS/data_smr/smr软件/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
               bfile = "F:/post-GWAS/data_smr/1kg.v3/EUR",
               gwas_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/ER-_Breast_cancer_BCACIEU.ma",
               qtls_path = i,
               probes_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/ER-_Breast_cancer/overlap_gene3_probe.txt",
               smr_multi_snp = T , # 基于多个SNP计算p值。
               smr_peqtl = 5e-08,
               smr_cis_wind = 1000,  # cis_wind值得是探针上游或者下游的碱基数量，即单侧碱基数量，单位为kb
               
               out_path = "E:/post_gwas_cancer/results/ER-_Breast_cancer_BCACIEU_smr/",
               out_prefix = basename(i))
}

# 结果处理
dir <- list.files("E:/post_gwas_cancer/results/ER-_Breast_cancer_BCACIEU_smr/results/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/results/ER-_Breast_cancer_BCACIEU_smr/results/"
n <- length(dir)
for (i in 1:n) {
  file = paste0("E:/post_gwas_cancer/results/ER-_Breast_cancer_BCACIEU_smr/results/",dir[i])
  tissues <- sub("\\.lite.*", "",basename(dir[i]))
  file_name = paste0("E:/post_gwas_cancer/results/ER-_Breast_cancer_BCACIEU_smr/result_tissue/",tissues,"smr_results.csv")
  eqtl2gwas <- data.table::fread(file = file,
                                 data.table = F)
  eqtl2gwas$p_SMR_fdr <- p.adjust(eqtl2gwas$p_SMR,method = "fdr")
  eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_HEIDI > 0.01,]
  eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_SMR_fdr < 0.05,]
  eqtl2gwas <- na.omit(eqtl2gwas)
  if (nrow(eqtl2gwas) > 0) {
    eqtl2gwas$tissues <- tissues
    write.csv(eqtl2gwas, file = file_name, row.names = FALSE)
  } else {
    message(paste("No data for", tissues, "after filtering."))
  }
}

dir <- list.files("E:/post_gwas_cancer/results/ER-_Breast_cancer_BCACIEU_smr/result_tissue/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/results/ER-_Breast_cancer_BCACIEU_smr/result_tissue/"
n <- length(dir)
data1 <- read.csv(file = paste0(data_dir, dir[1]), sep = ",", header = TRUE)
for (i in 2:n) {
  input_path = paste0(data_dir, dir[i])
  data2 <- read.csv(file = input_path, sep = ",", header = TRUE)
  data1 <- rbind(data1, data2)
}
write.csv(data1, file = "E:/post_gwas_cancer/results/ER-_Breast_cancer_BCACIEU_smr/ER-_Breast_cancer_BCACIEU_smr汇总结果.csv", row.names = FALSE)

######### ER+乳腺癌 ###########
## 数据预处理
# smr_prepare_ma(file_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/Melanoma_FinnGenR10.ma",
#                out_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/")

## 如果是多个组织 smr_qtl2gwas分析：
# ① 获取qtl文件前缀
qtl_besd <- dir("F:/post-GWAS/data_smr/GTEx_V8_cis_eqtl_summary_lite/",pattern = ".besd$",full.names = T)
qtl_besd <- gsub(".besd","",qtl_besd)
ts8 <- read.table(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/ER+_Breast_cancer/overlap_gene24.txt")
probe <- Ensembl_GRCh37[Ensembl_GRCh37$gene_name %in% ts8$V1,"gene_id"]
write.table(probe,file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/ER+_Breast_cancer/overlap_gene24_probe.txt",quote = F)

# ② 使用qtl文件前缀，进行批量smr分析
for (i in qtl_besd) {
  print(i)
  smr_qtl2gwas(smr_exe_path = "F:/post-GWAS/data_smr/smr软件/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
               bfile = "F:/post-GWAS/data_smr/1kg.v3/EUR",
               gwas_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/ER+_Breast_cancer_BCACIEU.ma",
               qtls_path = i,
               probes_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/ER+_Breast_cancer/overlap_gene24_probe.txt",
               smr_multi_snp = T , # 基于多个SNP计算p值。
               smr_peqtl = 5e-08,
               smr_cis_wind = 1000,  # cis_wind值得是探针上游或者下游的碱基数量，即单侧碱基数量，单位为kb
               
               out_path = "E:/post_gwas_cancer/results/ER+_Breast_cancer_BCACIEU_smr/",
               out_prefix = basename(i))
}

# 结果处理
dir <- list.files("E:/post_gwas_cancer/results/ER+_Breast_cancer_BCACIEU_smr/results/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/results/ER+_Breast_cancer_BCACIEU_smr/results/"
n <- length(dir)
for (i in 1:n) {
  file = paste0("E:/post_gwas_cancer/results/ER+_Breast_cancer_BCACIEU_smr/results/",dir[i])
  tissues <- sub("\\.lite.*", "",basename(dir[i]))
  file_name = paste0("E:/post_gwas_cancer/results/ER+_Breast_cancer_BCACIEU_smr/result_tissue/",tissues,"smr_results.csv")
  eqtl2gwas <- data.table::fread(file = file,
                                 data.table = F)
  eqtl2gwas$p_SMR_fdr <- p.adjust(eqtl2gwas$p_SMR,method = "fdr")
  eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_HEIDI > 0.01,]
  eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_SMR_fdr < 0.05,]
  eqtl2gwas <- na.omit(eqtl2gwas)
  if (nrow(eqtl2gwas) > 0) {
    eqtl2gwas$tissues <- tissues
    write.csv(eqtl2gwas, file = file_name, row.names = FALSE)
  } else {
    message(paste("No data for", tissues, "after filtering."))
  }
}

dir <- list.files("E:/post_gwas_cancer/results/ER+_Breast_cancer_BCACIEU_smr/result_tissue/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/results/ER+_Breast_cancer_BCACIEU_smr/result_tissue/"
n <- length(dir)
data1 <- read.csv(file = paste0(data_dir, dir[1]), sep = ",", header = TRUE)
for (i in 2:n) {
  input_path = paste0(data_dir, dir[i])
  data2 <- read.csv(file = input_path, sep = ",", header = TRUE)
  data1 <- rbind(data1, data2)
}
write.csv(data1, file = "E:/post_gwas_cancer/results/ER+_Breast_cancer_BCACIEU_smr/ER+_Breast_cancer_BCACIEU_smr汇总结果.csv", row.names = FALSE)


######### 黑色素瘤 ###########
## 数据预处理
# smr_prepare_ma(file_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/Melanoma_FinnGenR10.ma",
#                out_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/")

## 如果是多个组织 smr_qtl2gwas分析：
# ① 获取qtl文件前缀
qtl_besd <- dir("F:/post-GWAS/data_smr/GTEx_V8_cis_eqtl_summary_lite/",pattern = ".besd$",full.names = T)
qtl_besd <- gsub(".besd","",qtl_besd)
ts8 <- read.table(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/黑色素瘤FinnGen/overlap_gene3.txt")
probe <- Ensembl_GRCh37[Ensembl_GRCh37$gene_name %in% ts8$V1,"gene_id"]
write.table(probe,file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/黑色素瘤FinnGen/overlap_gene3_probe.txt",quote = F)

# ② 使用qtl文件前缀，进行批量smr分析
for (i in qtl_besd) {
  print(i)
  smr_qtl2gwas(smr_exe_path = "F:/post-GWAS/data_smr/smr软件/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
               bfile = "F:/post-GWAS/data_smr/1kg.v3/EUR",
               gwas_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/Melanoma_FinnGenR10_treated.ma",
               qtls_path = i,
               probes_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/黑色素瘤FinnGen/overlap_gene3_probe.txt",
               smr_multi_snp = T , # 基于多个SNP计算p值。
               smr_peqtl = 5e-08,
               smr_cis_wind = 1000,  # cis_wind值得是探针上游或者下游的碱基数量，即单侧碱基数量，单位为kb
               
               out_path = "E:/post_gwas_cancer/results/Melanoma_FinnGenR10_smr/",
               out_prefix = basename(i))
}

# 结果处理
dir <- list.files("E:/post_gwas_cancer/results/Melanoma_FinnGenR10_smr/results/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/results/Melanoma_FinnGenR10_smr/results/"
n <- length(dir)
for (i in 1:n) {
  file = paste0("E:/post_gwas_cancer/results/Melanoma_FinnGenR10_smr/results/",dir[i])
  tissues <- sub("\\.lite.*", "",basename(dir[i]))
  file_name = paste0("E:/post_gwas_cancer/results/Melanoma_FinnGenR10_smr/result_tissue/",tissues,"smr_results.csv")
  eqtl2gwas <- data.table::fread(file = file,
                                 data.table = F)
  eqtl2gwas$p_SMR_fdr <- p.adjust(eqtl2gwas$p_SMR,method = "fdr")
  eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_HEIDI > 0.01,]
  eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_SMR_fdr < 0.05,]
  eqtl2gwas <- na.omit(eqtl2gwas)
  if (nrow(eqtl2gwas) > 0) {
    eqtl2gwas$tissues <- tissues
    write.csv(eqtl2gwas, file = file_name, row.names = FALSE)
  } else {
    message(paste("No data for", tissues, "after filtering."))
  }
}

dir <- list.files("E:/post_gwas_cancer/results/Melanoma_FinnGenR10_smr/result_tissue/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/results/Melanoma_FinnGenR10_smr/result_tissue/"
n <- length(dir)
data1 <- read.csv(file = paste0(data_dir, dir[1]), sep = ",", header = TRUE)
for (i in 2:n) {
  input_path = paste0(data_dir, dir[i])
  data2 <- read.csv(file = input_path, sep = ",", header = TRUE)
  data1 <- rbind(data1, data2)
}
write.csv(data1, file = "E:/post_gwas_cancer/results/Melanoma_FinnGenR10_smr/Melanoma_FinnGenR10_smr汇总结果.csv", row.names = FALSE)


######### 前列腺癌 ###########
## 数据预处理
smr_prepare_ma(file_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/Prostate_cancer_FinnGenR10.ma",
               out_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/")

## 如果是多个组织 smr_qtl2gwas分析：
# ① 获取qtl文件前缀
qtl_besd <- dir("F:/post-GWAS/data_smr/GTEx_V8_cis_eqtl_summary_lite/",pattern = ".besd$",full.names = T)
qtl_besd <- gsub(".besd","",qtl_besd)
ts8 <- read.table(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/前列腺癌/overlap_gene16.txt")
probe <- Ensembl_GRCh37[Ensembl_GRCh37$gene_name %in% ts8$V1,"gene_id"]
write.table(probe,file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/前列腺癌/overlap_gene16_probe.txt",quote = F)

# ② 使用qtl文件前缀，进行批量smr分析
for (i in qtl_besd) {
  print(i)
  smr_qtl2gwas(smr_exe_path = "F:/post-GWAS/data_smr/smr软件/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
               bfile = "F:/post-GWAS/data_smr/1kg.v3/EUR",
               gwas_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/Prostate_cancer_FinnGenR10_treated.ma",
               qtls_path = i,
               probes_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/前列腺癌/overlap_gene16_probe.txt",
               smr_multi_snp = T , # 基于多个SNP计算p值。
               smr_peqtl = 5e-08,
               smr_cis_wind = 1000,  # cis_wind值得是探针上游或者下游的碱基数量，即单侧碱基数量，单位为kb
               
               out_path = "E:/post_gwas_cancer/results/Prostate_cancer_smr/",
               out_prefix = basename(i))
}

# 结果处理
dir <- list.files("E:/post_gwas_cancer/results/Prostate_cancer_smr/results/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/results/Prostate_cancer_smr/results/"
n <- length(dir)
for (i in 1:n) {
  file = paste0("E:/post_gwas_cancer/results/Prostate_cancer_smr/results/",dir[i])
  tissues <- sub("\\.lite.*", "",basename(dir[i]))
  file_name = paste0("E:/post_gwas_cancer/results/Prostate_cancer_smr/result_tissue/",tissues,"smr_results.csv")
  eqtl2gwas <- data.table::fread(file = file,
                                 data.table = F)
  eqtl2gwas$p_SMR_fdr <- p.adjust(eqtl2gwas$p_SMR,method = "fdr")
  eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_HEIDI > 0.01,]
  eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_SMR_fdr < 0.05,]
  eqtl2gwas <- na.omit(eqtl2gwas)
  if (nrow(eqtl2gwas) > 0) {
    eqtl2gwas$tissues <- tissues
    write.csv(eqtl2gwas, file = file_name, row.names = FALSE)
  } else {
    message(paste("No data for", tissues, "after filtering."))
  }
}

dir <- list.files("E:/post_gwas_cancer/results/Prostate_cancer_smr/result_tissue/")# GWAS所在文件夹
data_dir <- "E:/post_gwas_cancer/results/Prostate_cancer_smr/result_tissue/"
n <- length(dir)
data1 <- read.csv(file = paste0(data_dir, dir[1]), sep = ",", header = TRUE)
for (i in 2:n) {
  input_path = paste0(data_dir, dir[i])
  data2 <- read.csv(file = input_path, sep = ",", header = TRUE)
  data1 <- rbind(data1, data2)
}
write.csv(data1, file = "E:/post_gwas_cancer/results/Prostate_cancer_smr/Prostate_cancer_FinnGenR10_smr汇总结果.csv", row.names = FALSE)

