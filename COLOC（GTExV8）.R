### FUSION后Coloc分析  ###
library(DrugTargetMR)
library(openxlsx)
library(data.table)
rm(list = ls())

################ Breast_cancer_BCACIEU ####################
######## 输入文件和样本量信息准备
## 对应的gtex文件，增加样本量信息
gtex_ss <- read.xlsx("F:/post-GWAS/data_smr/GTEx_v8_science_2020 附图和附表-附表2.xlsx")
gtex_ss$Samples

gtex_folder <- dir("F:/post-GWAS/data_smr/GTEx_V8_cis_eqtl_summary_hg19_48G/",full.names = T)
gtex_list <- data.frame(folder = gtex_folder,
                        Tissue = basename(gtex_folder))

gtex_list <- merge(gtex_list,gtex_ss[,c("Tissue","Samples")],by = "Tissue")


######## 共定位分析
### 准备用于查询的探针：对应.epi文件的v2列。
ts8 <- read.table(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/Breast_cancer_BCACIEU/overlap_gene49.txt")
probe <- Ensembl_GRCh37[Ensembl_GRCh37$gene_name %in% ts8$V1,"gene_id"]

### 使用for循环，完成批量共定位
for (i in 1:nrow(gtex_list)) {
  print(i)
  coloc_batch_smr(probes = probe,
                  smr_exe_path = "E:/药物靶点肺癌/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
                  qtl_path = gtex_list$folder[i], # 这里传对应文件夹的路径即可
                  cis_wind = 1000,
                  
                  # prepare_method=2,按照SNP的rsid选择SNP，以GWAS1的染色体坐标为参考，不使用GWAS2的chr和pos信息。
                  prepare_method = 2,
                  type1 = "quant",
                  SS1 = gtex_list$Samples[i],
                  gwas_path2 = "F:/post-GWAS/data_gwas_list/data_prepare/coloc/Breast_cancer_BCACIEU_hg38.txt",
                  type2 = "cc",
                  SS2 = 228951, # 对应sample.size
                  NC2 = 122977,        # 对应ncase
                  bfile = "E:/药物靶点肺癌/配套/EUR",
                  coloc_plot = T,
                  coloc_plot_pph4 = 0.6,
                  coloc_plot_genome = "hg37",
                  out_path = paste0("E:/post_gwas_cancer/results/tmp_coloc/",gtex_list$Tissue[i]))
}


######## 结果整理
# 获取组织对应的文件夹
sum_folder <- dir("E:/post_gwas_cancer/results/tmp_coloc/",full.names = T,include.dirs = T)

# 汇总组织和对应的结果文件
sum_dat <- data.frame()
for (f in sum_folder) {
  
  # 匹配文件夹和组织信息
  f_fn <- dir(f,pattern = "coloc_summary_coloc_abf.csv",recursive = T,full.names = T)
  if(length(f_fn)==0){
    next
  }
  tissues <- basename(f)
  
  # 获取coloc结果，并增加组织信息
  temp_dat <- read.csv(f_fn)
  temp_dat$tissues <- tissues
  sum_dat <- rbind(sum_dat,temp_dat)
}

# 将汇总后的结果写出
write.csv(sum_dat,file = "E:/post_gwas_cancer/results/tmp_coloc/Breast_cancer_BCACIEU汇总的共定位结果.csv",row.names = F)

################ ER-_Breast_cancer_BCACIEU ####################
######## 输入文件和样本量信息准备
## 对应的gtex文件，增加样本量信息
gtex_ss <- read.xlsx("F:/post-GWAS/data_smr/GTEx_v8_science_2020 附图和附表-附表2.xlsx")
gtex_ss$Samples

gtex_folder <- dir("F:/post-GWAS/data_smr/GTEx_V8_cis_eqtl_summary_hg19_48G/",full.names = T)
gtex_list <- data.frame(folder = gtex_folder,
                        Tissue = basename(gtex_folder))

gtex_list <- merge(gtex_list,gtex_ss[,c("Tissue","Samples")],by = "Tissue")


######## 共定位分析
### 准备用于查询的探针：对应.epi文件的v2列。
ts8 <- read.table(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/ER-_Breast_cancer/overlap_gene3.txt")
probe <- Ensembl_GRCh37[Ensembl_GRCh37$gene_name %in% ts8$V1,"gene_id"]

### 使用for循环，完成批量共定位
for (i in 1:nrow(gtex_list)) {
  print(i)
  coloc_batch_smr(probes = probe,
                  smr_exe_path = "E:/药物靶点肺癌/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
                  qtl_path = gtex_list$folder[i], # 这里传对应文件夹的路径即可
                  cis_wind = 1000,
                  
                  # prepare_method=2,按照SNP的rsid选择SNP，以GWAS1的染色体坐标为参考，不使用GWAS2的chr和pos信息。
                  prepare_method = 2,
                  type1 = "quant",
                  SS1 = gtex_list$Samples[i],
                  gwas_path2 = "F:/post-GWAS/data_gwas_list/data_prepare/coloc/ER-_Breast_cancer_BCACIEU_hg38.txt",
                  type2 = "cc",
                  SS2 = 127442, # 对应sample.size
                  NC2 = 21468,        # 对应ncase
                  bfile = "E:/药物靶点肺癌/配套/EUR",
                  coloc_plot = T,
                  coloc_plot_pph4 = 0.6,
                  coloc_plot_genome = "hg37",
                  out_path = paste0("E:/post_gwas_cancer/results/ER-_Breast_cancer_BCACIEU_coloc/",gtex_list$Tissue[i]))
}


######## 结果整理
# 获取组织对应的文件夹
sum_folder <- dir("E:/post_gwas_cancer/results/ER-_Breast_cancer_BCACIEU_coloc/",full.names = T,include.dirs = T)

# 汇总组织和对应的结果文件
sum_dat <- data.frame()
for (f in sum_folder) {
  
  # 匹配文件夹和组织信息
  f_fn <- dir(f,pattern = "coloc_summary_coloc_abf.csv",recursive = T,full.names = T)
  if(length(f_fn)==0){
    next
  }
  tissues <- basename(f)
  
  # 获取coloc结果，并增加组织信息
  temp_dat <- read.csv(f_fn)
  temp_dat$tissues <- tissues
  sum_dat <- rbind(sum_dat,temp_dat)
}

# 将汇总后的结果写出
write.csv(sum_dat,file = "E:/post_gwas_cancer/results/ER-_Breast_cancer_BCACIEU_coloc/ER-_Breast_cancer_BCACIEU汇总的共定位结果.csv",row.names = F)

################ ER+_Breast_cancer_BCACIEU ####################
######## 输入文件和样本量信息准备
## 对应的gtex文件，增加样本量信息
gtex_ss <- read.xlsx("F:/post-GWAS/data_smr/GTEx_v8_science_2020 附图和附表-附表2.xlsx")
gtex_ss$Samples

gtex_folder <- dir("F:/post-GWAS/data_smr/GTEx_V8_cis_eqtl_summary_hg19_48G/",full.names = T)
gtex_list <- data.frame(folder = gtex_folder,
                        Tissue = basename(gtex_folder))

gtex_list <- merge(gtex_list,gtex_ss[,c("Tissue","Samples")],by = "Tissue")


######## 共定位分析
### 准备用于查询的探针：对应.epi文件的v2列。
ts8 <- read.table(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/ER+_Breast_cancer/overlap_gene24.txt")
probe <- Ensembl_GRCh37[Ensembl_GRCh37$gene_name %in% ts8$V1,"gene_id"]

### 使用for循环，完成批量共定位
for (i in 1:nrow(gtex_list)) {
  print(i)
  coloc_batch_smr(probes = probe,
                  smr_exe_path = "E:/药物靶点肺癌/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
                  qtl_path = gtex_list$folder[i], # 这里传对应文件夹的路径即可
                  cis_wind = 1000,
                  
                  # prepare_method=2,按照SNP的rsid选择SNP，以GWAS1的染色体坐标为参考，不使用GWAS2的chr和pos信息。
                  prepare_method = 2,
                  type1 = "quant",
                  SS1 = gtex_list$Samples[i],
                  gwas_path2 = "F:/post-GWAS/data_gwas_list/data_prepare/coloc/ER+_Breast_cancer_BCACIEU_hg38.txt",
                  type2 = "cc",
                  SS2 = 175475, # 对应sample.size
                  NC2 = 69501,        # 对应ncase
                  bfile = "E:/药物靶点肺癌/配套/EUR",
                  coloc_plot = T,
                  coloc_plot_pph4 = 0.6,
                  coloc_plot_genome = "hg37",
                  out_path = paste0("E:/post_gwas_cancer/results/ER+_Breast_cancer_BCACIEU_coloc/",gtex_list$Tissue[i]))
}


######## 结果整理
# 获取组织对应的文件夹
sum_folder <- dir("E:/post_gwas_cancer/results/ER+_Breast_cancer_BCACIEU_coloc/",full.names = T,include.dirs = T)

# 汇总组织和对应的结果文件
sum_dat <- data.frame()
for (f in sum_folder) {
  
  # 匹配文件夹和组织信息
  f_fn <- dir(f,pattern = "coloc_summary_coloc_abf.csv",recursive = T,full.names = T)
  if(length(f_fn)==0){
    next
  }
  tissues <- basename(f)
  
  # 获取coloc结果，并增加组织信息
  temp_dat <- read.csv(f_fn)
  temp_dat$tissues <- tissues
  sum_dat <- rbind(sum_dat,temp_dat)
}

# 将汇总后的结果写出
write.csv(sum_dat,file = "E:/post_gwas_cancer/results/ER+_Breast_cancer_BCACIEU_coloc/ER+_Breast_cancer_BCACIEU汇总的共定位结果.csv",row.names = F)


################ Melanoma_FinnGenR10 ####################
##########################  如果是跑全部组织的共定位分析，可以使用如下代码：

######## 输入文件和样本量信息准备
## 对应的gtex文件，增加样本量信息
gtex_ss <- read.xlsx("F:/post-GWAS/data_smr/GTEx_v8_science_2020 附图和附表-附表2.xlsx")
gtex_ss$Samples

gtex_folder <- dir("F:/post-GWAS/data_smr/GTEx_V8_cis_eqtl_summary_hg19_48G/",full.names = T)
gtex_list <- data.frame(folder = gtex_folder,
                        Tissue = basename(gtex_folder))

gtex_list <- merge(gtex_list,gtex_ss[,c("Tissue","Samples")],by = "Tissue")


######## 共定位分析
### 准备用于查询的探针：对应.epi文件的v2列。
ts8 <- read.table(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/黑色素瘤FinnGen/overlap_gene3.txt")
probe <- Ensembl_GRCh37[Ensembl_GRCh37$gene_name %in% ts8$V1,"gene_id"]

### 使用for循环，完成批量共定位
for (i in 1:nrow(gtex_list)) {
  print(i)
  coloc_batch_smr(probes = probe,
                  smr_exe_path = "E:/药物靶点肺癌/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
                  qtl_path = gtex_list$folder[i], # 这里传对应文件夹的路径即可
                  cis_wind = 1000,
                  
                  # prepare_method=2,按照SNP的rsid选择SNP，以GWAS1的染色体坐标为参考，不使用GWAS2的chr和pos信息。
                  prepare_method = 2,
                  type1 = "quant",
                  SS1 = gtex_list$Samples[i],
                  gwas_path2 = "F:/post-GWAS/data_gwas_list/data_prepare/coloc/Melanoma_FinnGenR10.txt",
                  type2 = "cc",
                  SS2 = 318158, # 对应sample.size
                  NC2 = 4261,        # 对应ncase
                  bfile = "E:/药物靶点肺癌/配套/EUR",
                  coloc_plot = T,
                  coloc_plot_pph4 = 0.6,
                  coloc_plot_genome = "hg37",
                  out_path = paste0("E:/post_gwas_cancer/results/tmp_coloc/",gtex_list$Tissue[i]))
}


######## 结果整理
# 获取组织对应的文件夹
sum_folder <- dir("E:/post_gwas_cancer/results/tmp_coloc/",full.names = T,include.dirs = T)

# 汇总组织和对应的结果文件
sum_dat <- data.frame()
for (f in sum_folder) {
  
  # 匹配文件夹和组织信息
  f_fn <- dir(f,pattern = "coloc_summary_coloc_abf.csv",recursive = T,full.names = T)
  if(length(f_fn)==0){
    next
  }
  tissues <- basename(f)
  
  # 获取coloc结果，并增加组织信息
  temp_dat <- read.csv(f_fn)
  temp_dat$tissues <- tissues
  sum_dat <- rbind(sum_dat,temp_dat)
}

# 将汇总后的结果写出
write.csv(sum_dat,file = "E:/post_gwas_cancer/results/tmp_coloc/Melanoma_FinnGenR10汇总的共定位结果.csv",row.names = F)

################ Prostate_cancer_FinnGenR10 ####################
##########################  如果是跑全部组织的共定位分析，可以使用如下代码：

######## 输入文件和样本量信息准备
## 对应的gtex文件，增加样本量信息
gtex_ss <- read.xlsx("F:/post-GWAS/data_smr/GTEx_v8_science_2020 附图和附表-附表2.xlsx")
gtex_ss$Samples

gtex_folder <- dir("F:/post-GWAS/data_smr/GTEx_V8_cis_eqtl_summary_hg19_48G/",full.names = T)
gtex_list <- data.frame(folder = gtex_folder,
                        Tissue = basename(gtex_folder))

gtex_list <- merge(gtex_list,gtex_ss[,c("Tissue","Samples")],by = "Tissue")


######## 共定位分析
### 准备用于查询的探针：对应.epi文件的v2列。
ts8 <- read.table(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/前列腺癌/overlap_gene16.txt")
probe <- Ensembl_GRCh37[Ensembl_GRCh37$gene_name %in% ts8$V1,"gene_id"]

### 使用for循环，完成批量共定位
for (i in 1:nrow(gtex_list)) {
  print(i)
  coloc_batch_smr(probes = probe,
                  smr_exe_path = "E:/药物靶点肺癌/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
                  qtl_path = gtex_list$folder[i], # 这里传对应文件夹的路径即可
                  cis_wind = 1000,
                  
                  # prepare_method=2,按照SNP的rsid选择SNP，以GWAS1的染色体坐标为参考，不使用GWAS2的chr和pos信息。
                  prepare_method = 2,
                  type1 = "quant",
                  SS1 = gtex_list$Samples[i],
                  gwas_path2 = "F:/post-GWAS/data_gwas_list/data_prepare/coloc/Prostate_cancer_FinnGenR10.txt",
                  type2 = "cc",
                  SS2 = 146465, # 对应sample.size
                  NC2 = 15199,        # 对应ncase
                  bfile = "E:/药物靶点肺癌/配套/EUR",
                  coloc_plot = T,
                  coloc_plot_pph4 = 0.6,
                  coloc_plot_genome = "hg37",
                  out_path = paste0("E:/post_gwas_cancer/results/tmp_coloc/",gtex_list$Tissue[i]))
}


######## 结果整理
# 获取组织对应的文件夹
sum_folder <- dir("E:/post_gwas_cancer/results/tmp_coloc/",full.names = T,include.dirs = T)

# 汇总组织和对应的结果文件
sum_dat <- data.frame()
for (f in sum_folder) {
  
  # 匹配文件夹和组织信息
  f_fn <- dir(f,pattern = "coloc_summary_coloc_abf.csv",recursive = T,full.names = T)
  if(length(f_fn)==0){
    next
  }
  tissues <- basename(f)
  
  # 获取coloc结果，并增加组织信息
  temp_dat <- read.csv(f_fn)
  temp_dat$tissues <- tissues
  sum_dat <- rbind(sum_dat,temp_dat)
}

# 将汇总后的结果写出
write.csv(sum_dat,file = "E:/post_gwas_cancer/results/tmp_coloc/Prostate_cancer_FinnGenR10汇总的共定位结果.csv",row.names = F)

