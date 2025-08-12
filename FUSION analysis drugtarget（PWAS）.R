### FUSION分析 ###
library(DrugTargetMR)

#### EA #####
# step1. -----------------------------------------------------------
prepare_sumstat <- data.table::fread("F:/post-GWAS/data_gwas_list/data_prepare/fusion/Breast_cancer_BCACIEU_hg38.sumstats",
                                     data.table = F)

fusion_assoc(sumstat_data = prepare_sumstat,
             weights = "F:/post-GWAS/data_fusion/权重文件/EA和AA血液蛋白权重/PWAS_EA/Plasma_Protein_EA_hg38.pos",
             weights_dir = "F:/post-GWAS/data_fusion/权重文件/EA和AA血液蛋白权重/PWAS_EA/Plasma_Protein_weights_EA/",
             ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
             ref_ld_chr_num = 1:22,
             out_path = "E:/post_gwas_cancer/tmp/PWAS/",
             out_prefix = "Breast_cancer_PWAS" )

### 结果整理
res_fusion_files <- dir("E:/post_gwas_cancer/tmp/PWAS/",pattern = ".csv",recursive = F,full.names = T)
res_fusion <- data.frame()

for (i in res_fusion_files) {
  
  if(!file.exists(i)){
    next
  }
  
  temp_res <- read.csv(i)
  temp_res <- temp_res[!is.na(temp_res$TWAS.P),]
  temp_res$TWAS.P.fdr <- p.adjust(temp_res$TWAS.P,method = "fdr")
  
  if(nrow(temp_res)>0){
    res_fusion <- rbind(res_fusion,temp_res)
  }
}

# 保留TWAS.P.fdr<0.05对应的行，并将结果写出
res_fusion_f <- res_fusion[res_fusion$TWAS.P.fdr < 0.05,]
res_fusion_f$gene_id <- substr(res_fusion_f$ID,1,regexpr("\\.",res_fusion_f$ID)-1)
res_fusion_f <- merge(res_fusion_f,Ensembl_GRCh38,by="gene_id",all.x=T)

write.csv(res_fusion,file = "E:/post_gwas_cancer/tmp/PWAS/Breast_cancer_PWAS_summary_p_fdr.csv")
write.csv(res_fusion_f,file = "E:/post_gwas_cancer/tmp/PWAS/Breast_cancer_PWAS_summary_p_fdr_filtered.csv")

# step2. -----------------------------------------------------------
prepare_sumstat <- data.table::fread("F:/post-GWAS/data_gwas_list/data_prepare/fusion/ER-_Breast_cancer_BCACIEU_hg38.sumstats",
                                     data.table = F)

fusion_assoc(sumstat_data = prepare_sumstat,
             weights = "F:/post-GWAS/data_fusion/权重文件/EA和AA血液蛋白权重/PWAS_EA/Plasma_Protein_EA_hg38.pos",
             weights_dir = "F:/post-GWAS/data_fusion/权重文件/EA和AA血液蛋白权重/PWAS_EA/Plasma_Protein_weights_EA/",
             ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
             ref_ld_chr_num = 1:22,
             out_path = "E:/post_gwas_cancer/tmp/PWAS/",
             out_prefix = "ER-_Breast_cancer_PWAS" )

### 结果整理
res_fusion_files <- dir("E:/post_gwas_cancer/tmp/PWAS/",pattern = ".csv",recursive = F,full.names = T)
res_fusion <- data.frame()

for (i in res_fusion_files) {
  
  if(!file.exists(i)){
    next
  }
  
  temp_res <- read.csv(i)
  temp_res <- temp_res[!is.na(temp_res$TWAS.P),]
  temp_res$TWAS.P.fdr <- p.adjust(temp_res$TWAS.P,method = "fdr")
  
  if(nrow(temp_res)>0){
    res_fusion <- rbind(res_fusion,temp_res)
  }
}

# 保留TWAS.P.fdr<0.05对应的行，并将结果写出
res_fusion_f <- res_fusion[res_fusion$TWAS.P.fdr < 0.05,]
res_fusion_f$gene_id <- substr(res_fusion_f$ID,1,regexpr("\\.",res_fusion_f$ID)-1)
res_fusion_f <- merge(res_fusion_f,Ensembl_GRCh38,by="gene_id",all.x=T)

write.csv(res_fusion,file = "E:/post_gwas_cancer/tmp/PWAS/ER-_Breast_cancer_PWAS_summary_p_fdr.csv")
write.csv(res_fusion_f,file = "E:/post_gwas_cancer/tmp/PWAS/ER-_Breast_cancer_PWAS_summary_p_fdr_filtered.csv")

# step3. -----------------------------------------------------------
prepare_sumstat <- data.table::fread("F:/post-GWAS/data_gwas_list/data_prepare/fusion/ER+_Breast_cancer_BCACIEU_hg38.sumstats",
                                     data.table = F)

fusion_assoc(sumstat_data = prepare_sumstat,
             weights = "F:/post-GWAS/data_fusion/权重文件/EA和AA血液蛋白权重/PWAS_EA/Plasma_Protein_EA_hg38.pos",
             weights_dir = "F:/post-GWAS/data_fusion/权重文件/EA和AA血液蛋白权重/PWAS_EA/Plasma_Protein_weights_EA/",
             ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
             ref_ld_chr_num = 1:22,
             out_path = "E:/post_gwas_cancer/tmp/PWAS/",
             out_prefix = "ER+_Breast_cancer_PWAS" )

### 结果整理
res_fusion_files <- dir("E:/post_gwas_cancer/tmp/PWAS/",pattern = ".csv",recursive = F,full.names = T)
res_fusion <- data.frame()

for (i in res_fusion_files) {
  
  if(!file.exists(i)){
    next
  }
  
  temp_res <- read.csv(i)
  temp_res <- temp_res[!is.na(temp_res$TWAS.P),]
  temp_res$TWAS.P.fdr <- p.adjust(temp_res$TWAS.P,method = "fdr")
  
  if(nrow(temp_res)>0){
    res_fusion <- rbind(res_fusion,temp_res)
  }
}

# 保留TWAS.P.fdr<0.05对应的行，并将结果写出
res_fusion_f <- res_fusion[res_fusion$TWAS.P.fdr < 0.05,]
res_fusion_f$gene_id <- substr(res_fusion_f$ID,1,regexpr("\\.",res_fusion_f$ID)-1)
res_fusion_f <- merge(res_fusion_f,Ensembl_GRCh38,by="gene_id",all.x=T)

write.csv(res_fusion,file = "E:/post_gwas_cancer/tmp/PWAS/ER+_Breast_cancer_PWAS_summary_p_fdr.csv")
write.csv(res_fusion_f,file = "E:/post_gwas_cancer/tmp/PWAS/ER+_Breast_cancer_PWAS_summary_p_fdr_filtered.csv")


#### AA #####
# step1. -----------------------------------------------------------
prepare_sumstat <- data.table::fread("F:/post-GWAS/data_gwas_list/data_prepare/fusion/Breast_cancer_BCACIEU_hg38.sumstats",
                                     data.table = F)

fusion_assoc(sumstat_data = prepare_sumstat,
             weights = "F:/post-GWAS/data_fusion/权重文件/EA和AA血液蛋白权重/PWAS_AA/Plasma_Protein_AA_hg38.pos",
             weights_dir = "F:/post-GWAS/data_fusion/权重文件/EA和AA血液蛋白权重/PWAS_AA/Plasma_Protein_weights_AA/",
             ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
             ref_ld_chr_num = 1:22,
             out_path = "E:/post_gwas_cancer/tmp/PWAS/",
             out_prefix = "Breast_cancer_PWAS" )

### 结果整理
res_fusion_files <- dir("E:/post_gwas_cancer/tmp/PWAS/",pattern = ".csv",recursive = F,full.names = T)
res_fusion <- data.frame()

for (i in res_fusion_files) {
  
  if(!file.exists(i)){
    next
  }
  
  temp_res <- read.csv(i)
  temp_res <- temp_res[!is.na(temp_res$TWAS.P),]
  temp_res$TWAS.P.fdr <- p.adjust(temp_res$TWAS.P,method = "fdr")
  
  if(nrow(temp_res)>0){
    res_fusion <- rbind(res_fusion,temp_res)
  }
}

# 保留TWAS.P.fdr<0.05对应的行，并将结果写出
res_fusion_f <- res_fusion[res_fusion$TWAS.P.fdr < 0.05,]
res_fusion_f$gene_id <- substr(res_fusion_f$ID,1,regexpr("\\.",res_fusion_f$ID)-1)
res_fusion_f <- merge(res_fusion_f,Ensembl_GRCh38,by="gene_id",all.x=T)

write.csv(res_fusion,file = "E:/post_gwas_cancer/tmp/PWAS/Breast_cancer_PWAS_summary_p_fdr.csv")
write.csv(res_fusion_f,file = "E:/post_gwas_cancer/tmp/PWAS/Breast_cancer_PWAS_summary_p_fdr_filtered.csv")

# step2. -----------------------------------------------------------
prepare_sumstat <- data.table::fread("F:/post-GWAS/data_gwas_list/data_prepare/fusion/ER-_Breast_cancer_BCACIEU_hg38.sumstats",
                                     data.table = F)

fusion_assoc(sumstat_data = prepare_sumstat,
             weights = "F:/post-GWAS/data_fusion/权重文件/EA和AA血液蛋白权重/PWAS_AA/Plasma_Protein_AA_hg38.pos",
             weights_dir = "F:/post-GWAS/data_fusion/权重文件/EA和AA血液蛋白权重/PWAS_AA/Plasma_Protein_weights_AA/",
             ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
             ref_ld_chr_num = 1:22,
             out_path = "E:/post_gwas_cancer/tmp/PWAS/",
             out_prefix = "ER-_Breast_cancer_PWAS" )

### 结果整理
res_fusion_files <- dir("E:/post_gwas_cancer/tmp/PWAS/",pattern = ".csv",recursive = F,full.names = T)
res_fusion <- data.frame()

for (i in res_fusion_files) {
  
  if(!file.exists(i)){
    next
  }
  
  temp_res <- read.csv(i)
  temp_res <- temp_res[!is.na(temp_res$TWAS.P),]
  temp_res$TWAS.P.fdr <- p.adjust(temp_res$TWAS.P,method = "fdr")
  
  if(nrow(temp_res)>0){
    res_fusion <- rbind(res_fusion,temp_res)
  }
}

# 保留TWAS.P.fdr<0.05对应的行，并将结果写出
res_fusion_f <- res_fusion[res_fusion$TWAS.P.fdr < 0.05,]
res_fusion_f$gene_id <- substr(res_fusion_f$ID,1,regexpr("\\.",res_fusion_f$ID)-1)
res_fusion_f <- merge(res_fusion_f,Ensembl_GRCh38,by="gene_id",all.x=T)

write.csv(res_fusion,file = "E:/post_gwas_cancer/tmp/PWAS/ER-_Breast_cancer_PWAS_summary_p_fdr.csv")
write.csv(res_fusion_f,file = "E:/post_gwas_cancer/tmp/PWAS/ER-_Breast_cancer_PWAS_summary_p_fdr_filtered.csv")

# step3. -----------------------------------------------------------
prepare_sumstat <- data.table::fread("F:/post-GWAS/data_gwas_list/data_prepare/fusion/ER+_Breast_cancer_BCACIEU_hg38.sumstats",
                                     data.table = F)

fusion_assoc(sumstat_data = prepare_sumstat,
             weights = "F:/post-GWAS/data_fusion/权重文件/EA和AA血液蛋白权重/PWAS_AA/Plasma_Protein_AA_hg38.pos",
             weights_dir = "F:/post-GWAS/data_fusion/权重文件/EA和AA血液蛋白权重/PWAS_AA/Plasma_Protein_weights_AA/",
             ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
             ref_ld_chr_num = 1:22,
             out_path = "E:/post_gwas_cancer/tmp/PWAS/",
             out_prefix = "ER+_Breast_cancer_PWAS" )

### 结果整理
res_fusion_files <- dir("E:/post_gwas_cancer/tmp/PWAS/",pattern = ".csv",recursive = F,full.names = T)
res_fusion <- data.frame()

for (i in res_fusion_files) {
  
  if(!file.exists(i)){
    next
  }
  
  temp_res <- read.csv(i)
  temp_res <- temp_res[!is.na(temp_res$TWAS.P),]
  temp_res$TWAS.P.fdr <- p.adjust(temp_res$TWAS.P,method = "fdr")
  
  if(nrow(temp_res)>0){
    res_fusion <- rbind(res_fusion,temp_res)
  }
}

# 保留TWAS.P.fdr<0.05对应的行，并将结果写出
res_fusion_f <- res_fusion[res_fusion$TWAS.P.fdr < 0.05,]
res_fusion_f$gene_id <- substr(res_fusion_f$ID,1,regexpr("\\.",res_fusion_f$ID)-1)
res_fusion_f <- merge(res_fusion_f,Ensembl_GRCh38,by="gene_id",all.x=T)

write.csv(res_fusion,file = "E:/post_gwas_cancer/tmp/PWAS/ER+_Breast_cancer_PWAS_summary_p_fdr.csv")
write.csv(res_fusion_f,file = "E:/post_gwas_cancer/tmp/PWAS/ER+_Breast_cancer_PWAS_summary_p_fdr_filtered.csv")
