### UTMOST分析 ###
library(DrugTargetMR)

# 循环语句进行 

dir <- list.files("F:/post-GWAS/data_gwas_list/data_prepare/coloc/")# GWAS所在文件夹
data_dir <- "F:/post-GWAS/data_gwas_list/data_prepare/coloc/"
n <- length(dir)
file_path <- paste0(data_dir,"",dir)


for (i in 15:n) {

  prefix <- sub("\\.txt$", "", dir[i])
  out_prefix1 <- paste0(prefix,"_utmost")
  out_prefix2 <- paste0("utmost_",prefix)
  out_prefix3 <- paste0("E:/post_gwas_cancer/results/utmost/Final_results_wtih_padj/utmost_",prefix,"_joint_GBJ_filter_p.csv")
  out_prefix4 <- paste0("E:/post_gwas_cancer/results/utmost/Final_results_wtih_padj/utmost_",prefix,"_joint_GBJ_filter_adj_p.csv")
  results_file <- paste0("E:/post_gwas_cancer/results/utmost/Final_results_GBJ/utmost_",prefix,"_joint_GBJ.csv")
  gene_info_file <- paste0("E:/post_gwas_cancer/results/utmost/",prefix,"_utmost_summary_gene_info.txt")
  input_file <- paste0("E:/post_gwas_cancer/results/utmost/",prefix,"_utmost_summary_results.csv")
  
  utmost <- utmost_single_tissue(gwas_file = file_path[i],
                               cov_folder  = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/covariance_GTEx8_normalized_pruned/",
                               weight_db_folder  = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/database_normalized_pruned/",
                               # workers = 12,
                               # single_tissue = c("Adipose_Subcutaneous","Artery_Tibial"),  # 默认使用指定文件夹下的所有.db文件
                               out_path = "E:/post_gwas_cancer/results/utmost/",
                               out_prefix = out_prefix1)

  utmost_joint_GBJ(gene_info_file = gene_info_file,
                 input_file = input_file,
                 cov_folder = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/covariance_joint_GTEx8_normalized_pruned/",
                 weight_db_folder = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/database_normalized_pruned/",
                 out_path = "E:/post_gwas_cancer/results/utmost/Final_results_GBJ/",
                 out_prefix = out_prefix2)
  
  
  cross_res <- data.table::fread(results_file,data.table = F)
  cross_res <- cross_res[!is.na(cross_res$p_value),]
  cross_res$ENSE <- substr(cross_res$gene,1,regexpr("\\.",cross_res$gene)-1)
  cross_res <- merge(cross_res,Ensembl_GRCh38,by.x="ENSE",by.y="gene_id",all.x = T)
  cross_res <- cross_res[!duplicated(cross_res$ENSE),]
  cross_res$p_fdr <- p.adjust(cross_res$p_value,method = "fdr")
  cross_res_f1  <- cross_res[cross_res$p_value < 0.05,]
  cross_res_f2  <- cross_res[cross_res$p_fdr < 0.05,]
  write.csv(cross_res_f1,file = out_prefix3,row.names = F)
  write.csv(cross_res_f2,file = out_prefix4,row.names = F)
  
}





# step1: 执行单组织utmost分析 #
# 说明1：gwas_file，传入对应的GWAS文件，为MR标准格式。程序已经内置了数据清洗，无需再次清洗数据。
# 说明2：utmost_single_tissue执行单组织utmost分析：这一步将按照组织进行结果数据，并汇总展示合并后的结果和基因信息
# 说明3：查询函数用法，请执行：?utmost_single_tissue
utmost <- utmost_single_tissue(gwas_file = "F:/post-GWAS/data_gwas_list/data_prepare/mr/Bladder_cancer_GWASCatalog.txt",
                               cov_folder  = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/covariance_GTEx8_normalized_pruned/",
                               weight_db_folder  = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/database_normalized_pruned/",
                               # single_tissue = c("Adipose_Subcutaneous","Artery_Tibial"),  # 默认使用指定文件夹下的所有.db文件
                               out_path = "E:/post_gwas_cancer/results/utmost/",
                               out_prefix = "Bladder_cancer_GWASCatalog_utmost")

# step2: 执行跨组织utmost分析 #
# 说明1：utmost_joint_GBJ执行跨组织utmost分析，这一步使用到的输入数据，为单组织utmost的输出文件
utmost_joint_GBJ(gene_info_file = "E:/post_gwas_cancer/results/utmost/Bladder_cancer_GWASCatalog_utmost_summary_gene_info.txt",
                 input_file = "E:/post_gwas_cancer/results/utmost/Bladder_cancer_GWASCatalog_utmost_summary_results.csv",
                 cov_folder = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/covariance_joint_GTEx8_normalized_pruned/",
                 weight_db_folder = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/database_normalized_pruned/",
                 out_path = "E:/post_gwas_cancer/results/utmost/GBJ/",
                 out_prefix ="utmost_Bladder_cancer_GWASCatalog" )


############ step3: 结果整理 ############
### 读取跨组织utmost分析结果
cross_res <- data.table::fread("E:/post_gwas_cancer/results/utmost/Final_results_GBJ/utmost_Breast_cancer_BCACIEU_joint_GBJ.csv",data.table = F)

## 去除缺失统计信息的行
cross_res <- cross_res[!is.na(cross_res$p_value),]

## 增加列：ENSE(ensemble id)
cross_res$ENSE <- substr(cross_res$gene,1,regexpr("\\.",cross_res$gene)-1)

## 匹配gene symbol，去重
cross_res <- merge(cross_res,Ensembl_GRCh37,by.x="ENSE",by.y="gene_id",all.x = T)
cross_res <- cross_res[!duplicated(cross_res$ENSE),]

## p值校正
cross_res$p_fdr <- p.adjust(cross_res$p_value,method = "fdr")

## 未校正p值，具有显著性的gene
cross_res_f1  <- cross_res[cross_res$p_value < 0.05,]

## 校正p值，具有显著性的gene
cross_res_f2  <- cross_res[cross_res$p_fdr < 0.05,]

# 将结果写出
write.csv(cross_res_f1,file = out_prefix3,row.names = F)
write.csv(cross_res_f2,file = out_prefix4,row.names = F)
