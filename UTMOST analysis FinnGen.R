library(DrugTargetMR)

#### Breast Cancer ####
############ step1: 执行单组织utmost分析 ############
# 说明1：gwas_file，传入对应的GWAS文件，为MR标准格式。程序已经内置了数据清洗，无需再次清洗数据。
# 说明2：utmost_single_tissue执行单组织utmost分析：这一步将按照组织进行结果数据，并汇总展示合并后的结果和基因信息
# 说明3：查询函数用法，请执行：?utmost_single_tissue

utmost <- utmost_single_tissue(gwas_file = "F:/tmp/finngen_R10_C3_BREAST_EXALLC.txt",
                               cov_folder  = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/covariance_GTEx8_normalized_pruned/",
                               weight_db_folder  = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/database_normalized_pruned/",
                               # workers = 12,
                               # single_tissue = c("Adipose_Subcutaneous","Artery_Tibial"),  # 默认使用指定文件夹下的所有.db文件
                               out_path = "E:/post_gwas_cancer/tmp1/",
                               out_prefix = "utmost")

############ step2: 执行跨组织utmost分析 ############
# 说明1：utmost_joint_GBJ执行跨组织utmost分析，这一步使用到的输入数据，为单组织utmost的输出文件
utmost_joint_GBJ(gene_info_file = "E:/post_gwas_cancer/tmp1/utmost_summary_gene_info.txt",
                 input_file = "E:/post_gwas_cancer/tmp1/utmost_summary_results.csv",
                 cov_folder = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/covariance_joint_GTEx8_normalized_pruned/",
                 weight_db_folder = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/database_normalized_pruned/",
                 out_path = "E:/post_gwas_cancer/tmp1/Breast_Cancer_FinnGen/",
                 out_prefix ="utmost_Breast_Cancer_FinnGen" )

############ step3: 结果整理 ############
### 读取跨组织utmost分析结果
cross_res <- data.table::fread("E:/post_gwas_cancer/tmp1/Breast_Cancer_FinnGen/utmost_Breast_Cancer_FinnGen_joint_GBJ.csv",data.table = F)

## 去除缺失统计信息的行
cross_res <- cross_res[!is.na(cross_res$p_value),]

## 增加列：ENSE(ensemble id)
cross_res$ENSE <- substr(cross_res$gene,1,regexpr("\\.",cross_res$gene)-1)

## 匹配gene symbol，去重
cross_res <- merge(cross_res,Ensembl_GRCh38,by.x="ENSE",by.y="gene_id",all.x = T)
cross_res <- cross_res[!duplicated(cross_res$ENSE),]

## p值校正
cross_res$p_fdr <- p.adjust(cross_res$p_value,method = "fdr")

## 未校正p值，具有显著性的gene
cross_res_f1  <- cross_res[cross_res$p_value < 0.05,]

## 校正p值，具有显著性的gene
cross_res_f2  <- cross_res[cross_res$p_fdr < 0.05,]

# 将结果写出
write.csv(cross_res_f1,file = "E:/post_gwas_cancer/tmp1/Breast_Cancer_FinnGen/utmost_Breast_Cancer_FinnGen_joint_GBJ_filter_p.csv",row.names = F)
write.csv(cross_res_f2,file = "E:/post_gwas_cancer/tmp1/Breast_Cancer_FinnGen/utmost_Breast_Cancer_FinnGen_joint_GBJ_filter_adj_p.csv",row.names = F)


#### ER- Breast Cancer ####
############ step1: 执行单组织utmost分析 ############
# 说明1：gwas_file，传入对应的GWAS文件，为MR标准格式。程序已经内置了数据清洗，无需再次清洗数据。
# 说明2：utmost_single_tissue执行单组织utmost分析：这一步将按照组织进行结果数据，并汇总展示合并后的结果和基因信息
# 说明3：查询函数用法，请执行：?utmost_single_tissue

utmost <- utmost_single_tissue(gwas_file = "F:/tmp/finngen_R10_C3_BREAST_ERNEG_EXALLC.txt",
                               cov_folder  = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/covariance_GTEx8_normalized_pruned/",
                               weight_db_folder  = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/database_normalized_pruned/",
                               # workers = 12,
                               # single_tissue = c("Adipose_Subcutaneous","Artery_Tibial"),  # 默认使用指定文件夹下的所有.db文件
                               out_path = "E:/post_gwas_cancer/tmp1/",
                               out_prefix = "utmost")

############ step2: 执行跨组织utmost分析 ############
# 说明1：utmost_joint_GBJ执行跨组织utmost分析，这一步使用到的输入数据，为单组织utmost的输出文件
utmost_joint_GBJ(gene_info_file = "E:/post_gwas_cancer/tmp1/utmost_summary_gene_info.txt",
                 input_file = "E:/post_gwas_cancer/tmp1/utmost_summary_results.csv",
                 cov_folder = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/covariance_joint_GTEx8_normalized_pruned/",
                 weight_db_folder = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/database_normalized_pruned/",
                 out_path = "E:/post_gwas_cancer/tmp1/ER-_Breast_Cancer_FinnGen/",
                 out_prefix ="utmost_ER-_Breast_Cancer_FinnGen" )

############ step3: 结果整理 ############
### 读取跨组织utmost分析结果
cross_res <- data.table::fread("E:/post_gwas_cancer/tmp1/ER-_Breast_Cancer_FinnGen/utmost_ER-_Breast_Cancer_FinnGen_joint_GBJ.csv",data.table = F)

## 去除缺失统计信息的行
cross_res <- cross_res[!is.na(cross_res$p_value),]

## 增加列：ENSE(ensemble id)
cross_res$ENSE <- substr(cross_res$gene,1,regexpr("\\.",cross_res$gene)-1)

## 匹配gene symbol，去重
cross_res <- merge(cross_res,Ensembl_GRCh38,by.x="ENSE",by.y="gene_id",all.x = T)
cross_res <- cross_res[!duplicated(cross_res$ENSE),]

## p值校正
cross_res$p_fdr <- p.adjust(cross_res$p_value,method = "fdr")

## 未校正p值，具有显著性的gene
cross_res_f1  <- cross_res[cross_res$p_value < 0.05,]

## 校正p值，具有显著性的gene
cross_res_f2  <- cross_res[cross_res$p_fdr < 0.05,]

# 将结果写出
write.csv(cross_res_f1,file = "E:/post_gwas_cancer/tmp1/ER-_Breast_Cancer_FinnGen/utmost_ER-_Breast_Cancer_FinnGen_joint_GBJ_filter_p.csv",row.names = F)
write.csv(cross_res_f2,file = "E:/post_gwas_cancer/tmp1/ER-_Breast_Cancer_FinnGen/utmost_ER-_Breast_Cancer_FinnGen_joint_GBJ_filter_adj_p.csv",row.names = F)

#### ER+ Breast Cancer ####
############ step1: 执行单组织utmost分析 ############
# 说明1：gwas_file，传入对应的GWAS文件，为MR标准格式。程序已经内置了数据清洗，无需再次清洗数据。
# 说明2：utmost_single_tissue执行单组织utmost分析：这一步将按照组织进行结果数据，并汇总展示合并后的结果和基因信息
# 说明3：查询函数用法，请执行：?utmost_single_tissue

utmost <- utmost_single_tissue(gwas_file = "F:/tmp/finngen_R10_C3_BREAST_ERPLUS_EXALLC.txt",
                               cov_folder  = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/covariance_GTEx8_normalized_pruned/",
                               weight_db_folder  = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/database_normalized_pruned/",
                               # workers = 12,
                               # single_tissue = c("Adipose_Subcutaneous","Artery_Tibial"),  # 默认使用指定文件夹下的所有.db文件
                               out_path = "E:/post_gwas_cancer/tmp1/",
                               out_prefix = "utmost")

############ step2: 执行跨组织utmost分析 ############
# 说明1：utmost_joint_GBJ执行跨组织utmost分析，这一步使用到的输入数据，为单组织utmost的输出文件
utmost_joint_GBJ(gene_info_file = "E:/post_gwas_cancer/tmp1/utmost_summary_gene_info.txt",
                 input_file = "E:/post_gwas_cancer/tmp1/utmost_summary_results.csv",
                 cov_folder = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/covariance_joint_GTEx8_normalized_pruned/",
                 weight_db_folder = "F:/post-GWAS/data_utmost/GTEX_v8_20230609/database_normalized_pruned/",
                 out_path = "E:/post_gwas_cancer/tmp1/ER+_Breast_Cancer_FinnGen/",
                 out_prefix ="utmost_ER+_Breast_Cancer_FinnGen" )

############ step3: 结果整理 ############
### 读取跨组织utmost分析结果
cross_res <- data.table::fread("E:/post_gwas_cancer/tmp1/ER+_Breast_Cancer_FinnGen/utmost_ER+_Breast_Cancer_FinnGen_joint_GBJ.csv",data.table = F)

## 去除缺失统计信息的行
cross_res <- cross_res[!is.na(cross_res$p_value),]

## 增加列：ENSE(ensemble id)
cross_res$ENSE <- substr(cross_res$gene,1,regexpr("\\.",cross_res$gene)-1)

## 匹配gene symbol，去重
cross_res <- merge(cross_res,Ensembl_GRCh38,by.x="ENSE",by.y="gene_id",all.x = T)
cross_res <- cross_res[!duplicated(cross_res$ENSE),]

## p值校正
cross_res$p_fdr <- p.adjust(cross_res$p_value,method = "fdr")

## 未校正p值，具有显著性的gene
cross_res_f1  <- cross_res[cross_res$p_value < 0.05,]

## 校正p值，具有显著性的gene
cross_res_f2  <- cross_res[cross_res$p_fdr < 0.05,]

# 将结果写出
write.csv(cross_res_f1,file = "E:/post_gwas_cancer/tmp1/ER+_Breast_Cancer_FinnGen/utmost_ER+_Breast_Cancer_FinnGen_joint_GBJ_filter_p.csv",row.names = F)
write.csv(cross_res_f2,file = "E:/post_gwas_cancer/tmp1/ER+_Breast_Cancer_FinnGen/utmost_ER+_Breast_Cancer_FinnGen_joint_GBJ_filter_adj_p.csv",row.names = F)
