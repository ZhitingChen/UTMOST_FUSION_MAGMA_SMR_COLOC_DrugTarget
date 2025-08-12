library(DrugTargetMR)

# dir <- list.files("F:/post-GWAS/data_gwas_list/data_prepare/ldsc/")# GWAS所在文件夹
# data_dir <- "F:/post-GWAS/data_gwas_list/data_prepare/ldsc/"
# file_path <- paste0(data_dir,"",dir)

## 方法1，使用ancestry参数，默认使用来自Pan-UKB对应人群的参考文件
postgwas_ldsc_h2(file_path =c("F:/post-GWAS/data_gwas_list/data_prepare/coloc/Breast_cancer_BCACIEU_hg38.txt",
                              "F:/post-GWAS/data_gwas_list/data_prepare/coloc/ER-_Breast_cancer_BCACIEU_hg38.txt", 
                              "F:/post-GWAS/data_gwas_list/data_prepare/coloc/ER+_Breast_cancer_BCACIEU_hg38.txt",
                              "F:/post-GWAS/data_gwas_list/data_prepare/coloc/finngen_R10_C3_BREAST_EXALLC.txt",
                              "F:/post-GWAS/data_gwas_list/data_prepare/coloc/finngen_R10_C3_BREAST_ERNEG_EXALLC.txt", 
                              "F:/post-GWAS/data_gwas_list/data_prepare/coloc/finngen_R10_C3_BREAST_ERPLUS_EXALLC.txt"),
                 ancestry = "EUR",
                 out_path = "E:/post_gwas_cancer/tmp2/")

## 方法2，使用ld和wld参数，指定对应人群的参考文件，包括.l2.ldscore.gz文件。
postgwas_ldsc_h2(file_path =c("F:/post-GWAS/data_gwas_list/data_prepare/coloc/Breast_cancer_BCACIEU_hg38.txt",
                              "F:/post-GWAS/data_gwas_list/data_prepare/coloc/ER-_Breast_cancer_BCACIEU_hg38.txt", 
                              "F:/post-GWAS/data_gwas_list/data_prepare/coloc/ER+_Breast_cancer_BCACIEU_hg38.txt",
                              "F:/post-GWAS/data_gwas_list/data_prepare/coloc/finngen_R10_C3_BREAST_EXALLC.txt",
                              "F:/post-GWAS/data_gwas_list/data_prepare/coloc/finngen_R10_C3_BREAST_ERNEG_EXALLC.txt", 
                              "F:/post-GWAS/data_gwas_list/data_prepare/coloc/finngen_R10_C3_BREAST_ERPLUS_EXALLC.txt"),
                 ld = "F:/post-GWAS/data_ldsc/eur_w_ld_chr/",
                 wld = "F:/post-GWAS/data_ldsc/eur_w_ld_chr/",
                 out_path = "E:/post_gwas_cancer/tmp3/")



### Gene 遗传力计算 ###
file_path <- "F:/28849TLR1.vcf.gz"
prepare_ieu(
  file_path,
  out_path = "F:/tmp/",
  generate_mr = T,
  generate_smr = F,
  sample_size = 28849
)

## 方法1，使用ancestry参数，默认使用来自Pan-UKB对应人群的参考文件
postgwas_ldsc_h2(file_path =c("E:/post_gwas_cancer/tmp3/26609ADCY3.txt",
                              "E:/post_gwas_cancer/tmp3/31470CASP8.txt", 
                              "E:/post_gwas_cancer/tmp3/25860GRHL1.txt",
                              "E:/post_gwas_cancer/tmp3/14263HELQ.txt",
                              "E:/post_gwas_cancer/tmp3/28849TLR1.txt"),
                 ancestry = "EUR",
                 out_path = "E:/post_gwas_cancer/tmp2/")

## 方法2，使用ld和wld参数，指定对应人群的参考文件，包括.l2.ldscore.gz文件。
postgwas_ldsc_h2(file_path =c("E:/post_gwas_cancer/tmp3/26609ADCY3.txt",
                              "E:/post_gwas_cancer/tmp3/31470CASP8.txt", 
                              "E:/post_gwas_cancer/tmp3/25860GRHL1.txt",
                              "E:/post_gwas_cancer/tmp3/14263HELQ.txt",
                              "E:/post_gwas_cancer/tmp3/28849TLR1.txt"),
                 ld = "F:/post-GWAS/data_ldsc/eur_w_ld_chr/",
                 wld = "F:/post-GWAS/data_ldsc/eur_w_ld_chr/",
                 out_path = "E:/post_gwas_cancer/tmp3/")

