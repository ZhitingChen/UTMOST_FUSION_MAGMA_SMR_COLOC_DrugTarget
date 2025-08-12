### FUSION后Coloc分析  ###
library(DrugTargetMR)
library(openxlsx)
library(data.table)
library(clusterProfiler)
library(org.Mm.eg.db)

### 探针数据处理 ###


### eqtl-GWAS ###
### eqtlGene 验证 ###
smr_qtl2gwas(smr_exe_path = "E:/药物靶点肺癌/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
             bfile = "E:/药物靶点肺癌/配套/EUR",
             gwas_path = "F:/post-GWAS/data_gwas_list/data_prepare/smr/Breast_cancer_BCACIEU.ma",
             qtls_path = "E:/药物靶点肺癌/暴露数据/eqtl/cis-eQTL-SMR_20191212/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense",
             probes_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/Breast_cancer_BCACIEU/overlap_gene49.txt",
             smr_multi_snp = F,
             smr_cis_wind = 2000,
             maf = 0.2,
             smr_peqtl = 5.0e-8,
             out_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/SMR结果/",
             out_prefix = "eqtlgen_Breast_cancer_BCACIEU")
### GTEx 验证（血+组织） ###
smr_qtl2gwas(smr_exe_path = "E:/药物靶点肺癌/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
             bfile = "E:/药物靶点肺癌/配套/EUR",
             gwas_path = c("E:/多组学MR肺癌/IEU肺癌数据/NSCLC/NSCLC_LEU.ma",
                           "E:/多组学MR肺癌/IEU肺癌数据/NSCLC/NSCLC_UKB.ma",
                           "E:/多组学MR肺癌/IEU肺癌数据/LUAD/LUAD.ma",
                           "E:/多组学MR肺癌/IEU肺癌数据/LUSC/LUSC.ma",
                           "E:/多组学MR肺癌/IEU肺癌数据/SCLC/SCLC.ma"),
             qtls_path = "E:/药物靶点肺癌/暴露数据/eqtl/GTEx_V8_cis_eqtl_summary_lite/Whole_Blood.lite",
             probes_path = "E:/多组学MR肺癌/基因集选择/ICD基因集/ICD.txt",
             smr_multi_snp = F,
             smr_cis_wind = 2000,
             maf = 0.2,
             smr_peqtl = 5.0e-8,
             out_path = "E:/多组学MR肺癌/结果/ICD基因集/smr-eqtl/全血验证/",
             out_prefix = "GTXv8_WB")
smr_qtl2gwas(smr_exe_path = "E:/药物靶点肺癌/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
             bfile = "E:/药物靶点肺癌/配套/EUR",
             gwas_path = c("E:/多组学MR肺癌/IEU肺癌数据/NSCLC/NSCLC_LEU.ma",
                           "E:/多组学MR肺癌/IEU肺癌数据/NSCLC/NSCLC_UKB.ma",
                           "E:/多组学MR肺癌/IEU肺癌数据/LUAD/LUAD.ma",
                           "E:/多组学MR肺癌/IEU肺癌数据/LUSC/LUSC.ma",
                           "E:/多组学MR肺癌/IEU肺癌数据/SCLC/SCLC.ma"),
             qtls_path = "E:/药物靶点肺癌/暴露数据/eqtl/GTEx_V8_cis_eqtl_summary_lite/Lung.lite",
             probes_path = "E:/多组学MR肺癌/基因集选择/ICD基因集/ICD.txt",
             smr_multi_snp = F,
             smr_cis_wind = 2000,
             maf = 0.2,
             smr_peqtl = 5.0e-8,
             out_path = "E:/多组学MR肺癌/结果/ICD基因集/smr-eqtl/肺组织验证/",
             out_prefix = "GTXv8_Lung")

# 结果处理
eqtl2gwas <- data.table::fread("E:/多组学MR肺癌/结果/ICD基因集/smr-eqtl/eqtlgen_NSCLC_LEU.smr",data.table = F)
eqtl2gwas$p_SMR_fdr <- p.adjust(eqtl2gwas$p_SMR,method = "fdr")
# eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_HEIDI > 0.01,]
# eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_SMR_fdr < 0.05,]
write.csv(eqtl2gwas,"E:/多组学MR肺癌/结果/ICD基因集/smr-eqtl/NSCLC/NSCLC_IEUeqtlgen汇总.csv")

eqtl2gwas <- data.table::fread("E:/多组学MR肺癌/结果/ICD基因集/smr-eqtl/eqtlgen_NSCLC_UKB.smr",data.table = F)
eqtl2gwas$p_SMR_fdr <- p.adjust(eqtl2gwas$p_SMR,method = "fdr")
# eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_HEIDI > 0.01,]
# eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_SMR_fdr < 0.05,]
write.csv(eqtl2gwas,"E:/多组学MR肺癌/结果/ICD基因集/smr-eqtl/NSCLC/NSCLC_UKBeqtlgen汇总.csv")

eqtl2gwas <- data.table::fread("E:/多组学MR肺癌/结果/ICD基因集/smr-eqtl/eqtlgen_LUAD.smr",data.table = F)
eqtl2gwas$p_SMR_fdr <- p.adjust(eqtl2gwas$p_SMR,method = "fdr")
# eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_HEIDI > 0.01,]
# eqtl2gwas <- eqtl2gwas[eqtl2gwas$p_SMR_fdr < 0.05,]
write.csv(eqtl2gwas,"E:/多组学MR肺癌/结果/ICD基因集/smr-eqtl/NSCLC/LUAD_IEUeqtlgen汇总.csv")





### coloc ###

ts8 <- read.table(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/ER+_Breast_cancer/overlap_gene24.txt")
ts8_sig_probe <- Ensembl_GRCh38[Ensembl_GRCh38$gene_name %in% ts8$V1,]

coloc_batch_smr(probes = ts8_sig_probe$gene_id,
                smr_exe_path = "E:/药物靶点肺癌/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
                qtl_path = "E:/药物靶点肺癌/暴露数据/eqtl/cis-eQTL-SMR_20191212/",
                cis_wind = 2000,
                prepare_method = 2,
                type1 ="quant",
                SS1 =31684,
                NC1=0,
                coloc_plot_pph4 = 0.5,
                gwas_path2 = "F:/post-GWAS/data_gwas_list/data_prepare/coloc/ER+_Breast_cancer_BCACIEU_hg38.txt",
                type2 = "cc",
                SS2 = 175475,# 结局数据样本总数
                NC2 = 69501,# 结局数量病例数
                bfile = "E:/药物靶点肺癌/配套/EUR",
                out_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/共定位分析/")

## 结果汇总整理
coloc_files <- dir("E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/共定位分析/",pattern = "coloc_abf_summary.csv",full.names = T,recursive = T)
coloc_summary <- data.frame()
for (i in coloc_files) {
  coloc_summary_temp <- read.csv(i)
  coloc_summary <- rbind(coloc_summary,coloc_summary_temp)
}
write.csv(coloc_summary,file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost,fusion,magma三者交集结果/乳腺癌/共定位分析/ER+_Breast_cancer_BCACIEU_hg38_coloc_summary.csv",row.names = F)







### smr-plot ###
ensem_id <- ENSE2Gene[ENSE2Gene$Gene == "PDIA3", ]
ensem_id

smr_data <- smr_plot_prepare(smr_exe_path = "E:/药物靶点肺癌/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
                             bfile  ="E:/药物靶点肺癌/配套/EUR" ,
                             qtls_path = "E:/药物靶点肺癌/暴露数据/eqtl/cis-eQTL-SMR_20191212/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense",
                             gwas_path = "E:/多组学MR肺癌/IEU肺癌数据/LUSC/LUSC.ma",
                             queried_probe = "ENSG00000167004",
                             smr_cis_wind = 2000,#这必须和我们smr分析的原始的一致
                             glist_path ="E:/药物靶点肺癌/配套2/smr-plot-data/glist-hg19",
                             out_path = "E:/多组学MR肺癌/结果/ICD基因集/smr-eqtl/smr散点图和locus图/",
                             out_prefix = "PDIA3-LUSC")

smr_plot_effect(
  # smr_plot_prepare 处理之后的数据
  smr_data = smr_data,
  trait_name = "LUSC",
  pdf = T,
  out_path = "E:/多组学MR肺癌/结果/ICD基因集/smr-eqtl/smr散点图和locus图/",
  out_prefix = "PDIA3-LUSC"
)

smr_plot_locus(
  # smr_plot_prepare 处理之后的数据
  smr_data = smr_data,
  pdf = T,
  smr_thresh = 0.001,max_anno_probe = 3,
  smr_thresh_plot = F,
  probeNEARBY = "ENSG00000167004",
  out_path = "E:/多组学MR肺癌/结果/ICD基因集/smr-eqtl/smr散点图和locus图/",
  out_prefix = "PDIA3-LUSC",
  anno_selfdef = F
)