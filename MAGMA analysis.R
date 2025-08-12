### MAGMA分析 ###
library(DrugTargetMR)

## 坐标轴转换：将hg38版坐标轴，转换为hg19版本。
dir <- list.files("F:/post-GWAS/data_gwas_list/data_prepare/coloc/")# GWAS所在文件夹
data_dir <- "F:/post-GWAS/data_gwas_list/data_prepare/coloc/"
n <- length(dir)

for (i in 1:n) {
  file_path = paste0(data_dir,dir[i])
  transform_hg38ToHg19(file_path = file_path,
                     out_path = "F:/post-GWAS/data_gwas_list/data_prepare/magma_hg37/")
}

## 执行magma分析
dir <- list.files("F:/post-GWAS/tmp3_magma/")# GWAS所在文件夹
data_dir <- "F:/post-GWAS/tmp3_magma/"
n <- length(dir)

for (i in 1:n) {
  #路径设置
  gwas_path = paste0(data_dir,dir[i])
  prefix = sub("\\.txt$", "", dir[i])
  dir.create(paste0("E:/post_gwas_cancer/results/utmost和fusion初筛/magma/",prefix))
  out_path = paste0("E:/post_gwas_cancer/results/utmost和fusion初筛/magma/",prefix,"/")
  readcsv_path = paste0(out_path,"magma_step2_snp_loc.genes.out.mapped.csv")
  writecsv_path = paste0(out_path,"magma_snp2genes_filter_adj_p.csv")

  # 开始分析
  postgwas_magma(magma_exe = "F:/post-GWAS/data_magma/MAGMA数据包/magma软件/magma.exe",
               gwas_path = gwas_path,
               bfile = "F:/post-GWAS/data_magma/MAGMA数据包/LD参考文件/g1000_eur",
               geneloc_file = "F:/post-GWAS/data_magma/MAGMA数据包/MAGMA_gene_boundary_files/ENSGv110.coding.genes.GRCh37.txt",  # 基因注释为hg37版本
               geneset_anno = T,
               geneset_file = "F:/post-GWAS/data_magma/MAGMA数据包/Gene_set_files/MSigDB_20231Hs_MAGMA.txt",
               genecovar_file = "F:/post-GWAS/data_magma/MAGMA数据包/Gene_expression_files/gtex_v8_ts_avg_log2TPM_54tissues.txt",
               out_path = out_path,
               out_prefix = "magma")

  # 结果整理
  magma_genes <- read.csv(readcsv_path)
  magma_genes$p_fdr <- p.adjust(magma_genes$P,method = "fdr")
  magma_genes_f <- magma_genes[magma_genes$p_fdr < 0.05,]
  write.csv(magma_genes_f,file = writecsv_path,row.names = F)

}
