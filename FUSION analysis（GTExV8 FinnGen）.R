### FUSION分析 ###
library(DrugTargetMR)

# step1.结局文件准备 -----------------------------------------------------------
## 说明1： 将GWAS标准文件，再次处理为.sumstats文件。
## 说明2： fusion_prepare_sumstat()函数返回的sumstat_file，为输出文件的路径
dir <- list.files("F:/tmp/")# GWAS所在文件夹
data_dir <- "F:/tmp/"
n <- length(dir)
out_path <- "F:/tmp1/"

for (i in 1:n) {
  file_path <- paste0(data_dir,"",dir[i])
  sumstat_file <- fusion_prepare_sumstat(file_path = file_path,
                                         out_path = out_path)
}

# step2.权重文件准备 -----------------------------------------------------------
## 权重文件准备 
## 获取weights和weights_dir对应的参数，制作成data.frame
dir_untar <- "F:/post-GWAS/data_fusion/权重文件/Mancuso-GTEx多组织基因权重/GTExv8.ALL2"
dir_weight <- list.files(dir_untar,recursive = FALSE, full.names = TRUE)
file_weights <-  unlist(lapply(dir_weight, FUN = function(x){ dir(x,pattern = ".pos$",recursive = F,full.names = T)}))
file_weights <- file_weights[!grepl(".nofilter.pos",file_weights)]
df_weights <- data.frame(file_weights = file_weights,weights_dir = dirname(file_weights))
df_weights$tissues <- basename(df_weights$file_weights)

## 查看
View(df_weights)

## 提前导入用于fusion分析的sumstat_data数据，类型为data.frame，可减少循环过程中的运算量
sumstat_data <- data.table::fread("F:/tmp1/finngen_R10_C3_BREAST_ERNEG_EXALLC.sumstats",data.table=F)
## 使用循环，挨个组织执行FUSION分析
for (i in 1:nrow(df_weights)) {
  fusion_assoc(sumstat_data = sumstat_data,
               weights = df_weights$file_weights[i],
               weights_dir =df_weights$weights_dir[i],
               ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
               ref_ld_chr_num = 1:22,
               out_path ="E:/post_gwas_cancer/tmp1/ER-_Breast_Cancer_FinnGen/",
               out_prefix = df_weights$tissues[i])
}

sumstat_data <- data.table::fread("F:/tmp1/finngen_R10_C3_BREAST_ERPLUS_EXALLC.sumstats",data.table=F)
## 使用循环，挨个组织执行FUSION分析
for (i in 1:nrow(df_weights)) {
  fusion_assoc(sumstat_data = sumstat_data,
               weights = df_weights$file_weights[i],
               weights_dir =df_weights$weights_dir[i],
               ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
               ref_ld_chr_num = 1:22,
               out_path ="E:/post_gwas_cancer/tmp1/ER+_Breast_Cancer_FinnGen/",
               out_prefix = df_weights$tissues[i])
}

sumstat_data <- data.table::fread("F:/tmp1/finngen_R10_C3_BREAST_EXALLC.sumstats",data.table=F)
## 使用循环，挨个组织执行FUSION分析
for (i in 1:nrow(df_weights)) {
  fusion_assoc(sumstat_data = sumstat_data,
               weights = df_weights$file_weights[i],
               weights_dir =df_weights$weights_dir[i],
               ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
               ref_ld_chr_num = 1:22,
               out_path ="E:/post_gwas_cancer/tmp1/Breast_Cancer_FinnGen/",
               out_prefix = df_weights$tissues[i])
}

### 结果整理
res_fusion_files <- dir("E:/post_gwas_cancer/tmp1/ER+_Breast_Cancer_FinnGen/",pattern = ".csv",recursive = F,full.names = T)
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
write.csv(res_fusion,file = "E:/post_gwas_cancer/tmp1/ER+_Breast_Cancer_FinnGen/GTExv8.ALL.summary.gene_p_fdr.csv")
write.csv(res_fusion_f,file = "E:/post_gwas_cancer/tmp1/ER+_Breast_Cancer_FinnGen/GTExv8.ALL.summary.gene_p_fdr_filtered.csv")

# 统计每个组织中的基因数量
sum1 <- as.data.frame(table(res_fusion_f$PANEL))

# 统计一共发现多少个基因
sum2 <- length(unique(res_fusion_f$gene_id))

