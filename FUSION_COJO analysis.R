### FUSION分析 ###
library(data.table)
library(DrugTargetMR)
library(tidyverse)

# Breast_cancer ----------------------------------------------------------------
b <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost/能用的结果/utmost_Breast_cancer_BCACIEU_hg38_joint_GBJ_filter_adj_p.csv")

aaa_dat <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/GTExV8_Breast_cancer_BCACIEU_hg38/GTExv8.ALL.summary.gene_p_fdr.csv",
           data.table = F)
aaa_dat <- aaa_dat %>% select(-V1)
aaa_dat$ENSG <- substr(aaa_dat$ID,1,regexpr("\\.",aaa_dat$ID)-1)

# 查看aaa_dat的ID是否是gene symbol。如果是gene_id，需要将ID转换成gene symbol
View(aaa_dat)

# Ensembl_GRCh38；Ensembl_GRCh37，均可。gene_id与gene_name，与不同参考基因组版本无关
aaa_dat1 <- merge(aaa_dat,Ensembl_GRCh38[,c("gene_id","gene_name")],by.x="ENSG",by.y="gene_id")

# 将ID改成gene_name
aaa_dat1$ID <- aaa_dat1$gene_name
candidate <- b$gene_name
aaa_dat1 <- aaa_dat1[aaa_dat1$ID %in% candidate,]
aaa_dat1 <- aaa_dat1[aaa_dat1$TWAS.P.fdr < 0.05,]
aaa_dat1[aaa_dat1$gene_name=="",] <- NA
aaa_dat1 <- na.omit(aaa_dat1)

# 染色体多分开跑
aaa_dat1_chr2_4 <- aaa_dat1 %>% filter(CHR==2|CHR==4)
aaa_dat1_chr18 <- aaa_dat1 %>% filter(CHR==18)
aaa_dat1_chr22 <- aaa_dat1 %>% filter(CHR==22)

file_genesymbol <- "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo_tmp/GTExv8.ALL.summary.gene_p_fdr1.csv"
data.table::fwrite(aaa_dat1_chr2_4,file = file_genesymbol,quote = F,sep = "\t",row.names = F)
fusion_cojo(input = file_genesymbol,
            out_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo/Breast_cancer_BCACIEU_hg38/chr2_4/", # 要改的地方
            out_prefix = "candidate",
            sumstats = "F:/post-GWAS/tmp2_fusion1/Breast_cancer_BCACIEU_hg38.sumstats", # 要改的地方
            ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
            ref_ld_chr_num =  c(2,4), # 要改的地方
            locus_win = 5e+05,
            p_adj_method = "fdr",
            plot = T,
            plot_legend = "joint",
            glist_hg19 = "F:/post-GWAS/data_fusion/glist-hg19")


aaa_dat1_chr18 <- aaa_dat1_chr18[-3,]
file_genesymbol <- "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo_tmp/GTExv8.ALL.summary.gene_p_fdr2.csv"
data.table::fwrite(aaa_dat1_chr18,file = file_genesymbol,quote = F,sep = "\t",row.names = F)
fusion_cojo(input = file_genesymbol,
            out_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo/Breast_cancer_BCACIEU_hg38/chr18/", # 要改的地方
            out_prefix = "candidate",
            sumstats = "F:/post-GWAS/tmp2_fusion1/Breast_cancer_BCACIEU_hg38.sumstats", # 要改的地方
            ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
            ref_ld_chr_num =  18, # 要改的地方
            locus_win = 5e+05,
            p_adj_method = "fdr",
            plot = T,
            plot_legend = "joint",
            glist_hg19 = "F:/post-GWAS/data_fusion/glist-hg19")


file_genesymbol <- "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo_tmp/GTExv8.ALL.summary.gene_p_fdr3.csv"
data.table::fwrite(aaa_dat1_chr22,file = file_genesymbol,quote = F,sep = "\t",row.names = F)
fusion_cojo(input = file_genesymbol,
            out_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo/Breast_cancer_BCACIEU_hg38/chr22/", # 要改的地方
            out_prefix = "candidate",
            sumstats = "F:/post-GWAS/tmp2_fusion1/Breast_cancer_BCACIEU_hg38.sumstats", # 要改的地方
            ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
            ref_ld_chr_num =  22, # 要改的地方
            locus_win = 5e+05,
            p_adj_method = "fdr",
            plot = T,
            plot_legend = "joint",
            glist_hg19 = "F:/post-GWAS/data_fusion/glist-hg19")

# ER-_Breast_cancer ----------------------------------------------------------------
b <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost/能用的结果/utmost_ER-_Breast_cancer_BCACIEU_hg38_joint_GBJ_filter_adj_p.csv")

aaa_dat <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/GTExV8_ER-_Breast_cancer_BCACIEU_hg38/GTExv8.ALL.summary.gene_p_fdr.csv",
                 data.table = F)
aaa_dat <- aaa_dat %>% select(-V1)
aaa_dat$ENSG <- substr(aaa_dat$ID,1,regexpr("\\.",aaa_dat$ID)-1)

# 查看aaa_dat的ID是否是gene symbol。如果是gene_id，需要将ID转换成gene symbol
View(aaa_dat)

# Ensembl_GRCh38；Ensembl_GRCh37，均可。gene_id与gene_name，与不同参考基因组版本无关
aaa_dat1 <- merge(aaa_dat,Ensembl_GRCh38[,c("gene_id","gene_name")],by.x="ENSG",by.y="gene_id")

# 将ID改成gene_name
aaa_dat1$ID <- aaa_dat1$gene_name
candidate <- b$gene_name
aaa_dat1 <- aaa_dat1[aaa_dat1$ID %in% candidate,]
aaa_dat1 <- aaa_dat1[aaa_dat1$TWAS.P.fdr < 0.05,]
aaa_dat1[aaa_dat1$gene_name=="",] <- NA
aaa_dat1 <- na.omit(aaa_dat1)

file_genesymbol <- "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo_tmp/GTExv8.ALL.summary.gene_p_fdr1.csv"
data.table::fwrite(aaa_dat1,file = file_genesymbol,quote = F,sep = "\t",row.names = F)
fusion_cojo(input = file_genesymbol,
            out_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo/ER-_Breast_cancer_BCACIEU_hg38/", # 要改的地方
            out_prefix = "candidate",
            sumstats = "F:/post-GWAS/tmp2_fusion1/ER-_Breast_cancer_BCACIEU_hg38.sumstats", # 要改的地方
            ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
            ref_ld_chr_num =  c(2,4), # 要改的地方
            locus_win = 5e+05,
            p_adj_method = "fdr",
            plot = T,
            plot_legend = "joint",
            glist_hg19 = "F:/post-GWAS/data_fusion/glist-hg19")

# ER+_Breast_cancer ----------------------------------------------------------------
b <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost/能用的结果/utmost_ER+_Breast_cancer_BCACIEU_hg38_joint_GBJ_filter_adj_p.csv")

aaa_dat <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/GTExV8_ER+_Breast_cancer_BCACIEU_hg38/GTExv8.ALL.summary.gene_p_fdr.csv",
                 data.table = F)
aaa_dat <- aaa_dat %>% select(-V1)
aaa_dat$ENSG <- substr(aaa_dat$ID,1,regexpr("\\.",aaa_dat$ID)-1)

# 查看aaa_dat的ID是否是gene symbol。如果是gene_id，需要将ID转换成gene symbol
View(aaa_dat)

# Ensembl_GRCh38；Ensembl_GRCh37，均可。gene_id与gene_name，与不同参考基因组版本无关
aaa_dat1 <- merge(aaa_dat,Ensembl_GRCh38[,c("gene_id","gene_name")],by.x="ENSG",by.y="gene_id")

# 将ID改成gene_name
aaa_dat1$ID <- aaa_dat1$gene_name
candidate <- b$gene_name
aaa_dat1 <- aaa_dat1[aaa_dat1$ID %in% candidate,]
aaa_dat1 <- aaa_dat1[aaa_dat1$TWAS.P.fdr < 0.05,]
aaa_dat1[aaa_dat1$gene_name=="",] <- NA
aaa_dat1 <- na.omit(aaa_dat1)

aaa_dat1_chr2_4_22 <- aaa_dat1 %>% filter(CHR==2|CHR==4|CHR==22)
aaa_dat1_chr18 <- aaa_dat1 %>% filter(CHR==18)

file_genesymbol <- "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo_tmp/GTExv8.ALL.summary.gene_p_fdr1.csv"
data.table::fwrite(aaa_dat1_chr2_4_22,file = file_genesymbol,quote = F,sep = "\t",row.names = F)
fusion_cojo(input = file_genesymbol,
            out_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo/ER+_Breast_cancer_BCACIEU_hg38/chr2_4_22/", # 要改的地方
            out_prefix = "candidate",
            sumstats = "F:/post-GWAS/tmp2_fusion1/ER+_Breast_cancer_BCACIEU_hg38.sumstats", # 要改的地方
            ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
            ref_ld_chr_num =  c(2,4,22), # 要改的地方
            locus_win = 5e+05,
            p_adj_method = "fdr",
            plot = T,
            plot_legend = "joint",
            glist_hg19 = "F:/post-GWAS/data_fusion/glist-hg19")


aaa_dat1_chr18 <- aaa_dat1_chr18[-2,]
file_genesymbol <- "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo_tmp/GTExv8.ALL.summary.gene_p_fdr2.csv"
data.table::fwrite(aaa_dat1_chr18,file = file_genesymbol,quote = F,sep = "\t",row.names = F)
fusion_cojo(input = file_genesymbol,
            out_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo/ER+_Breast_cancer_BCACIEU_hg38/chr18/", # 要改的地方
            out_prefix = "candidate",
            sumstats = "F:/post-GWAS/tmp2_fusion1/ER+_Breast_cancer_BCACIEU_hg38.sumstats", # 要改的地方
            ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
            ref_ld_chr_num =  18, # 要改的地方
            locus_win = 5e+05,
            p_adj_method = "fdr",
            plot = T,
            plot_legend = "joint",
            glist_hg19 = "F:/post-GWAS/data_fusion/glist-hg19")

# LUSC ----------------------------------------------------------------
b <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost/能用的结果/utmost_LUSC_TRICL_hg38_joint_GBJ_filter_adj_p.csv")

aaa_dat <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/GTExV8_LUSC_TRICL_hg38/GTExv8.ALL.summary.gene_p_fdr.csv",
                 data.table = F)
aaa_dat <- aaa_dat %>% select(-V1)
aaa_dat$ENSG <- substr(aaa_dat$ID,1,regexpr("\\.",aaa_dat$ID)-1)

# Ensembl_GRCh38；Ensembl_GRCh37，均可。gene_id与gene_name，与不同参考基因组版本无关
aaa_dat1 <- merge(aaa_dat,Ensembl_GRCh38[,c("gene_id","gene_name")],by.x="ENSG",by.y="gene_id")

# 将ID改成gene_name
aaa_dat1$ID <- aaa_dat1$gene_name
candidate <- b$gene_name
aaa_dat1 <- aaa_dat1[aaa_dat1$ID %in% candidate,]
aaa_dat1 <- aaa_dat1[aaa_dat1$TWAS.P.fdr < 0.05,]
aaa_dat1[aaa_dat1$gene_name=="",] <- NA
aaa_dat1 <- na.omit(aaa_dat1)

file_genesymbol <- "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo_tmp/GTExv8.ALL.summary.gene_p_fdr1.csv"
data.table::fwrite(aaa_dat1,file = file_genesymbol,quote = F,sep = "\t",row.names = F)
fusion_cojo(input = file_genesymbol,
            out_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo/LUSC_TRICL_hg38/", # 要改的地方
            out_prefix = "candidate",
            sumstats = "F:/post-GWAS/tmp2_fusion1/LUSC_TRICL_hg38.sumstats", # 要改的地方
            ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
            ref_ld_chr_num =  2, # 要改的地方
            locus_win = 5e+05,
            p_adj_method = "fdr",
            plot = T,
            plot_legend = "joint",
            glist_hg19 = "F:/post-GWAS/data_fusion/glist-hg19")

# Melanoma_FinnGenR10 ----------------------------------------------------------------
b <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost/能用的结果/utmost_Melanoma_FinnGenR10_joint_GBJ_filter_adj_p.csv")

aaa_dat <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/GTExV8_Melanoma_FinnGenR10/GTExv8.ALL.summary.gene_p_fdr.csv",
                 data.table = F)
aaa_dat <- aaa_dat %>% select(-V1)
aaa_dat$ENSG <- substr(aaa_dat$ID,1,regexpr("\\.",aaa_dat$ID)-1)

# Ensembl_GRCh38；Ensembl_GRCh37，均可。gene_id与gene_name，与不同参考基因组版本无关
aaa_dat1 <- merge(aaa_dat,Ensembl_GRCh38[,c("gene_id","gene_name")],by.x="ENSG",by.y="gene_id")

# 将ID改成gene_name
aaa_dat1$ID <- aaa_dat1$gene_name
candidate <- b$gene_name
aaa_dat1 <- aaa_dat1[aaa_dat1$ID %in% candidate,]
aaa_dat1 <- aaa_dat1[aaa_dat1$TWAS.P.fdr < 0.05,]
aaa_dat1[aaa_dat1$gene_name=="",] <- NA
aaa_dat1 <- na.omit(aaa_dat1)

file_genesymbol <- "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo_tmp/GTExv8.ALL.summary.gene_p_fdr1.csv"
data.table::fwrite(aaa_dat1,file = file_genesymbol,quote = F,sep = "\t",row.names = F)
fusion_cojo(input = file_genesymbol,
            out_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo/Melanoma_FinnGenR10/", # 要改的地方
            out_prefix = "candidate",
            sumstats = "F:/post-GWAS/tmp2_fusion1/Melanoma_FinnGenR10.sumstats", # 要改的地方
            ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
            ref_ld_chr_num =  22, # 要改的地方
            locus_win = 5e+05,
            p_adj_method = "fdr",
            plot = T,
            plot_legend = "joint",
            glist_hg19 = "F:/post-GWAS/data_fusion/glist-hg19")

# Melanoma_UKB_hg38 ----------------------------------------------------------------
b <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost/能用的结果/utmost_Melanoma_UKB_hg38_joint_GBJ_filter_adj_p.csv")

aaa_dat <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/GTExV8_Melanoma_UKB_hg38/GTExv8.ALL.summary.gene_p_fdr.csv",
                 data.table = F)
aaa_dat <- aaa_dat %>% select(-V1)
aaa_dat$ENSG <- substr(aaa_dat$ID,1,regexpr("\\.",aaa_dat$ID)-1)

# Ensembl_GRCh38；Ensembl_GRCh37，均可。gene_id与gene_name，与不同参考基因组版本无关
aaa_dat1 <- merge(aaa_dat,Ensembl_GRCh38[,c("gene_id","gene_name")],by.x="ENSG",by.y="gene_id")

# 将ID改成gene_name
aaa_dat1$ID <- aaa_dat1$gene_name
candidate <- b$gene_name
aaa_dat1 <- aaa_dat1[aaa_dat1$ID %in% candidate,]
aaa_dat1 <- aaa_dat1[aaa_dat1$TWAS.P.fdr < 0.05,]
aaa_dat1[aaa_dat1$gene_name=="",] <- NA
aaa_dat1 <- na.omit(aaa_dat1)

file_genesymbol <- "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo_tmp/GTExv8.ALL.summary.gene_p_fdr1.csv"
data.table::fwrite(aaa_dat1,file = file_genesymbol,quote = F,sep = "\t",row.names = F)
fusion_cojo(input = file_genesymbol,
            out_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo/Melanoma_UKB_hg38/", # 要改的地方
            out_prefix = "candidate",
            sumstats = "F:/post-GWAS/tmp2_fusion1/Melanoma_UKB_hg38.sumstats", # 要改的地方
            ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
            ref_ld_chr_num =  2, # 要改的地方
            locus_win = 5e+05,
            p_adj_method = "fdr",
            plot = T,
            plot_legend = "joint",
            glist_hg19 = "F:/post-GWAS/data_fusion/glist-hg19")

# Prostate_cancer ----------------------------------------------------------------
b <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost/能用的结果/utmost_Prostate_cancer_FinnGenR10_joint_GBJ_filter_adj_p.csv")

aaa_dat <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/GTExV8_Prostate_cancer_FinnGenR10/GTExv8.ALL.summary.gene_p_fdr.csv",
                 data.table = F)
aaa_dat <- aaa_dat %>% select(-V1)
aaa_dat$ENSG <- substr(aaa_dat$ID,1,regexpr("\\.",aaa_dat$ID)-1)

# Ensembl_GRCh38；Ensembl_GRCh37，均可。gene_id与gene_name，与不同参考基因组版本无关
aaa_dat1 <- merge(aaa_dat,Ensembl_GRCh38[,c("gene_id","gene_name")],by.x="ENSG",by.y="gene_id")

# 将ID改成gene_name
aaa_dat1$ID <- aaa_dat1$gene_name
candidate <- b$gene_name
aaa_dat1 <- aaa_dat1[aaa_dat1$ID %in% candidate,]
aaa_dat1 <- aaa_dat1[aaa_dat1$TWAS.P.fdr < 0.05,]
aaa_dat1[aaa_dat1$gene_name=="",] <- NA
aaa_dat1 <- na.omit(aaa_dat1)

file_genesymbol <- "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo_tmp/GTExv8.ALL.summary.gene_p_fdr1.csv"
data.table::fwrite(aaa_dat1,file = file_genesymbol,quote = F,sep = "\t",row.names = F)
fusion_cojo(input = file_genesymbol,
            out_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo/Prostate_cancer_FinnGenR10/", # 要改的地方
            out_prefix = "candidate",
            sumstats = "F:/post-GWAS/tmp2_fusion1/Prostate_cancer_FinnGenR10.sumstats", # 要改的地方
            ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
            ref_ld_chr_num =  c(2,4,18), # 要改的地方
            locus_win = 5e+05,
            p_adj_method = "fdr",
            plot = T,
            plot_legend = "joint",
            glist_hg19 = "F:/post-GWAS/data_fusion/glist-hg19")

# Thyroid_cancer ----------------------------------------------------------------
b <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/utmost/能用的结果/utmost_Thyroid_cancer_GWASCatalog_hg38_joint_GBJ_filter_adj_p.csv")

aaa_dat <- fread(file = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/GTExV8_Thyroid_cancer_GWASCatalog_hg38/GTExv8.ALL.summary.gene_p_fdr.csv",
                 data.table = F)
aaa_dat <- aaa_dat %>% select(-V1)
aaa_dat$ENSG <- substr(aaa_dat$ID,1,regexpr("\\.",aaa_dat$ID)-1)

# Ensembl_GRCh38；Ensembl_GRCh37，均可。gene_id与gene_name，与不同参考基因组版本无关
aaa_dat1 <- merge(aaa_dat,Ensembl_GRCh38[,c("gene_id","gene_name")],by.x="ENSG",by.y="gene_id")

# 将ID改成gene_name
aaa_dat1$ID <- aaa_dat1$gene_name
candidate <- b$gene_name
aaa_dat1 <- aaa_dat1[aaa_dat1$ID %in% candidate,]
aaa_dat1 <- aaa_dat1[aaa_dat1$TWAS.P.fdr < 0.05,]
aaa_dat1[aaa_dat1$gene_name=="",] <- NA
aaa_dat1 <- na.omit(aaa_dat1)

file_genesymbol <- "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo_tmp/GTExv8.ALL.summary.gene_p_fdr1.csv"
data.table::fwrite(aaa_dat1,file = file_genesymbol,quote = F,sep = "\t",row.names = F)
fusion_cojo(input = file_genesymbol,
            out_path = "E:/post_gwas_cancer/results/utmost和fusion初筛/fusion/fusion_cojo/Thyroid_cancer_GWASCatalog_hg38/", # 要改的地方
            out_prefix = "candidate",
            sumstats = "F:/post-GWAS/tmp2_fusion1/Thyroid_cancer_GWASCatalog_hg38.sumstats", # 要改的地方
            ref_ld_chr_prefix = "F:/post-GWAS/data_fusion/LDREF_hg38/1000G.EUR.",
            ref_ld_chr_num =  c(2,4,18), # 要改的地方
            locus_win = 5e+05,
            p_adj_method = "fdr",
            plot = T,
            plot_legend = "joint",
            glist_hg19 = "F:/post-GWAS/data_fusion/glist-hg19")