#### MR药物靶点+共定位 ####
library(DrugTargetMR)
library(openxlsx)
library(data.table)
library(clusterProfiler)
library(org.Mm.eg.db)
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(openxlsx)
library(qvalue)

#### 乳腺癌 ####
################# decode ##############
#step 1. read exposure data 
a <- read.csv(file = "F:/post-GWAS/data_mr/decode.csv",sep = ",",header = T)
colnames(a)
exposure_dat <- read_exposure_data("F:/post-GWAS/data_mr/decode.csv", 
                                   clump = F, sep = ",",
                                   phenotype_col = "gene.exposure",
                                   snp = "SNP", beta_col = "beta.exposure", 
                                   se_col = "se.exposure", 
                                   effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure", 
                                   pval_col = "pval.exposure", samplesize_col = "samplesize.exposure", 
                                   eaf_col = "eaf.exposure",
                                   chr_col = "chr.exposure",
                                   pos_col = "pos.exposure") 
exposure_dat <- exposure_dat %>% filter(mr_keep.exposure==TRUE)

#step 2. read outcome data 
a <- fread(file = "F:/post-GWAS/data_gwas_list/data_prepare/mr_with_phenotype/Breast_cancer_BCACIEU.txt")
head(a)
outcome_dat1<-read_outcome_data(filename="F:/post-GWAS/data_gwas_list/data_prepare/mr_with_phenotype/Breast_cancer_BCACIEU.txt", 
                                snps = exposure_dat$SNP,
                                sep = "\t",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",phenotype_col = "phenotype_col",
                                # ncase_col = "ncase_col",ncontrol_col = "ncontrol",
                                samplesize_col = "samplesize")
outcome_dat2<-read_outcome_data(filename="F:/post-GWAS/data_gwas_list/data_prepare/mr_with_phenotype/ER-_Breast_cancer_BCACIEU.txt", 
                                snps = exposure_dat$SNP,
                                sep = "\t",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",phenotype_col = "phenotype_col",
                                # ncase_col = "ncase_col",ncontrol_col = "ncontrol",
                                samplesize_col = "samplesize")
outcome_dat3<-read_outcome_data(filename="F:/post-GWAS/data_gwas_list/data_prepare/mr_with_phenotype/ER+_Breast_cancer_BCACIEU.txt", 
                                snps = exposure_dat$SNP,
                                sep = "\t",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",phenotype_col = "phenotype_col",
                                # ncase_col = "ncase_col",ncontrol_col = "ncontrol",
                                samplesize_col = "samplesize")


#step 4. harmonise 
outcome_dat<-rbind(outcome_dat1,outcome_dat2,outcome_dat3)
dat <- harmonise_data(exposure_dat, outcome_dat)
table(dat$exposure)

#step 5. caculate F-stat for each SNP 
dat$EAF2 <- (1 - dat$eaf.exposure) 
dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2) 
PVEfx <- function(BETA, MAF, SE, N){ 
  pve <- (2*(BETA^2)*MAF*(1 - MAF))/((2*(BETA^2)*MAF*(1 - MAF)) + ((SE^2)*2*N*MAF*(1 - MAF))) 
  return(pve) } 
dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N = dat$samplesize.exposure) 
dat$FSTAT <- ((dat$samplesize.exposure - 1 - 1)/1)*(dat$PVE/(1 - dat$PVE))
dat <- subset(dat,FSTAT>10)
dat <- subset(dat,MAF>0.01)

#step 6. MR analysis using Inverse variance weighted method 
res_ivw <- generate_odds_ratios(mr_res=mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_wald_ratio",
                                                              "mr_weighted_median","mr_simple_mode",
                                                              "mr_weighted_mode","mr_wald_ratio")
))
write.csv(res_ivw,file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/ivw_result.csv")
openxlsx::write.xlsx(res_ivw, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/ivw_result.xlsx",
                     sheetName = "Sheet1", colNames = TRUE)

#step 7. p-adj 
filtered_data <- res_ivw %>% filter(method==c("Inverse variance weighted","Wald ratio"))
name<-filtered_data$outcome[!duplicated(filtered_data$outcome)]
number_outcome <- filtered_data %>% group_by(outcome) %>% summarize(count = n()) 
number_exposure <- filtered_data %>% group_by(exposure) %>% summarize(count = n()) 
filtered_data1<-filtered_data
filtered_data1<-filtered_data%>%arrange(pval)
filtered_data1$adjusted_p <- NA
for (i in 1:length(name)) {
  b <- filtered_data1 %>% filter(outcome == name[i])
  a <- b[order(b$pval), ]
  plist <- a$pval
  p <- plist
  sorted_p_values <- sort(p)
  adjusted_p_values <- p.adjust(sorted_p_values, method = "fdr")
  filtered_data1$adjusted_p[filtered_data1$outcome == name[i]] <- adjusted_p_values
}
openxlsx::write.xlsx(filtered_data1, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/ivw_p_adj.xlsx",
                     sheetName = "Sheet1", colNames = TRUE) 


#step 8. Sensitive:pleiotropy analysis and heterogeneity test, heterogeneity (Inverse variance weighted), stiger-test, LOO 
# pval<0.05
dat_sen1<-res_ivw %>% filter(pval<0.05)
dat_sen11<-dat[which(dat$exposure%in%dat_sen1$exposure), ]
dat_sen<-dat_sen11
write.csv(dat_sen,file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/dat_sen.csv")

# heterogeneity test
mr_results_het <- mr_heterogeneity(dat_sen)
b1<-subset(res_ivw,pval<0.05)
b<-subset(b1,method==c("Inverse variance weighted","Wald ratio"))
mr_results_het_use<- merge(mr_results_het, b[, c("outcome", "exposure")], by = c("outcome", "exposure"))
openxlsx::write.xlsx(mr_results_het_use, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/het.xlsx",
                     sheetName = "Sheet1", Colnames = TRUE)

# pleiotropy analysis 
mr_results_ple <- mr_pleiotropy_test(dat_sen)
c1<-subset(res_ivw,pval<0.05)
c<-subset(c1,method==c("Inverse variance weighted","Wald ratio"))
mr_results_ple_use<- merge(mr_results_ple, c[, c("outcome", "exposure")], by = c("outcome", "exposure"))
openxlsx::write.xlsx(mr_results_ple_use, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/ple.xlsx",
                     sheetName = "Sheet1", Colnames = TRUE)

# steiger方向性检测 
steiger_whole_test <- directionality_test(dat =dat )
steiger_snp_test <- steiger_filtering(dat =dat )
write.csv(steiger_whole_test,file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/steiger.csv")

# MRPRESSO 
dat_sen1 <- res_ivw %>%
  filter(pval < 0.05, method == "Inverse variance weighted")
name <- dat_sen1 %>%
  count(exposure, outcome) %>%
  filter(n != 0) %>%
  select(exposure, outcome)
n <- nrow(name)

merged_data1 <- dat %>% filter(exposure==name$exposure[1]&outcome==name$outcome[1])
for(i in 2:n){new.data <- dat %>% filter(exposure==name$exposure[i]&outcome==name$outcome[i])
merged_data1 = rbind(merged_data1,new.data)}
dat_pro_tmp <- merged_data1 %>% group_by(exposure,outcome) %>% summarize(count = n()) %>% filter(count >= 4)
n <- nrow(dat_pro_tmp)

merged_data2 <- dat %>% filter(exposure==dat_pro_tmp$exposure[1]&outcome==dat_pro_tmp$outcome[1])
for(i in 2:n){new.data <- dat %>% filter(exposure==dat_pro_tmp$exposure[i]&outcome==dat_pro_tmp$outcome[i])
merged_data2 = rbind(merged_data2,new.data)}
dat_pro_tmp2 <- merged_data2 %>% group_by(exposure,outcome) %>% summarize(count = n())
dat_pro <- merged_data2

# mr_result_pro<-MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
#           BetaExposure = "beta.exposure",
#           SdOutcome = "se.outcome",
#           SdExposure = "se.exposure",
#           OUTLIERtest = TRUE,
#           DISTORTIONtest = TRUE,
#           data = dat_pro,
#           NbDistribution = 1000,
#           SignifThreshold = 0.05)
mr_result_pro<-run_mr_presso(dat_pro, NbDistribution = 1000, SignifThreshold = 0.05)
mr_results_pro <- mr_result_pro

significant_results <- list()
# 遍历每个结果块
for (i in seq_along(mr_results_pro)) {
  # 获取当前结果块的Global Test Pvalue
  p_value <- mr_results_pro[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  # 检查Pvalue是否小于0.05
  if (!is.na(p_value) && p_value < 0.05) {
    # 将符合条件的结果和位置添加到列表中
    significant_results[[length(significant_results) + 1]] <- list(result = mr_results_pro[[i]], position = i)
  }
}

Nosignificant_results <- list()
# 遍历每个结果块
for (i in seq_along(mr_results_pro)) {
  # 获取当前结果块的Global Test Pvalue
  p_value <- mr_results_pro[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  # 检查Pvalue是否小于0.05
  if (!is.na(p_value) && p_value > 0.05) {
    # 将符合条件的结果和位置添加到列表中
    Nosignificant_results[[length(Nosignificant_results) + 1]] <- list(result = mr_results_pro[[i]], position = i)
  }
}

extracted_data_sig <- data.frame(exposure = character(), outcome = character(), stringsAsFactors = FALSE)
for (i in seq_along(significant_results)){
  position <- significant_results[[i]]$position
  exposure <- dat_pro_tmp2$exposure[position]
  outcome <- dat_pro_tmp2$outcome[position]
  extracted_data_sig <- rbind(extracted_data_sig, data.frame(exposure = exposure, outcome = outcome, stringsAsFactors = FALSE))
}

extracted_data_nosig <- data.frame(exposure = character(), outcome = character(), stringsAsFactors = FALSE)
for (i in seq_along(Nosignificant_results)){
  position <- Nosignificant_results[[i]]$position
  exposure <- dat_pro_tmp2$exposure[position]
  outcome <- dat_pro_tmp2$outcome[position]
  extracted_data_nosig <- rbind(extracted_data_nosig, data.frame(exposure = exposure, outcome = outcome, stringsAsFactors = FALSE))
}

# 创建一个空的数据框来存储最终结果
final_result_sig <- data.frame(Estimate = numeric(0), Sd = numeric(0), CI_Up = numeric(0), CI_Down = numeric(0), 
                               OR_Estimate = character(0), OR_CI = character(0), stringsAsFactors = FALSE)
# 循环遍历significant_results中的每个元素
for (i in seq_along(significant_results)) {
  # 从每个子集中提取数据
  data_i <- significant_results[[i]]$result$`Main MR results`
  # 提取Causal Estimate和Sd
  estimate_i <- data_i$`Causal Estimate`[1]
  sd_i <- data_i$Sd[1]
  # 创建数据框并计算上下置信区间，排除NA值
  result_i <- data.frame(Estimate = estimate_i, Sd = sd_i)
  result_i$`CI_Up` <- ifelse(!is.na(result_i$Estimate), result_i$Estimate + 1.96 * result_i$Sd, NA)
  result_i$`CI_Down` <- ifelse(!is.na(result_i$Estimate), result_i$Estimate - 1.96 * result_i$Sd, NA)
  # 计算OR及其95%置信区间
  result_i$OR_Estimate <- ifelse(!is.na(result_i$Estimate), round(exp(result_i$Estimate), 5), NA)
  result_i$OR_CI <- ifelse(!is.na(result_i$Estimate), paste("(", round(exp(result_i$`CI_Down`), 5), " - ", round(exp(result_i$`CI_Up`), 5), ")", sep = ""), NA)
  # 排除包含NA值的行
  result_i <- result_i[complete.cases(result_i),]
  # 将结果合并到最终数据框中
  final_result_sig <- rbind(final_result_sig, result_i)
}
Pval<-c()
for (i in seq_along(significant_results)) {
  p_value_i <- significant_results[[i]]$result$`MR-PRESSO results`$`Global Test`$Pvalue
  Pval <- c(Pval,p_value_i)
}
final_result_sig$Pval <- Pval
final_result_sig$exposure<-extracted_data_sig$exposure
final_result_sig$outcome<-extracted_data_sig$outcome
openxlsx::write.xlsx(final_result_sig, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/PRO_sig_result.xlsx",
                     sheetName = "Sheet1", Colnames = TRUE)

# 创建一个空的数据框来存储最终结果
final_result_nosig <- data.frame(Estimate = numeric(0), Sd = numeric(0), CI_Up = numeric(0), CI_Down = numeric(0), 
                                 OR_Estimate = character(0), OR_CI = character(0), stringsAsFactors = FALSE)
# 循环遍历significant_results中的每个元素
for (i in seq_along(Nosignificant_results)) {
  # 从每个子集中提取数据
  data_i <- Nosignificant_results[[i]]$result$`Main MR results`
  # 提取Causal Estimate和Sd
  estimate_i <- data_i$`Causal Estimate`[1]
  sd_i <- data_i$Sd[1]
  # 创建数据框并计算上下置信区间，排除NA值
  result_i <- data.frame(Estimate = estimate_i, Sd = sd_i)
  result_i$`CI_Up` <- ifelse(!is.na(result_i$Estimate), result_i$Estimate + 1.96 * result_i$Sd, NA)
  result_i$`CI_Down` <- ifelse(!is.na(result_i$Estimate), result_i$Estimate - 1.96 * result_i$Sd, NA)
  # 计算OR及其95%置信区间
  result_i$OR_Estimate <- ifelse(!is.na(result_i$Estimate), round(exp(result_i$Estimate), 5), NA)
  result_i$OR_CI <- ifelse(!is.na(result_i$Estimate), paste("(", round(exp(result_i$`CI_Down`), 5), " - ", round(exp(result_i$`CI_Up`), 5), ")", sep = ""), NA)
  # 排除包含NA值的行
  result_i <- result_i[complete.cases(result_i),]
  # 将结果合并到最终数据框中
  final_result_nosig <- rbind(final_result_nosig, result_i)
}
Pval<-c()
for (i in seq_along(Nosignificant_results)) {
  p_value_i <- Nosignificant_results[[i]]$result$`MR-PRESSO results`$`Global Test`$Pvalue
  Pval <- c(Pval,p_value_i)
}
final_result_nosig$Pval <- Pval
final_result_nosig$exposure<-extracted_data_nosig$exposure
final_result_nosig$outcome<-extracted_data_nosig$outcome
openxlsx::write.xlsx(final_result_nosig, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/PRO_nosig_result.xlsx",
                     sheetName = "Sheet1", Colnames = TRUE)

#step . save result 
save(exposure_dat,outcome_dat,dat,
     dat_sen, dat_sen1, dat_pro, #dat_pro1,
     file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/all_har.RData")
save(res_ivw,#mr_results_het,mr_results_ple,
     mr_result_pro,final_result_sig,final_result_nosig,
     filtered_data1,
     #dat_padj,final_result,
     file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/all_result.RData")


################# ukbppp ##############
#step 1. read exposure data 
a <- read.csv(file = "F:/post-GWAS/data_mr/ukbppp.csv",sep = ",",header = T)
colnames(a)
exposure_dat <- read_exposure_data("F:/post-GWAS/data_mr/ukbppp.csv", 
                                   clump = F, sep = ",",
                                   phenotype_col = "gene.exposure",
                                   snp = "SNP", beta_col = "beta.exposure", 
                                   se_col = "se.exposure", 
                                   effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure", 
                                   pval_col = "pval.exposure", samplesize_col = "samplesize.exposure", 
                                   eaf_col = "eaf.exposure",
                                   chr_col = "chr.exposure",
                                   pos_col = "pos.exposure") 
exposure_dat <- exposure_dat %>% filter(mr_keep.exposure==TRUE)

#step 2. read outcome data 
# a <- fread(file = "F:/post-GWAS/data_gwas_list/data_prepare/mr_with_phenotype/Breast_cancer_BCACIEU.txt")
# head(a)
outcome_dat1<-read_outcome_data(filename="F:/post-GWAS/data_gwas_list/data_prepare/mr_with_phenotype/Breast_cancer_BCACIEU.txt", 
                                snps = exposure_dat$SNP,
                                sep = "\t",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",phenotype_col = "phenotype_col",
                                # ncase_col = "ncase_col",ncontrol_col = "ncontrol",
                                samplesize_col = "samplesize")
outcome_dat2<-read_outcome_data(filename="F:/post-GWAS/data_gwas_list/data_prepare/mr_with_phenotype/ER-_Breast_cancer_BCACIEU.txt", 
                                snps = exposure_dat$SNP,
                                sep = "\t",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",phenotype_col = "phenotype_col",
                                # ncase_col = "ncase_col",ncontrol_col = "ncontrol",
                                samplesize_col = "samplesize")
outcome_dat3<-read_outcome_data(filename="F:/post-GWAS/data_gwas_list/data_prepare/mr_with_phenotype/ER+_Breast_cancer_BCACIEU.txt", 
                                snps = exposure_dat$SNP,
                                sep = "\t",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",phenotype_col = "phenotype_col",
                                # ncase_col = "ncase_col",ncontrol_col = "ncontrol",
                                samplesize_col = "samplesize")


#step 4. harmonise 
outcome_dat<-rbind(outcome_dat1,outcome_dat2,outcome_dat3)
dat <- harmonise_data(exposure_dat, outcome_dat)
table(dat$exposure)

#step 5. caculate F-stat for each SNP 
dat$EAF2 <- (1 - dat$eaf.exposure) 
dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2) 
PVEfx <- function(BETA, MAF, SE, N){ 
  pve <- (2*(BETA^2)*MAF*(1 - MAF))/((2*(BETA^2)*MAF*(1 - MAF)) + ((SE^2)*2*N*MAF*(1 - MAF))) 
  return(pve) } 
dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N = dat$samplesize.exposure) 
dat$FSTAT <- ((dat$samplesize.exposure - 1 - 1)/1)*(dat$PVE/(1 - dat$PVE))
dat <- subset(dat,FSTAT>10)
dat <- subset(dat,MAF>0.01)

#step 6. MR analysis using Inverse variance weighted method 
res_ivw <- generate_odds_ratios(mr_res=mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_wald_ratio",
                                                              "mr_weighted_median","mr_simple_mode",
                                                              "mr_weighted_mode","mr_wald_ratio")
))
write.csv(res_ivw,file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/deocde/ivw_result.csv")
openxlsx::write.xlsx(res_ivw, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/deocde/ivw_result.xlsx",
                     sheetName = "Sheet1", colNames = TRUE)

#step 7. p-adj 
filtered_data <- res_ivw %>% filter(method==c("Inverse variance weighted","Wald ratio"))
name<-filtered_data$outcome[!duplicated(filtered_data$outcome)]
number_outcome <- filtered_data %>% group_by(outcome) %>% summarize(count = n()) 
number_exposure <- filtered_data %>% group_by(exposure) %>% summarize(count = n()) 
filtered_data1<-filtered_data
filtered_data1<-filtered_data%>%arrange(pval)
filtered_data1$adjusted_p <- NA
for (i in 1:length(name)) {
  b <- filtered_data1 %>% filter(outcome == name[i])
  a <- b[order(b$pval), ]
  plist <- a$pval
  p <- plist
  sorted_p_values <- sort(p)
  adjusted_p_values <- p.adjust(sorted_p_values, method = "fdr")
  filtered_data1$adjusted_p[filtered_data1$outcome == name[i]] <- adjusted_p_values
}
openxlsx::write.xlsx(filtered_data1, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/deocde/ivw_p_adj.xlsx",
                     sheetName = "Sheet1", colNames = TRUE) 


#step 8. Sensitive:pleiotropy analysis and heterogeneity test, heterogeneity (Inverse variance weighted), stiger-test, LOO 
# pval<0.05
dat_sen1<-res_ivw %>% filter(pval<0.05)
dat_sen11<-dat[which(dat$exposure%in%dat_sen1$exposure), ]
dat_sen<-dat_sen11
write.csv(dat_sen,file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/deocde/dat_sen.csv")

# heterogeneity test
mr_results_het <- mr_heterogeneity(dat_sen)
b1<-subset(res_ivw,pval<0.05)
b<-subset(b1,method==c("Inverse variance weighted","Wald ratio"))
mr_results_het_use<- merge(mr_results_het, b[, c("outcome", "exposure")], by = c("outcome", "exposure"))
openxlsx::write.xlsx(mr_results_het_use, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/deocde/het.xlsx",
                     sheetName = "Sheet1", Colnames = TRUE)

# pleiotropy analysis 
mr_results_ple <- mr_pleiotropy_test(dat_sen)
c1<-subset(res_ivw,pval<0.05)
c<-subset(c1,method==c("Inverse variance weighted","Wald ratio"))
mr_results_ple_use<- merge(mr_results_ple, c[, c("outcome", "exposure")], by = c("outcome", "exposure"))
openxlsx::write.xlsx(mr_results_ple_use, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/deocde/ple.xlsx",
                     sheetName = "Sheet1", Colnames = TRUE)

# steiger方向性检测 
steiger_whole_test <- directionality_test(dat =dat )
steiger_snp_test <- steiger_filtering(dat =dat )
write.csv(steiger_whole_test,file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/deocde/steiger.csv")

# MRPRESSO 
dat_sen1 <- res_ivw %>%
  filter(pval < 0.05, method == "Inverse variance weighted")
name <- dat_sen1 %>%
  count(exposure, outcome) %>%
  filter(n != 0) %>%
  select(exposure, outcome)
n <- nrow(name)

merged_data1 <- dat %>% filter(exposure==name$exposure[1]&outcome==name$outcome[1])
for(i in 2:n){new.data <- dat %>% filter(exposure==name$exposure[i]&outcome==name$outcome[i])
merged_data1 = rbind(merged_data1,new.data)}
dat_pro_tmp <- merged_data1 %>% group_by(exposure,outcome) %>% summarize(count = n()) %>% filter(count >= 4)
n <- nrow(dat_pro_tmp)

merged_data2 <- dat %>% filter(exposure==dat_pro_tmp$exposure[1]&outcome==dat_pro_tmp$outcome[1])
for(i in 2:n){new.data <- dat %>% filter(exposure==dat_pro_tmp$exposure[i]&outcome==dat_pro_tmp$outcome[i])
merged_data2 = rbind(merged_data2,new.data)}
dat_pro_tmp2 <- merged_data2 %>% group_by(exposure,outcome) %>% summarize(count = n())
dat_pro <- merged_data2

# mr_result_pro<-MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
#           BetaExposure = "beta.exposure",
#           SdOutcome = "se.outcome",
#           SdExposure = "se.exposure",
#           OUTLIERtest = TRUE,
#           DISTORTIONtest = TRUE,
#           data = dat_pro,
#           NbDistribution = 1000,
#           SignifThreshold = 0.05)
mr_result_pro<-run_mr_presso(dat_pro, NbDistribution = 1000, SignifThreshold = 0.05)
mr_results_pro <- mr_result_pro

significant_results <- list()
# 遍历每个结果块
for (i in seq_along(mr_results_pro)) {
  # 获取当前结果块的Global Test Pvalue
  p_value <- mr_results_pro[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  # 检查Pvalue是否小于0.05
  if (!is.na(p_value) && p_value < 0.05) {
    # 将符合条件的结果和位置添加到列表中
    significant_results[[length(significant_results) + 1]] <- list(result = mr_results_pro[[i]], position = i)
  }
}

Nosignificant_results <- list()
# 遍历每个结果块
for (i in seq_along(mr_results_pro)) {
  # 获取当前结果块的Global Test Pvalue
  p_value <- mr_results_pro[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  # 检查Pvalue是否小于0.05
  if (!is.na(p_value) && p_value > 0.05) {
    # 将符合条件的结果和位置添加到列表中
    Nosignificant_results[[length(Nosignificant_results) + 1]] <- list(result = mr_results_pro[[i]], position = i)
  }
}

extracted_data_sig <- data.frame(exposure = character(), outcome = character(), stringsAsFactors = FALSE)
for (i in seq_along(significant_results)){
  position <- significant_results[[i]]$position
  exposure <- dat_pro_tmp2$exposure[position]
  outcome <- dat_pro_tmp2$outcome[position]
  extracted_data_sig <- rbind(extracted_data_sig, data.frame(exposure = exposure, outcome = outcome, stringsAsFactors = FALSE))
}

extracted_data_nosig <- data.frame(exposure = character(), outcome = character(), stringsAsFactors = FALSE)
for (i in seq_along(Nosignificant_results)){
  position <- Nosignificant_results[[i]]$position
  exposure <- dat_pro_tmp2$exposure[position]
  outcome <- dat_pro_tmp2$outcome[position]
  extracted_data_nosig <- rbind(extracted_data_nosig, data.frame(exposure = exposure, outcome = outcome, stringsAsFactors = FALSE))
}

# 创建一个空的数据框来存储最终结果
final_result_sig <- data.frame(Estimate = numeric(0), Sd = numeric(0), CI_Up = numeric(0), CI_Down = numeric(0), 
                               OR_Estimate = character(0), OR_CI = character(0), stringsAsFactors = FALSE)
# 循环遍历significant_results中的每个元素
for (i in seq_along(significant_results)) {
  # 从每个子集中提取数据
  data_i <- significant_results[[i]]$result$`Main MR results`
  # 提取Causal Estimate和Sd
  estimate_i <- data_i$`Causal Estimate`[1]
  sd_i <- data_i$Sd[1]
  # 创建数据框并计算上下置信区间，排除NA值
  result_i <- data.frame(Estimate = estimate_i, Sd = sd_i)
  result_i$`CI_Up` <- ifelse(!is.na(result_i$Estimate), result_i$Estimate + 1.96 * result_i$Sd, NA)
  result_i$`CI_Down` <- ifelse(!is.na(result_i$Estimate), result_i$Estimate - 1.96 * result_i$Sd, NA)
  # 计算OR及其95%置信区间
  result_i$OR_Estimate <- ifelse(!is.na(result_i$Estimate), round(exp(result_i$Estimate), 5), NA)
  result_i$OR_CI <- ifelse(!is.na(result_i$Estimate), paste("(", round(exp(result_i$`CI_Down`), 5), " - ", round(exp(result_i$`CI_Up`), 5), ")", sep = ""), NA)
  # 排除包含NA值的行
  result_i <- result_i[complete.cases(result_i),]
  # 将结果合并到最终数据框中
  final_result_sig <- rbind(final_result_sig, result_i)
}
Pval<-c()
for (i in seq_along(significant_results)) {
  p_value_i <- significant_results[[i]]$result$`MR-PRESSO results`$`Global Test`$Pvalue
  Pval <- c(Pval,p_value_i)
}
final_result_sig$Pval <- Pval
final_result_sig$exposure<-extracted_data_sig$exposure
final_result_sig$outcome<-extracted_data_sig$outcome
openxlsx::write.xlsx(final_result_sig, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/deocde/PRO_sig_result.xlsx",
                     sheetName = "Sheet1", Colnames = TRUE)

# 创建一个空的数据框来存储最终结果
final_result_nosig <- data.frame(Estimate = numeric(0), Sd = numeric(0), CI_Up = numeric(0), CI_Down = numeric(0), 
                                 OR_Estimate = character(0), OR_CI = character(0), stringsAsFactors = FALSE)
# 循环遍历significant_results中的每个元素
for (i in seq_along(Nosignificant_results)) {
  # 从每个子集中提取数据
  data_i <- Nosignificant_results[[i]]$result$`Main MR results`
  # 提取Causal Estimate和Sd
  estimate_i <- data_i$`Causal Estimate`[1]
  sd_i <- data_i$Sd[1]
  # 创建数据框并计算上下置信区间，排除NA值
  result_i <- data.frame(Estimate = estimate_i, Sd = sd_i)
  result_i$`CI_Up` <- ifelse(!is.na(result_i$Estimate), result_i$Estimate + 1.96 * result_i$Sd, NA)
  result_i$`CI_Down` <- ifelse(!is.na(result_i$Estimate), result_i$Estimate - 1.96 * result_i$Sd, NA)
  # 计算OR及其95%置信区间
  result_i$OR_Estimate <- ifelse(!is.na(result_i$Estimate), round(exp(result_i$Estimate), 5), NA)
  result_i$OR_CI <- ifelse(!is.na(result_i$Estimate), paste("(", round(exp(result_i$`CI_Down`), 5), " - ", round(exp(result_i$`CI_Up`), 5), ")", sep = ""), NA)
  # 排除包含NA值的行
  result_i <- result_i[complete.cases(result_i),]
  # 将结果合并到最终数据框中
  final_result_nosig <- rbind(final_result_nosig, result_i)
}
Pval<-c()
for (i in seq_along(Nosignificant_results)) {
  p_value_i <- Nosignificant_results[[i]]$result$`MR-PRESSO results`$`Global Test`$Pvalue
  Pval <- c(Pval,p_value_i)
}
final_result_nosig$Pval <- Pval
final_result_nosig$exposure<-extracted_data_nosig$exposure
final_result_nosig$outcome<-extracted_data_nosig$outcome
openxlsx::write.xlsx(final_result_nosig, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/deocde/PRO_nosig_result.xlsx",
                     sheetName = "Sheet1", Colnames = TRUE)

#step . save result 
save(exposure_dat,outcome_dat,dat,
     dat_sen, dat_sen1, dat_pro, #dat_pro1,
     file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/deocde/all_har.RData")
save(res_ivw,#mr_results_het,mr_results_ple,
     mr_result_pro,final_result_sig,final_result_nosig,
     filtered_data1,
     #dat_padj,final_result,
     file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/deocde/all_result.RData")


############ ukb-coloc #############
coloc_batch_ukbppp(
  gwas_path1 = "E:/coloc-ukbppp/protein/tlr1_q15399_oid30532_v1_inflammation_ii.tar",
  snp_annotation_dir="E:/coloc-ukbppp/snp/",
  pos_col = "POS38",
  gwas_annotation = gwas_annotation_ukbppp,
  gwas_path2="F:/post-GWAS/data_gwas_list/data_prepare/coloc/ER+_Breast_cancer_BCACIEU_hg38.txt",
  bfile = "E:/coloc-ukbppp/peitao/EUR",
  prepare_method = 1,
  SS2 = 175475,
  coloc_plot = T,
  coloc_plot_pph4 = 0.38,
  out_path = "E:/coloc-ukbppp/result/",
  NC2 = 69501)

# coloc_files <- dir("E:/coloc-ukbppp/result/PCD基因集/IEU结局/NSCLC/",pattern = "coloc_abf_summary.csv",full.names = T,recursive = T)
# coloc_summary <- data.frame()
# for (i in coloc_files) {
#   coloc_summary_temp <- read.csv(i)
#   coloc_summary <- rbind(coloc_summary,coloc_summary_temp)
# }
# write.csv(coloc_summary,file = "E:/coloc-ukbppp/result/PCD基因集/IEU结局/NSCLC/NSCLCcoloc_summary.csv",row.names = F)

############ decode-coloc #############
coloc_batch_decode(
  gwas_path1 = "E:/coloc-decode/protein/16324_38_TLR1_TLR1.txt.gz",
  gwas_path2="F:/post-GWAS/data_gwas_list/data_prepare/coloc/ER+_Breast_cancer_BCACIEU_hg38.txt",
  type1 = "quant",
  SS1 = 35559,
  bfile = "E:/药物靶点肺癌/配套/EUR",
  prepare_method = 1,
  type2 = "cc",
  SS2 = 175475,
  NC2 = 69501,
  coloc_plot = T,
  coloc_plot_pph4 = 0.38,
  out_path = "E:/coloc-decode/result/",
  )

# coloc_files <- dir("E:/coloc-ukbppp/result/PCD基因集/IEU结局/NSCLC/",pattern = "coloc_abf_summary.csv",full.names = T,recursive = T)
# coloc_summary <- data.frame()
# for (i in coloc_files) {
#   coloc_summary_temp <- read.csv(i)
#   coloc_summary <- rbind(coloc_summary,coloc_summary_temp)
# }
# write.csv(coloc_summary,file = "E:/coloc-ukbppp/result/PCD基因集/IEU结局/NSCLC/NSCLCcoloc_summary.csv",row.names = F)







#### 黑色素瘤（UKB、decode、FinnGen都没有对应蛋白的cis-pqtl） #### 
################# decode ##############
#step 1. read exposure data 
a <- read.csv(file = "F:/post-GWAS/data_mr/decode.csv",sep = ",",header = T)
colnames(a)
exposure_dat <- read_exposure_data("F:/post-GWAS/data_mr/decode.csv", 
                                   clump = F, sep = ",",
                                   phenotype_col = "gene.exposure",
                                   snp = "SNP", beta_col = "beta.exposure", 
                                   se_col = "se.exposure", 
                                   effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure", 
                                   pval_col = "pval.exposure", samplesize_col = "samplesize.exposure", 
                                   eaf_col = "eaf.exposure",
                                   chr_col = "chr.exposure",
                                   pos_col = "pos.exposure") 
exposure_dat <- exposure_dat %>% filter(mr_keep.exposure==TRUE)

#step 2. read outcome data 
a <- fread(file = "F:/post-GWAS/data_gwas_list/data_prepare/mr_with_phenotype/Breast_cancer_BCACIEU.txt")
head(a)
outcome_dat1<-read_outcome_data(filename="F:/post-GWAS/data_gwas_list/data_prepare/mr_with_phenotype/Breast_cancer_BCACIEU.txt", 
                                snps = exposure_dat$SNP,
                                sep = "\t",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",phenotype_col = "phenotype_col",
                                # ncase_col = "ncase_col",ncontrol_col = "ncontrol",
                                samplesize_col = "samplesize")
outcome_dat2<-read_outcome_data(filename="F:/post-GWAS/data_gwas_list/data_prepare/mr_with_phenotype/ER-_Breast_cancer_BCACIEU.txt", 
                                snps = exposure_dat$SNP,
                                sep = "\t",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",phenotype_col = "phenotype_col",
                                # ncase_col = "ncase_col",ncontrol_col = "ncontrol",
                                samplesize_col = "samplesize")
outcome_dat3<-read_outcome_data(filename="F:/post-GWAS/data_gwas_list/data_prepare/mr_with_phenotype/ER+_Breast_cancer_BCACIEU.txt", 
                                snps = exposure_dat$SNP,
                                sep = "\t",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "se",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",
                                pval_col = "pval",phenotype_col = "phenotype_col",
                                # ncase_col = "ncase_col",ncontrol_col = "ncontrol",
                                samplesize_col = "samplesize")


#step 4. harmonise 
outcome_dat<-rbind(outcome_dat1,outcome_dat2,outcome_dat3)
dat <- harmonise_data(exposure_dat, outcome_dat)
table(dat$exposure)

#step 5. caculate F-stat for each SNP 
dat$EAF2 <- (1 - dat$eaf.exposure) 
dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2) 
PVEfx <- function(BETA, MAF, SE, N){ 
  pve <- (2*(BETA^2)*MAF*(1 - MAF))/((2*(BETA^2)*MAF*(1 - MAF)) + ((SE^2)*2*N*MAF*(1 - MAF))) 
  return(pve) } 
dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N = dat$samplesize.exposure) 
dat$FSTAT <- ((dat$samplesize.exposure - 1 - 1)/1)*(dat$PVE/(1 - dat$PVE))
dat <- subset(dat,FSTAT>10)
dat <- subset(dat,MAF>0.01)

#step 6. MR analysis using Inverse variance weighted method 
res_ivw <- generate_odds_ratios(mr_res=mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_wald_ratio",
                                                              "mr_weighted_median","mr_simple_mode",
                                                              "mr_weighted_mode","mr_wald_ratio")
))
write.csv(res_ivw,file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/ivw_result.csv")
openxlsx::write.xlsx(res_ivw, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/ivw_result.xlsx",
                     sheetName = "Sheet1", colNames = TRUE)

#step 7. p-adj 
filtered_data <- res_ivw %>% filter(method==c("Inverse variance weighted","Wald ratio"))
name<-filtered_data$outcome[!duplicated(filtered_data$outcome)]
number_outcome <- filtered_data %>% group_by(outcome) %>% summarize(count = n()) 
number_exposure <- filtered_data %>% group_by(exposure) %>% summarize(count = n()) 
filtered_data1<-filtered_data
filtered_data1<-filtered_data%>%arrange(pval)
filtered_data1$adjusted_p <- NA
for (i in 1:length(name)) {
  b <- filtered_data1 %>% filter(outcome == name[i])
  a <- b[order(b$pval), ]
  plist <- a$pval
  p <- plist
  sorted_p_values <- sort(p)
  adjusted_p_values <- p.adjust(sorted_p_values, method = "fdr")
  filtered_data1$adjusted_p[filtered_data1$outcome == name[i]] <- adjusted_p_values
}
openxlsx::write.xlsx(filtered_data1, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/ivw_p_adj.xlsx",
                     sheetName = "Sheet1", colNames = TRUE) 


#step 8. Sensitive:pleiotropy analysis and heterogeneity test, heterogeneity (Inverse variance weighted), stiger-test, LOO 
# pval<0.05
dat_sen1<-res_ivw %>% filter(pval<0.05)
dat_sen11<-dat[which(dat$exposure%in%dat_sen1$exposure), ]
dat_sen<-dat_sen11
write.csv(dat_sen,file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/dat_sen.csv")

# heterogeneity test
mr_results_het <- mr_heterogeneity(dat_sen)
b1<-subset(res_ivw,pval<0.05)
b<-subset(b1,method==c("Inverse variance weighted","Wald ratio"))
mr_results_het_use<- merge(mr_results_het, b[, c("outcome", "exposure")], by = c("outcome", "exposure"))
openxlsx::write.xlsx(mr_results_het_use, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/het.xlsx",
                     sheetName = "Sheet1", Colnames = TRUE)

# pleiotropy analysis 
mr_results_ple <- mr_pleiotropy_test(dat_sen)
c1<-subset(res_ivw,pval<0.05)
c<-subset(c1,method==c("Inverse variance weighted","Wald ratio"))
mr_results_ple_use<- merge(mr_results_ple, c[, c("outcome", "exposure")], by = c("outcome", "exposure"))
openxlsx::write.xlsx(mr_results_ple_use, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/ple.xlsx",
                     sheetName = "Sheet1", Colnames = TRUE)

# steiger方向性检测 
steiger_whole_test <- directionality_test(dat =dat )
steiger_snp_test <- steiger_filtering(dat =dat )
write.csv(steiger_whole_test,file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/steiger.csv")

# MRPRESSO 
dat_sen1 <- res_ivw %>%
  filter(pval < 0.05, method == "Inverse variance weighted")
name <- dat_sen1 %>%
  count(exposure, outcome) %>%
  filter(n != 0) %>%
  select(exposure, outcome)
n <- nrow(name)

merged_data1 <- dat %>% filter(exposure==name$exposure[1]&outcome==name$outcome[1])
for(i in 2:n){new.data <- dat %>% filter(exposure==name$exposure[i]&outcome==name$outcome[i])
merged_data1 = rbind(merged_data1,new.data)}
dat_pro_tmp <- merged_data1 %>% group_by(exposure,outcome) %>% summarize(count = n()) %>% filter(count >= 4)
n <- nrow(dat_pro_tmp)

merged_data2 <- dat %>% filter(exposure==dat_pro_tmp$exposure[1]&outcome==dat_pro_tmp$outcome[1])
for(i in 2:n){new.data <- dat %>% filter(exposure==dat_pro_tmp$exposure[i]&outcome==dat_pro_tmp$outcome[i])
merged_data2 = rbind(merged_data2,new.data)}
dat_pro_tmp2 <- merged_data2 %>% group_by(exposure,outcome) %>% summarize(count = n())
dat_pro <- merged_data2

# mr_result_pro<-MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
#           BetaExposure = "beta.exposure",
#           SdOutcome = "se.outcome",
#           SdExposure = "se.exposure",
#           OUTLIERtest = TRUE,
#           DISTORTIONtest = TRUE,
#           data = dat_pro,
#           NbDistribution = 1000,
#           SignifThreshold = 0.05)
mr_result_pro<-run_mr_presso(dat_pro, NbDistribution = 1000, SignifThreshold = 0.05)
mr_results_pro <- mr_result_pro

significant_results <- list()
# 遍历每个结果块
for (i in seq_along(mr_results_pro)) {
  # 获取当前结果块的Global Test Pvalue
  p_value <- mr_results_pro[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  # 检查Pvalue是否小于0.05
  if (!is.na(p_value) && p_value < 0.05) {
    # 将符合条件的结果和位置添加到列表中
    significant_results[[length(significant_results) + 1]] <- list(result = mr_results_pro[[i]], position = i)
  }
}

Nosignificant_results <- list()
# 遍历每个结果块
for (i in seq_along(mr_results_pro)) {
  # 获取当前结果块的Global Test Pvalue
  p_value <- mr_results_pro[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  # 检查Pvalue是否小于0.05
  if (!is.na(p_value) && p_value > 0.05) {
    # 将符合条件的结果和位置添加到列表中
    Nosignificant_results[[length(Nosignificant_results) + 1]] <- list(result = mr_results_pro[[i]], position = i)
  }
}

extracted_data_sig <- data.frame(exposure = character(), outcome = character(), stringsAsFactors = FALSE)
for (i in seq_along(significant_results)){
  position <- significant_results[[i]]$position
  exposure <- dat_pro_tmp2$exposure[position]
  outcome <- dat_pro_tmp2$outcome[position]
  extracted_data_sig <- rbind(extracted_data_sig, data.frame(exposure = exposure, outcome = outcome, stringsAsFactors = FALSE))
}

extracted_data_nosig <- data.frame(exposure = character(), outcome = character(), stringsAsFactors = FALSE)
for (i in seq_along(Nosignificant_results)){
  position <- Nosignificant_results[[i]]$position
  exposure <- dat_pro_tmp2$exposure[position]
  outcome <- dat_pro_tmp2$outcome[position]
  extracted_data_nosig <- rbind(extracted_data_nosig, data.frame(exposure = exposure, outcome = outcome, stringsAsFactors = FALSE))
}

# 创建一个空的数据框来存储最终结果
final_result_sig <- data.frame(Estimate = numeric(0), Sd = numeric(0), CI_Up = numeric(0), CI_Down = numeric(0), 
                               OR_Estimate = character(0), OR_CI = character(0), stringsAsFactors = FALSE)
# 循环遍历significant_results中的每个元素
for (i in seq_along(significant_results)) {
  # 从每个子集中提取数据
  data_i <- significant_results[[i]]$result$`Main MR results`
  # 提取Causal Estimate和Sd
  estimate_i <- data_i$`Causal Estimate`[1]
  sd_i <- data_i$Sd[1]
  # 创建数据框并计算上下置信区间，排除NA值
  result_i <- data.frame(Estimate = estimate_i, Sd = sd_i)
  result_i$`CI_Up` <- ifelse(!is.na(result_i$Estimate), result_i$Estimate + 1.96 * result_i$Sd, NA)
  result_i$`CI_Down` <- ifelse(!is.na(result_i$Estimate), result_i$Estimate - 1.96 * result_i$Sd, NA)
  # 计算OR及其95%置信区间
  result_i$OR_Estimate <- ifelse(!is.na(result_i$Estimate), round(exp(result_i$Estimate), 5), NA)
  result_i$OR_CI <- ifelse(!is.na(result_i$Estimate), paste("(", round(exp(result_i$`CI_Down`), 5), " - ", round(exp(result_i$`CI_Up`), 5), ")", sep = ""), NA)
  # 排除包含NA值的行
  result_i <- result_i[complete.cases(result_i),]
  # 将结果合并到最终数据框中
  final_result_sig <- rbind(final_result_sig, result_i)
}
Pval<-c()
for (i in seq_along(significant_results)) {
  p_value_i <- significant_results[[i]]$result$`MR-PRESSO results`$`Global Test`$Pvalue
  Pval <- c(Pval,p_value_i)
}
final_result_sig$Pval <- Pval
final_result_sig$exposure<-extracted_data_sig$exposure
final_result_sig$outcome<-extracted_data_sig$outcome
openxlsx::write.xlsx(final_result_sig, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/PRO_sig_result.xlsx",
                     sheetName = "Sheet1", Colnames = TRUE)

# 创建一个空的数据框来存储最终结果
final_result_nosig <- data.frame(Estimate = numeric(0), Sd = numeric(0), CI_Up = numeric(0), CI_Down = numeric(0), 
                                 OR_Estimate = character(0), OR_CI = character(0), stringsAsFactors = FALSE)
# 循环遍历significant_results中的每个元素
for (i in seq_along(Nosignificant_results)) {
  # 从每个子集中提取数据
  data_i <- Nosignificant_results[[i]]$result$`Main MR results`
  # 提取Causal Estimate和Sd
  estimate_i <- data_i$`Causal Estimate`[1]
  sd_i <- data_i$Sd[1]
  # 创建数据框并计算上下置信区间，排除NA值
  result_i <- data.frame(Estimate = estimate_i, Sd = sd_i)
  result_i$`CI_Up` <- ifelse(!is.na(result_i$Estimate), result_i$Estimate + 1.96 * result_i$Sd, NA)
  result_i$`CI_Down` <- ifelse(!is.na(result_i$Estimate), result_i$Estimate - 1.96 * result_i$Sd, NA)
  # 计算OR及其95%置信区间
  result_i$OR_Estimate <- ifelse(!is.na(result_i$Estimate), round(exp(result_i$Estimate), 5), NA)
  result_i$OR_CI <- ifelse(!is.na(result_i$Estimate), paste("(", round(exp(result_i$`CI_Down`), 5), " - ", round(exp(result_i$`CI_Up`), 5), ")", sep = ""), NA)
  # 排除包含NA值的行
  result_i <- result_i[complete.cases(result_i),]
  # 将结果合并到最终数据框中
  final_result_nosig <- rbind(final_result_nosig, result_i)
}
Pval<-c()
for (i in seq_along(Nosignificant_results)) {
  p_value_i <- Nosignificant_results[[i]]$result$`MR-PRESSO results`$`Global Test`$Pvalue
  Pval <- c(Pval,p_value_i)
}
final_result_nosig$Pval <- Pval
final_result_nosig$exposure<-extracted_data_nosig$exposure
final_result_nosig$outcome<-extracted_data_nosig$outcome
openxlsx::write.xlsx(final_result_nosig, file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/PRO_nosig_result.xlsx",
                     sheetName = "Sheet1", Colnames = TRUE)

#step . save result 
save(exposure_dat,outcome_dat,dat,
     dat_sen, dat_sen1, dat_pro, #dat_pro1,
     file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/all_har.RData")
save(res_ivw,#mr_results_het,mr_results_ple,
     mr_result_pro,final_result_sig,final_result_nosig,
     filtered_data1,
     #dat_padj,final_result,
     file = "E:/post_gwas_cancer/results/乳腺癌药物靶点+共定位结果/ukbppp/all_result.RData")

