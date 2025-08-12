#### DGIdb_PharmGKB_DrugCentral药物筛选####
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
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(reshape2)

# 通过三大数据库，获得基因-药物配对关系 ------------------------------------------
######### DGIdb
dgidb <- fread("F:/post-GWAS/data_DGIdb_PharmGKB_DrugCentral/DGIdb/latest-20250103获取/interactions.tsv",data.table = F)
dgidb$database  <- "DGIdb"
dgidb$drug_name
dgidb$gene_name

######### DrugCentral 
drugcenter <- fread("F:/post-GWAS/data_DGIdb_PharmGKB_DrugCentral/DrugCentral/drug.target.interaction.tsv.gz",data.table = F)
drugcenter <- drugcenter[!drugcenter$DRUG_NAME == "",]
drugcenter$DRUG_NAME <- toupper(drugcenter$DRUG_NAME)
drugcenter <- separate_rows(drugcenter, GENE, sep = "\\|")
colnames(drugcenter)[colnames(drugcenter) =="DRUG_NAME" ] <- "drug_name"
colnames(drugcenter)[colnames(drugcenter) =="GENE" ] <- "gene_name"
drugcenter$database <- "drugcenter"

######### pharmgkb
pharmgkb <- fread("F:/post-GWAS/data_DGIdb_PharmGKB_DrugCentral/PharmGKB/relationships.tsv",data.table = F)
pharmgkb <- pharmgkb[pharmgkb$Entity1_type == "Chemical" & pharmgkb$Entity2_type == "Gene",]
pharmgkb$drug_name <- toupper(pharmgkb$Entity1_name)
pharmgkb$gene_name <- toupper(pharmgkb$Entity2_name)
pharmgkb$database <- "pharmgkb"


# 三个数据库合并，得药物名-基因名-数据库得汇总列表。
common_col <- c("drug_name","gene_name","database")
drug_to_gene <- rbind(pharmgkb[,common_col],drugcenter[,common_col],dgidb[,common_col])
drug_to_gene <- drug_to_gene[!drug_to_gene$gene_name == "",]
drug_to_gene <- as.data.frame(apply(drug_to_gene, 2, toupper))
drug_to_gene <- drug_to_gene[!duplicated(drug_to_gene),]
write.csv(drug_to_gene,file = "E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/药靶匹配/药理基因drug_to_gene.csv",row.names = F)






# 选择满足4种方法均阳性的基因，用于匹配药物 ------------------------------------
######## 1. 药物→基因
drug_to_gene <- read.csv("E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/药靶匹配/药理基因drug_to_gene.csv")
drug_to_gene <- drug_to_gene[,c("drug_name","gene_name" )]
drug_to_gene <- drug_to_gene[!duplicated(drug_to_gene),]

######## 2. 基于基因进行匹配，研究性状→基因→药物
BC_gene_drug <- drug_to_gene %>% filter(gene_name %in% c("ADCY3","CASP8","CENPO","GRHL1","HELQ","TLR1"))


# 将匹配到的药物进行写出
write.csv(BC_gene_drug,file = "E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/乳腺癌-基因-药物.csv",row.names = F)




# 选取药物 ---------------------------------------------------------------------
Antihypertensive_drugs<- c("Amlodipine",
                           "Nifedipine",
                           "Atenolol",
                           "Metoprolol",
                           "Propranolo",
                           "Candesartan",
                           "Olmesartan",
                           "Telmisartan",
                           "Valsartan",
                           "Perindopri",
                           "Hydrochlorothiazide")
Lipid_lowering_drugs<- c("Atorvastatin",
                         "Fluvastatin",
                         "Lovastatin",
                         "Pitavastatin",
                         "Simvastatin",
                         "Bezafibrate",
                         "Fenofibrate",
                         "Gemfibrozi")
Glucose_lowering_drugs <- c("Metformin",
                            "Alogliptin",
                            "Sitagliptin")
Antiarrhythmics<- "Digoxin"
Antithrombotics <- c("Aspirin",
                     "Clopidogrel",
                     "Ticagrelor",
                     "Warfarin")
Antioxidants <- c("Vitamin E",
                  "Tocotrienol",
                  "Tanshinone",
                  "Resveratrol",
                  "Quercetin",
                  "Curcumin")
drug_choose <- data.frame(drug_group = c(rep("Antihypertensive_drugs",length(Antihypertensive_drugs)),
                                         rep("Lipid_lowering_drugs",length(Lipid_lowering_drugs)),
                                         rep("Glucose_lowering_drugs",length(Glucose_lowering_drugs)),
                                         rep("Antiarrhythmics",length(Antiarrhythmics)),
                                         rep("Antithrombotics",length(Antithrombotics)),
                                         rep("Antioxidants",length(Antioxidants))),
                          drug_name = c(Antihypertensive_drugs,
                                        Lipid_lowering_drugs,
                                        Glucose_lowering_drugs,
                                        Antiarrhythmics,
                                        Antithrombotics,
                                        Antioxidants))

# 将数据写出到csv文件
write.csv(drug_choose,file = "E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/选取33种药物.csv",row.names = F)


# 打分 -------------------------------------------------------------------------
#### 病理基因与病理term ####
### 病理基因
ts4 <- openxlsx::read.xlsx("E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/TWAS腹主动脉瘤和代谢特征 附表4.xlsx",
                           startRow = 2)
ts4 <- ts4[grepl("AAA",ts4$`Cross-traits`),]


#### 病理term
DEG_raw <- openxlsx::read.xlsx(xlsxFile = "E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌GO_KEGG富集分析结果/Genelist_enrich_analysis.xlsx",
                               sheet = "GeneMANIA")
Genes_BC1 <- DEG_raw$ADCY3
Genes_BC2 <- DEG_raw$CASP8
Genes_BC3 <- DEG_raw$CENPO
Genes_BC4 <- DEG_raw$GRHL1
Genes_BC5 <- DEG_raw$HELQ
Genes_BC6 <- DEG_raw$TLR1
Genes_BC <- rbind(Genes_BC1,Genes_BC2,Genes_BC3,Genes_BC4,Genes_BC5,Genes_BC6)
Genes <- bitr(Genes_BC,  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

KEGG <- enrichKEGG(gene = Genes$ENTREZID,
                   organism = "hsa", 
                   keyType = "kegg", #KEGG数据库
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.15,
                   qvalueCutoff = 0.2)

## 将对应的结果输出到文件
KEGG_raw_top <- KEGG@result[1:50,]
KEGG_raw_top1 <- KEGG@result
# 手动添加trait
write.csv(KEGG_raw_top,file = "E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/乳腺癌_富集分析_病理通路top50.csv",row.names = F)
write.csv(KEGG_raw_top1,file = "E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/乳腺癌_富集分析_病理通路top501.csv",row.names = F)


#### 药理基因与药理term
### 纳入的药物, 选取的drug列表
drugs1 <- read.csv("E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/可以用 的结果/乳腺癌选取21种药物.csv")
drugs1 <- drugs1$drug_name


### 药理基因
drug_to_gene <- read.csv("E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/药靶匹配/药理基因drug_to_gene.csv")


### 药理term
enrich_phar01 <- data.frame()
for (temp_drug in drugs1 ) {
  
  # 按照药物进行寻找，计算药物影响的term
  print(temp_drug)
  
  # 药物作用的基因
  temp_drug2gene <- drug_to_gene[grepl(temp_drug,drug_to_gene$drug_name,ignore.case = T),]
  
  if (nrow(temp_drug2gene) ==0) {
    next
  }
  
  # 药理term
  gene_phar <-   bitr(geneID = unique(temp_drug2gene$gene_name),fromType ="SYMBOL" ,toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  enrich_phar <- enrichKEGG(gene = unique(gene_phar$ENTREZID),
                            organism = "hsa",  #人
                            pAdjustMethod ="BH",  
                            keyType = "kegg",
                            pvalueCutoff =0.05,#设置pvalue界值
                            qvalueCutoff = 0.2 #设置qvalue界值(FDR校正后的p值）
  ) 
  
  # enrichKEGG()返回的富集通路就是按照p值升序排序的，选取top50 term
  enrich_phar <- enrich_phar@result[1:50,]
  enrich_phar$drug <- temp_drug
  
  enrich_phar01 <- rbind(enrich_phar01,enrich_phar)
}




## 将对应的结果输出到文件
write.csv(enrich_phar01,file = "E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/乳腺癌_富集分析_药理通路top50.csv",row.names = F)



#### 配对打分
### 按照性状对进行药物打分
enrich_phar01 <- read.csv("E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/乳腺癌_富集分析_药理通路top50.csv")
enrich_patho01 <- read.csv("E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/乳腺癌_富集分析_病理通路top50.csv")

enrichment_df1 <- data.frame()
for (temp_trait in unique(enrich_patho01$trait))  {
  print(temp_trait)
  
  for (temp_drug in unique(enrich_phar01$drug) ) {
    print(temp_drug)
    
    temp_gene_patho <- enrich_patho0[enrich_patho01$trait == temp_trait,]
    temp_drug_phar <- enrich_phar0[enrich_phar01$drug == temp_drug,]
    
    # 打分  
    df <- data.frame(Trait_Pair =temp_trait,
                     Drug = temp_drug,
                     pairing_score = length(intersect(temp_gene_patho$ID[1:20],temp_drug_phar$ID[1:20]))/20 + 
                       (sum( temp_gene_patho$ID[1:20] %in% temp_drug_phar$ID[21:50]))/20 +
                       sum((na.omit(temp_drug_phar$ID[1:20]) %in% na.omit(temp_gene_patho$ID[21:50])))/20
    )
    
    enrichment_df1 <- rbind(df,enrichment_df1)
  }
}

## 将对应的结果输出到文件
write.csv(enrichment_df1,file = "E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/乳腺癌_富集分析_pairing_score.csv",row.names = F)


#### 绘制热图-打分可视化
#############  加载R包
library(pheatmap)
library(reshape2)

#############  1. 长数据转宽数据
enrichment_df <- read.csv("E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/可以用 的结果/乳腺癌_富集分析_pairing_score.csv")
enrichment_df$pairing_score
data_matrix <- reshape2::dcast(enrichment_df, Drug ~ Trait_Pair, value.var = "pairing_score")
head(data_matrix)

# 设置行名为Gene.symbol，删除Gene.symbol列
row.names(data_matrix) <- data_matrix$Drug
data_matrix$Drug <- NULL

############# 2. pheatmap作图
#### 设置图例
# 定义颜色渐变，使用RColorBrewer提供的颜色。结合数据，确定图例的start、end和图例标记
# start为数据中值的最小值，为了展示好看，也可以设置为0；
# end为数据中值的最大值，为了展示好看，可以设置为整数。
start <- 0
end <- max(enrichment_df$pairing_score)

# 设置图块的颜色。colorRampPalette的作用是选择任意个颜色，比如这里选择20个颜色
breaksList <- seq(start, end, length.out = 20)
color <- colorRampPalette(c("white","royalblue"))(20)


# 设置图例标记，by是步长，可根据情况调整
legend_breaks <-  seq(start, end, by = 0.2)
legend_labels <-  seq(start, end, by = 0.2)


# 行注释
annotation_row <- read.csv("E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/可以用 的结果/乳腺癌选取21种药物.csv")
colnames(annotation_row) <- c("Group","Name")
row.names(annotation_row) <- annotation_row$Name
annotation_row$Name <- NULL

# 行注释对应的颜色设置
drug_group_color <- ggsci::pal_npg()(7)
names(drug_group_color) <- unique(annotation_row$Group)
annotation_colors <- list(Group=drug_group_color) #颜色设置


# data_matrix与annotation_row 行名保持一致
data_matrix <- data_matrix[row.names(annotation_row),]
identical(rownames(data_matrix),rownames(annotation_row))


# 图中数字展示
data_number <- data_matrix
data_number[data_number < 0.4] <- ""

# pheatmap绘制热图(横向)
pheatmap(data_matrix, 
         scale = "none",
         
         # # 输出pdf文件相关的参数。如果不输出文件，则删除或者注释filename参数
         filename = "E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/可以用 的结果/乳腺癌heatmap_2.pdf", # filename，定义输出文件的文件名
         height = 8, # height，定义输出文件的高度
         width = 20,   # width，定义输出文件的宽度
         
         color = color, 
         na_col = "grey",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = "grey60", # 单元格边框
         
         # 方块中数字显示
         display_numbers = data_number,
         fontsize_number = 6,
         
         # 标题
         main = "Heatmap of pairing score",
         
         # 行注释
         annotation_row  = annotation_row,
         annotation_colors = annotation_colors,
         
         # 图例
         legend = T,
         breaks = breaksList,
         legend_breaks = legend_breaks,
         legend_labels   = legend_labels,
         
         angle_col = 45,  # 设置列标签的旋转角度为45度
         fontsize_col = 8,  # 调整列标签字体大小
         cellwidth = 15,  # 增加单元格宽度
         cellheight = 15  # 增加单元格高度
) 


enrichment_df <- read.csv("E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/可以用 的结果/乳腺癌_富集分析_pairing_score.csv")
enrichment_df$pairing_score
data_matrix <- reshape2::dcast(enrichment_df, Drug ~ Trait_Pair, value.var = "pairing_score")
head(data_matrix)

# 设置行名为Gene.symbol，删除Gene.symbol列
row.names(data_matrix) <- data_matrix$Drug
data_matrix$Drug <- NULL

############# 2. pheatmap作图
#### 设置图例
# 定义颜色渐变，使用RColorBrewer提供的颜色。结合数据，确定图例的start、end和图例标记
# start为数据中值的最小值，为了展示好看，也可以设置为0；
# end为数据中值的最大值，为了展示好看，可以设置为整数。
start <- 0
end <- max(enrichment_df$pairing_score)

# 设置图块的颜色。colorRampPalette的作用是选择任意个颜色，比如这里选择20个颜色
breaksList <- seq(start, end, length.out = 20)
color <- colorRampPalette(c("white","royalblue"))(20)


# 设置图例标记，by是步长，可根据情况调整
legend_breaks <-  seq(start, end, by = 0.2)
legend_labels <-  seq(start, end, by = 0.2)


# 行注释
annotation_row <- read.csv("E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/可以用 的结果/乳腺癌选取21种药物.csv")
colnames(annotation_row) <- c("Group","Name")
row.names(annotation_row) <- annotation_row$Name
annotation_row$Name <- NULL

# 行注释对应的颜色设置
drug_group_color <- ggsci::pal_npg()(7)
names(drug_group_color) <- unique(annotation_row$Group)
annotation_colors <- list(Group=drug_group_color) #颜色设置


# data_matrix与annotation_row 行名保持一致
data_matrix <- data_matrix[row.names(annotation_row),]
identical(rownames(data_matrix),rownames(annotation_row))


# 图中数字展示
data_number <- data_matrix
data_number[data_number < 0.4] <- ""

# pheatmap绘制热图(横向)
pheatmap(data_matrix, 
         scale = "none",
         
         # # 输出pdf文件相关的参数。如果不输出文件，则删除或者注释filename参数
         filename = "E:/post_gwas_cancer/results/思路：跨组织分析/乳腺癌/乳腺癌结果/乳腺癌DGIdb_PharmGKB_DrugCentral药物筛选/可以用 的结果/乳腺癌heatmap_2.pdf", # filename，定义输出文件的文件名
         height = 8, # height，定义输出文件的高度
         width = 20,   # width，定义输出文件的宽度
         
         color = color, 
         na_col = "grey",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = "grey60", # 单元格边框
         
         # 方块中数字显示
         display_numbers = data_number,
         fontsize_number = 6,
         
         # 标题
         main = "Heatmap of pairing score",
         
         # 行注释
         annotation_row  = annotation_row,
         annotation_colors = annotation_colors,
         
         # 图例
         legend = T,
         breaks = breaksList,
         legend_breaks = legend_breaks,
         legend_labels   = legend_labels,
         
         angle_col = 45,  # 设置列标签的旋转角度为45度
         fontsize_col = 8,  # 调整列标签字体大小
         cellwidth = 15,  # 增加单元格宽度
         cellheight = 15  # 增加单元格高度
) 


