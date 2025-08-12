library(DrugTargetMR)

### BACA 功效计算 ###

# 定义TWAS功效计算函数
calculate_twas_power <- function(h2, h2_ge, N, alpha=5e-08) {
  # 参数说明:
  # h2: 表型遗传力 (来自BCAC GWAS)
  # h2_ge: 基因表达遗传力 (乳腺组织推荐值)
  # N: 样本量
  # alpha: 显著性阈值 (默认0.05)
  
  # 计算非中心参数 (NCP)
  lambda <- N * h2 * h2_ge
  
  # 计算卡方检验临界值 (df=1)
  crit_val <- qchisq(1 - alpha, df=1)
  
  # 计算统计功效
  power <- 1 - pchisq(crit_val, df=1, ncp=lambda)
  
  return(power)
}

# BCAC实际参数应用
bcac_power <- calculate_twas_power(
  h2 = 0.1344,       # BCAC GWAS遗传力估计
  h2_ge = 0.2,       # 乳腺组织基因表达遗传力(调整后)
  N = 228951         # BCAC样本量
)

cat(sprintf("BCAC TWAS预计功效: %.2f%%\n", bcac_power*100))

# 敏感性分析：样本量对功效的影响
sample_sizes <- seq(10000, 228951, by=1000)  # 更合理的样本量范围

# 计算不同遗传力场景下的功效曲线
power_curves <- list(
  "h2=0.1344 (BCAC)" = sapply(sample_sizes, function(n) 
    calculate_twas_power(h2=0.1344, h2_ge=0.2, N=n)),
  # 
  # "h2=0.2 (保守估计)" = sapply(sample_sizes, function(n) 
  #   calculate_twas_power(h2=0.2, h2_ge=0.2, N=n)),
  # 
  # "h2=0.3 (高遗传力)" = sapply(sample_sizes, function(n) 
  #   calculate_twas_power(h2=0.3, h2_ge=0.2, N=n))
)

# 高级可视化
library(ggplot2)
library(reshape2)

# 准备绘图数据
plot_data <- melt(data.frame(
  N = sample_sizes,
  power_curves
), id.vars = "N")

# 绘制功效曲线
pdf(file = "E:/post_gwas_cancer/tmp2/BACA_BC_TWAS功效图.pdf", 
    width = 15, height = 8,  # 宽幅设置
    pointsize = 18)  # 全局字体基准
p<-ggplot(plot_data, aes(x = N/1000, y = value, color = variable)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  annotate("text", x = max(sample_sizes)/1000, y = 0.82, 
           label = "80% Power Threshold", color = "red") +
  scale_x_continuous("Sample Size (thousands)", 
                     breaks = seq(0, 250, by = 50)) +
  scale_y_continuous("Power", labels = scales::percent) +
  ggtitle("Breast Cancer TWAS Power Analysis",
          subtitle = "Varying Heritability Scenarios (h2_ge=0.2)") +
  labs(color = "Heritability Scenario") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.8, 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
print(p)
dev.off()

# 输出当前参数下的功效表
power_table <- data.frame(
  N = c(50000, 100000, 150000, 228951),
  Power = sapply(c(50000, 100000, 150000, 228951), function(n) 
    calculate_twas_power(h2=0.1344, h2_ge=0.2, N=n))
)

knitr::kable(power_table, 
             col.names = c("Sample Size", "Power"),
             caption = "BCAC TWAS Power at Different Sample Sizes")
