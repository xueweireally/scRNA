library(IOBR)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)
library(readxl)
library(writexl)
library(data.table)
library(ggpubr)

# 设置工作目录，沿用之前
out_dir<- '~/Desktop/CRC_Macrophage_Figure5_6/'
# 取一个IOBR分析项目名
proj.name <- 'CRC'

# 1.读入数据 -----
# 1.1 样本分组数据 ----
lasso_risk <- read.csv('~/Desktop/CRC_Macrophage_Figure5_6/Rds_02_Lasso_Cox_04_risk_group_surv.txt', sep = "\t", h = T, check.names = F)
pset <- lasso_risk %>%
  dplyr::select(sample_id, risk_group, risk_score) %>%
  dplyr::rename(ID = sample_id,
                group = risk_group) %>%
  dplyr::select(ID, group, risk_score)

# 1.2 表达数据 -----
#log2(TPM+1) gene行样本列
expr_dats <- read_rds("~/Desktop/CRC_Macrophage_Figure5_6/TCGA-COAD.tpm_clinical.rds")
eset <- expr_dats %>%
  as.matrix() %>%
  t()
mode(eset) <- "numeric"
eset <-eset[,pset$ID]

# pset和eset为IOBR免疫分析输入数据
# 2. 计算各种功能特征的富集分数 -----
## 2.1 查看有哪些功能特征集 ----
# 包含3类：TME-associated, tumor-metabolism, and tumor-intrinsic signatures
names(signature_tme)[1:20]
names(signature_metabolism)[1:20]
names(signature_tumor)
#所有特征集的参考文献查看
signature_collection_citation %>% View()
# 特征集的分组
sig_group

## 2.2 查看有哪些计算特征富集分数的方法 ----
signature_score_calculation_methods
# 一般使用ssgsea

## 2.3 计算免疫相关特征富集分数  ------
sig_res<-calculate_sig_score(pdata           = NULL,
                             eset            = eset,
                             # signature       = signature_tme,
                             signature       = signature_collection,
                             method          = "ssgsea",
                             mini_gene_count = 2)
write_xlsx(sig_res, file.path(out_dir, 'IOBR_signature_all.xlsx'))

# 3. 免疫特征与表型分组的相关性 -----
pset$group <- factor(pset$group, levels = c('low_risk', 'high_risk'))
res<-iobr_cor_plot(pdata_group           = pset,
                   id1                   = "ID",
                   feature_data          = sig_res,
                   id2                   = "ID",
                   target                = NULL,
                   group                 = "group",
                   is_target_continuous  = FALSE,
                   padj_cutoff           = 1,
                   index                 = 1,
                   category              = "signature",
                   # 关注哪些特征组
                   signature_group       = sig_group[c(1,2,3,5,6,7,8,28,29,30,31,40)],
                   ProjectID             = proj.name,
                   # 选择两组的颜色,多种期刊配色可选如"npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty"
                   # palette_box           = "paired1",
                   palette_box           = 'jama',
                   palette_corplot       = "pheatmap",
                   palette_heatmap       = 2,
                   feature_limit         = 26,
                   character_limit       = 30,
                   show_heatmap_col_name = FALSE,
                   show_col              = FALSE,
                   show_plot             = TRUE,
                   path = file.path('./', paste0("1-",proj.name, "-relevant-signatures")))

# 4. 最相关特征中是哪个基因与表型相关 -----
# 查看上面输出结果中，横坐标哪个特征在两组间差别很大，需要进一步关注
# 示例
pos_sig <- c( "Pan_F_TBRs", "EMT1","EMT2","EMT3","Fatty_Acid_Degradation","Fatty_Acid_Elongation","Fatty_Acid_Biosynthesis", "Hu_hypoxia_signature", "Nature_metabolism_Hypoxia")
# 可查看免疫检查点基因在两组表达的差异
res<-iobr_cor_plot(pdata_group           = pset,
                   id1                   = "ID",
                   feature_data          = eset,
                   id2                   = "ID",
                   target                = NULL,
                   group                 = "group",
                   is_target_continuous  = FALSE,
                   padj_cutoff           = 1,
                   index                 = 2,
                   category              = "gene",
                   signature_group       = signature_collection[pos_sig],
                   ProjectID             = proj.name,
                   palette_box           = "paired1",
                   palette_corplot       = "pheatmap",
                   palette_heatmap       = 4,
                   feature_limit         = 26,
                   character_limit       = 30,
                   show_heatmap_col_name = FALSE,
                   show_col              = FALSE,
                   show_plot             = TRUE,
                   path                  = file.path('./', paste0("2-",proj.name,"-relevant-genes")))

# *免疫检查点基因与risk score相关性 -----
ICs <- c('CD274', 'CTLA4', 'PDCD1', 'LAG3', 'TIGIT', 'HAVCR2')
for (IC in ICs) {
  cor_df <- data.frame(
    gene = eset[IC,] %>% as.numeric(),
    riskScore = pset %>%
      arrange(factor(ID, level = colnames(eset))) %>%
      pull(risk_score)
  )
  corplot1 <- ggscatter(cor_df, x = 'gene', y = 'riskScore',
                        add = 'reg.line',
                        conf.int = TRUE,
                        add.params = list(color = '#e41a1c', fill = '#fb8072'))+
    stat_cor(method = 'pearson')+
    theme_bw() +
    theme(axis.title = element_text(face = 'bold', size = 14),
          axis.text = element_text(face = 'bold', size = 12))
  if (!dir.exists(file.path('./', paste0("3-",proj.name,"-riskScore-checkpoint")))) {
    dir.create(file.path('./', paste0("3-",proj.name,"-riskScore-checkpoint")))
  }
  ggsave(file.path('./',paste0("3-",proj.name,"-riskScore-checkpoint"), paste0(IC, '_riskScore_cor.pdf')), corplot1, width = 4, height = 4)
}

# 5. 解卷积计算免疫细胞组成 -----
# TIMER计算需要提供癌种名称,根据自己研究而定！！
cancer.name <- 'COAD'
# 各工具推断免疫组分
cibersort<-deconvo_tme(eset = eset, method = "cibersort", arrays = FALSE, perm = 200 )
epic<-deconvo_tme(eset = eset, method = "epic", arrays = FALSE)
mcp<-deconvo_tme(eset = eset, method = "mcpcounter")
xcell<-deconvo_tme(eset = eset, method = "xcell",arrays = FALSE)
estimate<-deconvo_tme(eset = eset, method = "estimate")

# TIMER需要安装sva
# BiocManager::install("sva")
library(sva)
timer<-deconvo_tme(eset = eset, method = "timer", group_list = rep(cancer.name,dim(eset)[2]))

quantiseq<-deconvo_tme(eset = eset, tumor = TRUE, arrays = FALSE, scale_mrna = TRUE, method = "quantiseq")
ips<-deconvo_tme(eset = eset, method = "ips", plot= FALSE)
tme_combine<-cibersort %>%
  inner_join(.,mcp,by       = "ID") %>%
  inner_join(.,xcell,by     = "ID") %>%
  inner_join(.,epic,by      = "ID") %>%
  inner_join(.,estimate,by  = "ID") %>%
  inner_join(.,timer,by     = "ID") %>%
  inner_join(.,quantiseq,by = "ID") %>%
  inner_join(.,ips,by       = "ID")
write_xlsx(sig_res, file.path(out_dir, 'IOBR_deconvolution_tme.xlsx'))

# 每个免疫细胞组成推断工具的结果画图
plot_immune_cell <- function(tme_score, prefix, cell_col, outdir, w = 14, h = 8,  color.pal = c("#0072B5FF", "#BC3C29FF")) {
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  plot_df <- tme_score %>%
    dplyr::select(cell_col) %>%
    rename_with(~str_remove(.x, prefix), ends_with(prefix)) %>%
    gather(key = 'cell_type',
           value = 'value',
           -c('ID', 'group')
    )

  p <- ggboxplot(plot_df,
                 x = "cell_type", y = "value",
                 color = "group", palette = color.pal,
                 add = "jitter", add.params = list(size = 1), size = 0.7) +
    stat_compare_means(aes(group = group, label = ..p.signif..)) +
    xlab('')+
    ylab('Cell Fraction') +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = 'top',
          axis.text = element_text(size = 12, face = 'bold', colour = 'black'),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title = element_text(size = 15, face = 'bold'))
  ggsave(file.path(outdir, paste0('immune',prefix, '.pdf')), p, width = w, height = h)
  write_xlsx(
    list(data = plot_df),
    file.path(outdir, paste0('immune',prefix, '.xlsx')))
}

# 画图数据
tme_score <- pset %>%
  dplyr::select(-risk_score) %>%
  left_join(tme_combine)
tme_score$group <- factor(tme_score$group, levels = c('low_risk', 'high_risk'))

# 调用函数针对每个工具画图
# 主要调整输出图的宽度w和高度h
plot_immune_cell(tme_score = tme_score, prefix = paste0('_CIBERSORT'), c(1,2,3:24), outdir = file.path(out_dir, paste0('4-',proj.name,'-relevant-TME-cells')), w = 10, h = 6)
plot_immune_cell(tme_score = tme_score, prefix = '_MCPcounter', c(1,2,28:37), outdir = file.path(out_dir, paste0('4-',proj.name,'-relevant-TME-cells')), w = 9, h = 6) #注意都可以换颜色 , color.pal = c('#008080', '#FFC0CB')
plot_immune_cell(tme_score = tme_score, prefix = '_xCell', c(1,2,38:101), outdir = file.path(out_dir, paste0('4-',proj.name,'-relevant-TME-cells')), w = 20, h = 6, color.pal = c('#2A9D8F', '#E76F51'))
plot_immune_cell(tme_score = tme_score, prefix = '_EPIC', c(1,2,105:111), outdir = file.path(out_dir, paste0('4-',proj.name,'-relevant-TME-cells')), w = 6, h = 6)
plot_immune_cell(tme_score = tme_score, prefix = '_estimate', c(1,2,113:115), outdir = file.path(out_dir, paste0('4-',proj.name,'-relevant-TME-cells')), w = 6, h = 6, color.pal = c('#2A9D8F', '#E76F51'))
 plot_immune_cell(tme_score = tme_score, prefix = '_TIMER', c(1,2,117:122), outdir = file.path(out_dir, paste0('4-',proj.name,'-relevant-TME-cells')), w = 8, h = 6)
plot_immune_cell(tme_score = tme_score, prefix = '_quantiseq', c(1,2,123:132), outdir = file.path(out_dir, paste0('4-',proj.name,'-relevant-TME-cells')), w = 8, h = 6)
plot_immune_cell(tme_score = tme_score, prefix = '_IPS', c(1,2,134:139), outdir = file.path(out_dir, paste0('4-',proj.name,'-relevant-TME-cells')), w = 8, h = 6)





