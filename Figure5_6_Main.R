# library ----
require(dplyr)
require(tibble)
require(tidyr)
require(dplyr)
require(stringr)
require(readr)
require(readxl)
require(writexl)
require(ggplot2)
require(ggpubr)

# work dir -----
dir1 <- '~/Desktop/CRC_Macrophage_Figure5_6/'
dir.create(dir1)

count_df <- read_rds("~/Desktop/CRC_Macrophage_Figure5_6/TCGA-COAD.count.rds")
tpm <- read_rds("~/Desktop/CRC_Macrophage_Figure5_6/TCGA-COAD.tpm_clinical.rds")
clin_df <- read_rds("~/Desktop/CRC_Macrophage_Figure5_6/TCGA-COAD.group_clinical.rds")
risk_df <- read.csv('~/Desktop/CRC_Macrophage_Figure5_6/Rds_02_Lasso_Cox_04_risk_group_surv.txt', sep = "\t", h = T, check.names = F)

source('~/Desktop/CRC_Macrophage_Figure5_6/My_Funcs.R')

sample.group <- risk_df %>%
  dplyr::rename(sample = sample_id) %>%
  dplyr::rename(group = risk_group)
risk_deseq2 <-
  run_deseq2(count_df, sample.group, control.group = 'low_risk')
write_rds(risk_deseq2, file.path(dir1, 'risk_deseq2.rds'))
# dysregulated genes
risk_sig_deseq2 <-
  plotutil.expr.volcano(
    deseq2.res = risk_deseq2,
    log2FC.cutoff = c(-1.58, 0.58), #1.5倍表达
    # log2FC.cutoff = c(-1, 1), #2倍差异表达
    padj.cutoff = 0.05,
    show.num = 5,
    color.pal = c("#B15B00", "#008C81")
  )
ggsave(file.path(dir1, 'risk_group_deseq2_volcano.pdf'), risk_sig_deseq2$volcano.plot, height = 4, width = 5)
risk_sig_deseq2$deseq2.sig %>% View()
risk_sig_deseq2$deseq2.sig %>% dim() #56014     8
write_xlsx(risk_sig_deseq2$deseq2.sig, file.path(dir1, 'risk_group_deseq2_FC1.5.xlsx'))

# 高低风险两组通路富集分析 ----
risk_deseq2 <- read_xlsx(file.path(dir1, 'risk_group_deseq2_FC1.5.xlsx'))
## GO和KEGG富集（超几何检验）-----
#高风险 vs 低风险差异表达基因列表,用于GO和KEGG富集分析输入
sig_genes <- risk_deseq2 %>%
  filter(group %in% c('up', 'down'))
sig_genes$group %>% table()
# down   up
# 16 2779
sig_genes <- sig_genes %>%
  pull(gene)
GO_KEGG_res <- enrich.GO.KEGG(sig.genes = sig_genes, GO.db = 'org.Hs.eg.db', KEGG.db = 'hsa',p.cutoff = 0.05, q.cutoff = 0.05, dot.color = c("#FF5733","#127475"))
# 保存GO和KEGG富集结果
write_xlsx(GO_KEGG_res[1:2], file.path(dir1, 'GO_KEGG_enrich_result_FC1.5.xlsx'))
# 保存GO和KEGG富集结果图
ggsave(filename = file.path(dir1, 'GO_barplot_FC1.5.pdf'), GO_KEGG_res$GO.barplot,
       height = 7.5, width = 7)

## GSEA富集（KS检验），计算通路富集分数ES-----
# GSEA分析函数输入为deseq2差异表
# 输出:1.GSEA对象，可用于后续挑选通路画图；2. GSEA富集结果表；3. 山脊图
# 图中横轴为富集到通路中的基因的log2FC分布，因此山脊大于0为通路在高风险中显著上调，小于0为通路在高风险中显著下调
h.gmt <- '~/Desktop/CRC_Macrophage_Figure5_6/h.all.v2023.2.Hs.symbols.gmt'
GSEA_res <- enrich.GSEA(deseqs.res = risk_deseq2, OrgDb = "org.Hs.eg.db", gmt.file = h.gmt, kegg.db = 'hsa',fill.color = c("#FF5733","#127475"),p.cutoff = 0.05)
write_rds(GSEA_res$GSEA.obj, file.path(dir1, 'GSEA_enrich_result_FC1.5.rds')) #后续可复用的GSEA结果对象，只能保存为rds形式
write_xlsx(GSEA_res$GSEA.df, file.path(dir1, 'GSEA_enrich_result_FC1.5.xlsx'))
ggsave(filename = file.path(dir1, 'GSEA_ridgeplot_FC2.pdf'), GSEA_res$GSEA.ridgeplot, height = 10, width = 8)

# 读入之前的GSEA富集结果，选择想展示的通路
GSEA_df <- read_xlsx(file.path(dir1, 'GSEA_enrich_result_FC1.5.xlsx'))
GSEA_res <- read_rds(file.path(dir1, 'GSEA_enrich_result_FC1.5.rds'))
# 调用函数画gseaplot
gseap <- plotutil.gseaplot(gsea.res = GSEA_res, plot.term = c(1,2,3,4), title.name = 'Hallmark')
ggsave(file.path(dir1, 'gseaplot.pdf'), gseap, height = 6, width = 7)

# after METAFlux_CRC.R ----
met_ptws <- read_rds('~/Desktop/CRC_Macrophage_Figure5_6/METAFlux_pathway_142.rds')
met_flux <- read_rds('~/Desktop/CRC_Macrophage_Figure5_6/METAFlux_flux_469.rds')

met_flux <- met_flux %>%
  as.data.frame() %>%
  rownames_to_column('flux')
met_flux <- data.table::transpose(met_flux, keep.names = 'sample_id', make.names = 1) %>%
  left_join(dplyr::select(risk_df,sample_id, risk_group))
flux_wilcox <- batch.wilcoxon.test(test.columns = colnames(met_flux)[2:(ncol(met_flux)-1)], group_col = 'risk_group', data = met_flux)
flux_wilcox <- flux_wilcox %>%
  rename(flux = drug)
write_xlsx(flux_wilcox, file.path(dir1, 'metabolic_flux_wilcox.xlsx')) #37个小于0.01, 166小于0.05
flux_wilcox <- read_xlsx(file.path(dir1, 'metabolic_flux_wilcox.xlsx'))
# flux对应反应注释文件
flux_anno <- read_rds('~/Desktop/CRC_Macrophage_Figure5_6/METAFlux_human_gem.rds')
sig_flux <- flux_wilcox %>%
  filter(p.adjust<0.05) %>%
  left_join(dplyr::select(flux_anno, ID, EQUATION, SUBSYSTEM), by = c('flux' = 'ID'))
write_xlsx(sig_flux, '~/Desktop/CRC_Macrophage_Figure5_6/metabolic_sig_flux.xlsx')

# 展示特定flux之间差异

## 代谢通路差异 ----
met_ptws <- met_ptws %>%
  rownames_to_column('pathway')
met_ptws <- data.table::transpose(met_ptws, keep.names = 'sample_id', make.names = 1) %>%
  left_join(dplyr::select(risk_df,sample_id, risk_group))
ptws_wilcox <- batch.wilcoxon.test(test.columns = colnames(met_ptws)[2:(ncol(met_ptws)-1)], group_col = 'risk_group', data = met_ptws)
ptws_wilcox <- ptws_wilcox %>%
  rename(pathway = drug)
write_xlsx(ptws_wilcox, file.path(dir1, 'metabolic_ptws_wilcox.xlsx'))
#富集叶酸代谢和脂质代谢
ptws_wilcox <- read_xlsx(file.path(dir1, 'metabolic_ptws_wilcox.xlsx'))

plot_df <- met_ptws
plot_df$`Beta oxidation of fatty acids` <- plot_df$`Beta oxidation of di-unsaturated fatty acids (n-6) (peroxisomal)`
# plot_df$`Beta oxidation of fatty acids` <- plot_df$`Beta oxidation of even-chain fatty acids (peroxisomal)`
# plot_df$`Beta oxidation of fatty acids` <- plot_df$`Beta oxidation of unsaturated fatty acids (n-9) (peroxisomal)`

plot_df <- plot_df %>%
  filter(!is.na(risk_group)) %>%
  mutate(risk_group = ifelse(risk_group == 'high_risk', 'high', 'low'))
plot_df$risk_group <- factor(plot_df$risk_group, levels = c('low', 'high'))
plotviolin <- function(plot_df, y){
  p <- ggviolin(plot_df, x = "risk_group", y = y,color = "risk_group", palette = get_palette('lancet',2),
                add = "boxplot",
                # facet.by = "type",
                short.panel.labs = TRUE)
  # Use only p.format as label. Remove method name.
  p1 <- p + stat_compare_means(label = "p.format", label.x = 1.3) +
    xlab('')+
    ylab('')+
    theme(axis.text = element_text(face = 'bold'),
          # axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(face = 'bold'),
          legend.position = 'none'
    )
  return(p1)
}

p1 <- plotviolin(plot_df, y = "Beta oxidation of fatty acids")
p2 <- plotviolin(plot_df, y = "Folate metabolism")
p3 <- plotviolin(plot_df, y = "Sulfur metabolism")
p4 <- plotviolin(plot_df, y = "Inositol phosphate metabolism")
p5 <- plotviolin(plot_df, y = "Ether lipid metabolism")

ggsave(
  file.path(dir1, 'Meta_ptws_violins.pdf'),
  ggarrange(
    plotlist = list(p1,p2,p5,p3,p4),
    nrow = 1,
    labels = c("Beta oxidation of fatty acids", "Folate metabolism", "Sulfur metabolism", "Inositol phosphate metabolism", "Ether lipid metabolism"),
    font.label = list(size = 12, color = "black", face = "bold", family = NULL)
  ),
  height = 4,
  width = 15
)

