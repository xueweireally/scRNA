

#' DESeq2 Differential Expression Analysis
#'
#' @param expr.count A count expression matrix with the first column as "gene" and subsequent columns as samples
#' @param sample.group A data frame containing a "sample" column and a "group" column
#' @param control.group Reference group name
#'
#' @return A data frame of DESeq2 results
#' @export
#' @import DESeq2
#'
#' @examples
run_deseq2 <- function(expr.count, sample.group, control.group) {
  require(DESeq2)
  countMatrix <- expr.count %>%
    column_to_rownames("gene") %>%
    as.matrix()
  countMatrix <- countMatrix[, sample.group[["sample"]]]
  # 重要：调整列顺序与sample表的行样本名一致
  sample.group <- sample.group %>%
    dplyr::rename(condition = group)
  dds <-
    DESeqDataSetFromMatrix(round(countMatrix),
                           colData = sample.group,
                           design = ~ condition)
  dds$condition <-
    relevel(dds$condition, ref = control.group) # 指定哪一组作为control
  dds <- dds[rowSums(counts(dds)) > 1, ] # 过滤掉那些 count 结果都为 0 的数据
  dds <- DESeq(dds)
  deseq2_res <- results(dds)
  deseq2_res <- as.data.frame(deseq2_res) %>%
    rownames_to_column("gene")
  return(deseq2_res)
}

#' Visualization of DESeq2 Differential Expression Analysis Results
#'
#' @param deseq2.res A data frame of DESeq2 differential expression result
#' @param log2FC.cutoff The log2FoldChange threshold of significant differentially expressed genes. Defult c(-1, 1) means values outside this range are considered significantly differentially expressed
#' @param padj.cutoff The 'padj' threshold of significant differentially expressed genes. Default 0.05
#' @param show.num Number of upregulated and downregulated genes labeled in the plot. Default n=3
#' @param color.pal colors of upregulated and downregulated genes. Default c("#fa6e57", "#4695d6")
#'
#' @return A list containing DESeq2 result including column of upregulated and downregulated groups, and a volcano plot
#' @export
#' @import ggrepel
#'
#' @examples
plotutil.expr.volcano <-
  function(deseq2.res,
           log2FC.cutoff = c(-1, 1),
           padj.cutoff = 0.05,
           show.num = 3,
           color.pal = c("#fa6e57", "#4695d6")) {
    require(ggrepel)
    plot_data <- deseq2.res %>%
      mutate(group = case_when(
        (padj < padj.cutoff) & (log2FoldChange > log2FC.cutoff[2]) ~ 'up',
        (padj < padj.cutoff) &
          (log2FoldChange < log2FC.cutoff[1]) ~ 'down',
        TRUE ~ 'none'
      )) %>%
      arrange(factor(group, level = c('up', 'down', 'none')),
              desc(baseMean), # 表达量由高到低排序
              desc(abs(log2FoldChange))) # 差异表达倍数由高到低排序
    plot_data$group <-
      factor(plot_data$group, level = c('up', 'down', 'none'))
    show_subset1 <- plot_data %>%
      filter(group %in% c('up')) %>%
      arrange(desc(abs(log2FoldChange))) %>%
      head(show.num)
    show_subset2 <- plot_data %>%
      filter(group %in% c('down')) %>%
      arrange(desc(abs(log2FoldChange))) %>%
      head(show.num)
    show_subset <- rbind.data.frame(show_subset1, show_subset2)
    p1 <-
      ggplot(data = plot_data, aes(
        x = log2FoldChange,
        y = -log10(padj),
        color = group
      )) + geom_point(alpha = 0.8, size = 1.5) +
      scale_color_manual(values = c(color.pal, "grey")) +
      labs(x = "log2 (fold change)", y = "-log10 (q-value)") +
      theme(plot.title = element_text(hjust = 0.4)) +
      geom_hline(
        yintercept = -log10(padj.cutoff),
        lty = 4,
        lwd = 0.6,
        alpha = 0.8,
        color = '#de4307'
      ) +
      geom_vline(
        xintercept = log2FC.cutoff ,
        lty = 4,
        lwd = 0.6,
        alpha = 0.8,
        color = '#f29c2b'
      ) +
      theme_bw() +
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 12, face = 'plain'),
        legend.title = element_text(size = 12, face = 'bold')
      ) +
      scale_x_continuous(limits = c(
        min(plot_data$log2FoldChange),
        max(plot_data$log2FoldChange)
      )) +
      geom_text_repel(
        data = show_subset,
        aes(label = gene),
        col = "black",
        alpha = 0.8
      )
    return(list(
      deseq2.sig = plot_data,
      volcano.plot = p1
    ))
  }

#' GO and KEGG Pathway Enrichment Analysis
#'
#' @param sig.genes A charater vector of diffrentially expressed genes
#' @param GO.db Database for GO enrichment. DEFULAT 'org.Hs.eg.db' for homo sapiens
#' @param KEGG.db Database for KEGG enrichment. set to "internal" to use internal data. DEFULAT 'hsa' for homo sapiens
#' @param p.cutoff pvalue cutoff for significantly enriched terms. DEFAULT 0.05
#' @param q.cutoff qvalue cutoff for significantly enriched terms. DEFAULT 0.05
#' @param dot.color character vector, color hex in order
#'
#' @return A list containing 2 tables and 2 plots
#' @import org.Hs.eg.db
#' @import clusterProfiler
#' @import enrichplot
#' @export
#'
#' @examples
enrich.GO.KEGG <- function(sig.genes, GO.db = 'org.Hs.eg.db', KEGG.db = c('hsa', 'internal')[1],p.cutoff = 0.05, q.cutoff = 0.05, dot.color = c("#C34A36","#008BC8")) {
  require(org.Hs.eg.db)
  require(clusterProfiler)
  require(enrichplot)
  gene <- bitr(sig.genes,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO.db)
  if (KEGG.db == 'internal') {
    tryCatch({
      library(KEGG.db)
      KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                       organism = 'hsa',
                       pvalueCutoff = p.cutoff,
                       qvalueCutoff = q.cutoff,
                       use_internal_data = TRUE
      )
    },error = function(e){
      print('Please run following command first...')
      cat('install.packages(system.file("extdata", "KEGG.db_1.0.tar.gz", package = "ProgMan"), type = "source")')
      stop()
    }
    )
  }else{
    KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                     organism = KEGG.db,
                     pvalueCutoff = p.cutoff,
                     qvalueCutoff = q.cutoff,
                     use_internal_data = FALSE
    )
  }
  GO<-enrichGO(gene$ENTREZID,#GO富集分析
               OrgDb = GO.db,
               keyType = "ENTREZID",#设定读取的gene ID类型
               ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
               pvalueCutoff = p.cutoff,#设定p值阈值
               qvalueCutoff = q.cutoff,#设定q值阈值
               readable = T)
  GO.df <- as.data.frame(GO)
  print(paste0(nrow(GO.df), ' GO terms were significantly enriched...'))
  p1 <- barplot(GO, showCategory = 5, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(face = 'bold'),
          axis.text = element_text(face = 'bold'),
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(face = 'bold', colour = '#F1592A', size = 16)
    )+
    scale_fill_gradient(low=dot.color[1],high=dot.color[2])
  KEGG.df <- as.data.frame(KEGG)
  print(paste0(nrow(KEGG.df), ' KEGG pathways were significantly enriched...'))
  if (is.na(unique(KEGG.df$Description)[1])) {
    ptw_name <- toTable(KEGGPATHID2NAME)
    KEGG.df <- KEGG.df %>%
      left_join(ptw_name, by = c('ID' = 'path_id')) %>%
      mutate(Description = path_name) %>%
      dplyr::select(-path_name)
    KEGG@result$Description[1:nrow(KEGG.df)] <- KEGG.df$Description
  }
  p2 <- dotplot(KEGG) + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.title = element_text(face = 'bold'),
                              axis.text = element_text(face = 'bold')
  )+
    scale_colour_gradient(low=dot.color[1],high=dot.color[2])
  return(list(GO.df = GO.df, KEGG.df= KEGG.df, GO.barplot = p1, KEGG.dotplot = p2))
}

#' GSEA Analysis and Ridgeplot
#'
#' @param deseqs.res A data frame of DESeq2 results
#' @param OrgDb Database for annotation. DEFULAT 'org.Hs.eg.db' for homo sapiens
#' @param gmt.file A gmt format file with a 'term' column and 'gene' column
#' @param kegg.db Mandatory when gmt.file=NULL. DEFULAT 'hsa' for homo sapiens
#' @param fill.color character vector, color hex in order
#' @param p.cutoff pvalue cutoff for significantly enriched terms. DEFAULT 0.05
#'
#' @return A list containing a GSEA result object, a GSEA result table and a ridgeplot
#' @import org.Hs.eg.db
#' @import clusterProfiler
#' @import enrichplot
#' @import ggridges
#' @export
#'
#' @examples
enrich.GSEA <- function(deseqs.res, OrgDb = "org.Hs.eg.db", gmt.file = NULL, kegg.db = 'hsa',fill.color = c("#C34A36","#008BC8"),p.cutoff = 0.05) {
  require(org.Hs.eg.db)
  require(clusterProfiler)
  require(enrichplot)
  require(ggridges)
  gsea_gene <-
    data.frame(gene = deseqs.res$gene, log2FC = deseqs.res$log2FoldChange)
  eid_gene <-
    bitr(
      gsea_gene$gene,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = OrgDb
    )
  gsea_gene1 <- eid_gene %>%
    dplyr::rename(gene = SYMBOL) %>%
    left_join(gsea_gene, by = 'gene')
  geneList <- gsea_gene1$log2FC
  names(geneList) <-  gsea_gene1$ENTREZID #使用转换好的ID
  genelist <- sort(geneList, decreasing = T) #从高到低排序
  if (is.null(gmt.file)) {
    gsea_res <- gseKEGG(genelist, organism = kegg.db, pvalueCutoff = p.cutoff)#GSEA富集分析
    gsea_df <-setReadable(gsea_res, OrgDb = OrgDb, keyType = 'ENTREZID') %>% as.data.frame()
  }else{
    genelist <- gsea_gene$log2FC
    names(genelist) <-  gsea_gene$gene
    genelist <- sort(genelist, decreasing = T)
    gmt <- read.gmt(gmt.file)
    # gmt$term <- lapply(gmt$term, function(x){
    #   str_split(x, pattern = '_')[[1]][-1] %>% paste(., collapse = ' ') %>% str_to_title()
    # }) %>% unlist()
    gsea_res <- GSEA(genelist, TERM2GENE = gmt, pvalueCutoff = p.cutoff)
    gsea_df <- as.data.frame(gsea_res)
  }
  p1 <- enrichplot::ridgeplot(gsea_res, showCategory = 20, orderBy = "NES") +
    scale_fill_gradient(low=fill.color[1],high=fill.color[2]) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(face = 'bold'),
          axis.text = element_text(face = 'bold'))+
    labs(x = "Enrichment Distribution")
  return(list(GSEA.obj = gsea_res, GSEA.df = gsea_df, GSEA.ridgeplot = p1))
}


#' GSEAplot for Specific Terms
#'
#' @param gsea.res A GSEA result object
#' @param plot.term A numeric vector of terms for plotting
#' @param title.name Name for plot title
#' @param color.pal character vector, color hex in order
#'
#' @return A gseaplot
#' @export
#'
#' @examples
plotutil.gseaplot <- function(gsea.res, plot.term = c(1,2,3), title.name = 'KEGG', color.pal = c("#E64B35FF", "#00A087FF", "#3C5488FF", "#7570B3", "#E7298A","#8f262b", "#81532f")) {
  if (length(plot.term)>5) {
    print('Not recommended to create gseaplot for more than 5 terms')
  }
  p1 <- enrichplot::gseaplot2(
    gsea.res,
    title = paste0('Enrichment of ', title.name, ' terms'),
    geneSetID = plot.term,
    pvalue_table = FALSE,
    subplots = 1:3,
    rel_heights = c(2,0.5,0.8),
    color = color.pal[1:length(plot.term)],
    ES_geom = 'line'
  )
  return(p1)
}

#' Batch Wilcoxon Test for Multiple Columns
#'
#' @param test.columns A charactor vector of column names
#' @param group_col A character of group column name
#' @param data A data frame containing 'group_col' and 'test.columns'
#' @param padj.method choose from c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")
#'
#' @return A data frame of wilcoxon test statistics
#' @export
#'
#' @examples
batch.wilcoxon.test <- function(test.columns, group_col, data, padj.method = 'fdr') {
  wilcox_comparison <- function(column, group_col, data) {
    # 分组后是否两组都有数据,这里代码主要针对于高低风险两组各种药物敏感性差异
    if (!all(table(data[[group_col]], !is.na(data[[column]])) [,'TRUE'])){
      return(data.frame(statistic = NA, p.value = NA, method = NA, alternative = NA))
    } else {
      fml_str <- paste0('`',column, '`~`', group_col,'`')
      result <- wilcox.test(as.formula(fml_str), data = data)
      # 整理结果
      tidy_result <- broom::tidy(result)
      return(tidy_result)
    }
  }
  # 批量应用函数进行 Wilcoxon 检验
  wilcox_results <- lapply(test.columns, wilcox_comparison, group_col = group_col, data = data)
  # 如果需要合并所有药物的检验结果到一个数据框中
  wilcox_results_combined <- wilcox_results %>%
    do.call(rbind.data.frame, .) %>%
    mutate(drug = test.columns) %>%
    dplyr::select(drug, everything()) %>%
    filter(!is.na(p.value))
  wilcox_results_combined$p.adjust <- p.adjust(wilcox_results_combined$p.value, method = padj.method)
  return(wilcox_results_combined)
}

