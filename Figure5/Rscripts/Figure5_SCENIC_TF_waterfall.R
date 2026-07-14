#######################################################################################
# Date:   03.12.2025 
# Author: Philipp Hofmann
# 
# This Rscript will be used to post-process LOOM-file output from pySCENIC.
#
#
#######################################################################################


#### Load Packages
# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(Seurat)

# For some of the plots:
#library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(SeuratExtend)
library(dplyr)

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------- (1) LOADING PRE-PROCESSED SEURAT-OBJECT ----------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
### LOAD METADATA WITH VELOCITY RESULTS
velo_data <- read.csv('./Figure5/metadata/CD_gut_CalbiRXT_metadata.csv')


### GENERATE SEURAT OBJECT
seu_h5 <- Read10X_h5('./Figure5/metadata/CD_GUT_CalbiRXT_h5/CD_GUT_CalbiRXT.h5')
seu_meta <- velo_data %>% tibble::column_to_rownames('cellID') %>% as.data.frame()
seu_calb <- CreateSeuratObject(counts = seu_h5,meta.data = seu_meta)
seu_calb <- NormalizeData(seu_calb)

## ADD UMAP CELL EMBEDDINGS
cell_embeddings <- velo_data %>% tibble::column_to_rownames('cellID') %>% dplyr::select(UMAP_1,UMAP_2)
colnames(cell_embeddings) <- c('emb_1','emb_2')
seu_calb[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(cell_embeddings),assay = "RNA",key = "UMAP_")


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ --------------- (2) LOADING pySCENIC OUTFILE (LOOM) ------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

### IMPORTING pySCENIC LOOM FILE
scenicLoomPath_scope <- "./Figure5/metadata/CDgut_CalbiTh17_SCENIC_SCope_output_CPUmulti_10runs.loom"

loom <- open_loom(scenicLoomPath_scope)
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
close_loom(loom)

## get cluster labels from seurat
cell_cl <- seu_calb@meta.data %>% dplyr::select(cluster)
cellCluster <- cell_cl[colnames(exprMat),, drop = FALSE]


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------------ (3) PLOT REGULON RSS-SCORES ------------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓


## Add scenic AUC-embedding clustered using UMAP
auc_emb <- embeddings$`SCENIC AUC UMAP`
colnames(auc_emb) <- c('auc_UMAP_1','auc_UMAP_2')
seu_calb[["auc_umap"]] <- CreateDimReducObject(embeddings = as.matrix(auc_emb),key = "AUC_UMAP_",assay = DefaultAssay(seu_calb))

## Add scenic AUC-embedding clustered using UMAP
auc_emb_tsne <- embeddings$`SCENIC AUC t-SNE`
colnames(auc_emb_tsne) <- c('auc_TSNE_1','auc_TSNE_2')
seu_calb[["auc_tsne"]] <- CreateDimReducObject(embeddings = as.matrix(auc_emb_tsne),key = "AUC_TSNE_",assay = DefaultAssay(seu_calb))


## Plots
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellCluster), cellCluster[,'cluster']) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

# plot:
options(repr.plot.width=8, repr.plot.height=10) # To set the figure size in Jupyter
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)
viewTable(topRegulators, options = list(pageLength = 10))

## RSS score
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellCluster[colnames(regulonAUC), 'cluster'])
rss_df <- reshape2::melt(rss)
colnames(rss_df) <- c("Regulon", "Cluster", "RSS")
head(rss_df)
# ---- PARAMETERS ----
top_n <- 20   # number of regulons to label per cluster

# ---- PREPARE DATA ----
rss_df2 <- rss_df %>%
  group_by(Cluster) %>%
  arrange(desc(RSS), .by_group = TRUE) %>%
  mutate(
    Rank = dplyr::row_number(),
    Label = ifelse(Rank <= top_n, as.character(Regulon), NA)
  ) %>%
  ungroup()
# ---- PLOT ----
ggplot(rss_df2, aes(x = Rank, y = RSS)) +
  geom_point(size = 0.8, color = ifelse(!is.na(rss_df2$Label),"#D51506","#0076B4")) +
  ggrepel::geom_text_repel(
    aes(label = Label),
    size = 3,
    color = "#D51506",
    segment.color = 'grey80',
    bg.colour = "white",
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.3,
    point.padding = 0.2,
    na.rm = TRUE
  ) + 
  facet_wrap(~ Cluster, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Ranked regulons",
    y = "Regulon specificity score (RSS)",
    title = "RSS per cluster"
  ) +
  BuenColors::pretty_plot() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 10)
  ) -> p_RSSscore
p_RSSscore

#################################
# -------- SAVE PLOT  --------- #
#################################
plt_dir <- './Figure5/Plots/'
cowplot::ggsave2(p_RSSscore, file = paste0(plt_dir,'regulon_RSSscore_perCluster.pdf'), width = 8, height = 8)


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ --------- (4) DIFFERENTIAL TF ANALYSIS (WATERFALL-PLOT) --------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# NOTE: This Analysis is inspired by the WaterfallPlot-function from SeuratExtend R package

set.seed(123)

## Manual p-values and t-scores [done for each regulon]
auc = getAUC(regulonAUC)
## Get cluster
cell_cl <- seu_calb@meta.data %>% dplyr::select(cluster)
cluster <- cell_cl[colnames(auc),, drop = FALSE] %>%
           dplyr::mutate(cell_id=rownames(.)) %>%
           dplyr::rename(Cluster=cluster)

## Set identities to test
ident1='Th17_2'
ident2='Th17_1'
cell.1 = cluster %>% filter(Cluster==ident1) %>% pull('cell_id')
cell.2 = cluster %>% filter(Cluster==ident2) %>% pull('cell_id')

# safe Wilcoxon returning BOTH statistic and p-value
safe_wilcox_full <- function(x, idx1, idx2) {
  x1 <- x[idx1]; x2 <- x[idx2]
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  if (length(x1) < 1 || length(x2) < 1) return(list(statistic = NA_real_, p.value = NA_real_))
  if (length(unique(c(x1, x2))) <= 1) return(list(statistic = 0, p.value = 1))
  out <- suppressWarnings(stats::wilcox.test(x1, x2, exact = FALSE))
  list(statistic = unname(out$statistic), p.value = out$p.value)
}

# One Wilcoxon per regulon
wt <- apply(auc, 1, function(x) safe_wilcox_full(x, cell.1, cell.2))

# Extract statistic + p
w_stat  <- sapply(wt, function(res) res$statistic)
p_value_w <- sapply(wt, function(res) res$p.value)

## assemble test df
wx_df <- data.frame(regulon=names(w_stat),wx_stat=unname(w_stat)) 
wx_p_df <- data.frame(regulon=names(p_value_w),pvalue_wx=unname(p_value_w))

stat_df <- data.frame(regulon=rownames(auc)) %>% 
            left_join(.,wx_df,by='regulon') %>%
            left_join(.,wx_p_df,by='regulon') 

stat_df$p_adj_wx <- p.adjust(stat_df$pvalue_wx, method = "BH")

## re-order
stat_df <- stat_df %>% dplyr::select(regulon,wx_stat,pvalue_wx,p_adj_wx)

## Calculating logFC values (AUC ranges from 0-1 and can be very close_toggle zero, pseudocount is added to stabilize ratios and avoids division by zero)
pc <- 0.01  # pseudocount for AUC data

logFC_e <- apply(auc, 1, function(x) {
  mean1 <- mean(x[cell.1], na.rm = TRUE)
  mean2 <- mean(x[cell.2], na.rm = TRUE)
  log((mean1 + pc) / (mean2 + pc))
})

logFC_2 <- apply(auc, 1, function(x) {
  mean1 <- mean(x[cell.1], na.rm = TRUE)
  mean2 <- mean(x[cell.2], na.rm = TRUE)
  log2((mean1 + pc) / (mean2 + pc))
})


logFCe_df <- data.frame(
  regulon = names(logFC_e),
  logFCe = unname(logFC_e)
)

log2FC_df <- data.frame(
  regulon = names(logFC_2),
  log2FC= unname(logFC_2)
)


stat_df_final <- stat_df %>%
                  left_join(logFCe_df, by = "regulon") %>%
                  left_join(log2FC_df, by = "regulon") %>%
                  mutate(signed_log10_padj_wx = -log10(p_adj_wx) * sign(logFCe))

# ---- WATERFALL Plot ----
stat_df_final_plt <- stat_df_final %>% 
                      filter(regulon %in% c('BHLHE40_(+)','TBX21_(+)','IRF4_(+)','IRF1_(+)','ATF3_(+)','JUNB_(+)','ETS1_(+)','FOS_(+)','LEF1_(+)')) %>%
                      dplyr::mutate(regulon=gsub('_','',regulon))

stat_df_final_plt$regulon <- factor(stat_df_final_plt$regulon,levels = stat_df_final_plt$regulon[order(stat_df_final_plt$log2FC, decreasing = TRUE)])

lim <- max(abs(stat_df_final_plt$signed_log10_padj_wx), na.rm = TRUE)

ggplot(stat_df_final_plt, 
       aes(x = regulon, y = log2FC, color = signed_log10_padj_wx)) +
  geom_segment(aes(xend = regulon, y = 0, yend = log2FC), linewidth = 0.5) +
  geom_hline(yintercept = 0,color = "grey90") +
  geom_point(size = 1) +
  #geom_point(aes(size = signed_log10_padj)) +
  scale_color_gradient2(
    low = "#0076B4", mid = "white", high = "#D51506",
    midpoint = 0,
    limits = c(-ceiling(lim), ceiling(lim)),
    oob = scales::squish
  ) +
  guides(color = guide_colorbar(frame.colour = "black",
                                ticks.colour = "black",
                                ticks = F, 
                                label = T,
                                label.theme =  element_text(size=4),
                                #direction = "horizontal",
                                title.position = "left",
                                title = 'Signed -log10(FDR)',
                                title.hjust = .5,
                                title.theme = element_text(angle = 90,size=4))) +
  labs(
    x = "",
    y = "log2(FC)",
    color = "signed(-log10(p))",
    title = paste0(ident1," vs. ",ident2)
  ) +
  ylim(-0.6,0.6) + 
  BuenColors::pretty_plot(fontsize = 4) + 
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.key.size = unit(0.2, 'cm'),
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    plot.title = element_blank()
  ) -> p_legend


legend_waterfall <- cowplot::get_legend(p_legend)


ggplot(stat_df_final_plt, 
       aes(x = regulon, y = log2FC, color = signed_log10_padj_wx)) +
  geom_segment(aes(xend = regulon, y = 0, yend = log2FC), linewidth = 0.4) +
  geom_hline(yintercept = 0,color = "grey90") +
  geom_point(size = 0.8) +
  #geom_point(aes(size = signed_log10_padj)) +
  scale_color_gradient2(
    low = "#0076B4", mid = "white", high = "#D51506",
    midpoint = 0,
    limits = c(-ceiling(lim), ceiling(lim)),
    oob = scales::squish
  ) +
  labs(
    x = "",
    y = "log2(FC)",
    color = "signed(-log10(p))",
    title = paste0(ident1," vs. ",ident2)
  ) +
  ylim(-0.6,0.6) + 
  BuenColors::pretty_plot(fontsize = 4) + 
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.position = 'none',
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5,size=3),
    plot.title = element_blank()
  ) -> p_waterfall
p_waterfall


#################################
# -------- SAVE PLOT  --------- #
#################################
plt_dir <- './Figure5/Plots/'
cowplot::ggsave2(p_waterfall, file = paste0(plt_dir,'waterfallPlot_Regulon_GutCalbiTCells_Th17-1_vs_Th17-2.pdf'), width = 0.6, height = 0.8)
cowplot::ggsave2(legend_waterfall, file = paste0(plt_dir,'legend_waterfallPlot.pdf'), width = 0.6, height = 0.8)


