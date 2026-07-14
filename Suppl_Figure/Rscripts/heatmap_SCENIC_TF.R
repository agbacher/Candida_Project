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
# ▓ ------------ (2) PLOTTING TF HEATMAP (SEURATEXTEND) ------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

### IMPORTING pySCENIC LOOM FILE
scenicLoomPath_scope <- "./Figure5/metadata/CDgut_CalbiTh17_SCENIC_SCope_output_CPUmulti_10runs.loom"
scenic_SCOPEoutput <- ImportPyscenicLoom(scenicLoomPath_scope)

# Importing SCENIC Loom Files into Seurat
identical(colnames(seu_calb),rownames(scenic_SCOPEoutput$RegulonsAUC))

## Identifying Top Activated TFs in Each Cluster
tf_auc <- scenic_SCOPEoutput$RegulonsAUC
tf_gene_list <- scenic_SCOPEoutput$Regulons

mosaic::set.rseed(456)
top_n=30
tf_zscore <- CalcStats(tf_auc, f = seu_calb$cluster, order = "p", n = top_n, t = TRUE)
SeuratExtend::Heatmap(tf_zscore, lab_fill = "zscore")

rownames(tf_zscore) <- paste0(rownames(tf_zscore),'(+)')

## generate heatmap annotation
head(tf_zscore)
df_tf_zscore <- tf_zscore %>% mutate(regulon=rownames(.))
ha <- HeatmapAnnotation(
  Celltype = colnames(tf_zscore),
  col = list(Celltype = c(
    "Tcm"     = "#94c8daff",
    "Tem"     = "#1d77b4ff",
    "Th17_1"  = "#d890a2ff",
    "Th17_2"  = "#9e0142ff"
  )),
  gp = gpar(fontsize = 8, col = NA),
  border = T
)

# Example: genes to highlight
highlight_genes <- c('TBX21','BHLHE40','EOMES','LEF1','KLF2','JUNB','ETS1','IRF4','STAT4','IRF1','REL','FOSB','FOS','EGR1','JUN','TCF7','ATF3','NR4A1','BATF','EOMES','FOSL2','FOXO1','BACH2')
highlight_genes <- paste0(highlight_genes,'(+)')
# Row annotation for highlighting
ra = rowAnnotation(foo = anno_mark(at = which(rownames(tf_zscore) %in% highlight_genes),
                                   labels = rownames(tf_zscore)[rownames(tf_zscore)%in%highlight_genes],
                                   labels_gp = gpar(fontsize = 8, col='#3b3b3b')))

## Generate heatmap
my_palette_fn <- circlize::colorRamp2(
  breaks = c(min(tf_zscore), 0, max(tf_zscore)),
  colors = c("#0076B4", "#ffffff", "#D51506")
)

ComplexHeatmap::Heatmap(
  as.matrix(tf_zscore),
  name = "zscore",
  col = my_palette_fn,
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = T,
  row_names_gp = grid::gpar(fontsize = 8),
  show_column_names = F,
  top_annotation = ha,
  #right_annotation = ra,
  border = T,
  use_raster = T,
  #heatmap_legend_param = list(
  #  at = c("0","1"),
  #  labels = c("0","1"),
  #  title = "zscore"
  #)
) -> hm
hm

#################################
# -------- SAVE PLOT  --------- #
#################################
pdf("./Suppl_Figure/Plots/Heatmap_Regulons_Zscore_TF_perCellType.pdf",width = 5,height=12)
draw(hm,
     annotation_legend_side = "right",
     heatmap_legend_side = "right",
     column_title_gp=grid::gpar(fontsize=10))
dev.off()


