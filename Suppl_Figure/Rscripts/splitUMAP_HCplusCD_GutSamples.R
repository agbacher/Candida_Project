#######################################################################################
# Author: Philipp Hofmann
# 
# This script will be used to reproduce Supplementary Figure 4 - Panel A/B
#
#######################################################################################

library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrastr)
library(cowplot)

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ -------------- 1) LOADING PROCESSED SINGLE-CELL DATA -------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

## Create SeuratObject
seu_h5 <- Read10X_h5('./Figure5/metadata/HCplusCD_GUT_h5/HCplusCD_GUT.h5')
CDplusHC_gut <- CreateSeuratObject(counts = seu_h5)
CDplusHC_gut <- NormalizeData(CDplusHC_gut)

## Add metadata
meta <- read.csv('./Figure5/metadata/HCplusCD_Gut_metadata.csv') %>%
        tibble::column_to_rownames('cellID')
CDplusHC_gut <- AddMetaData(CDplusHC_gut,meta)

## Add embedding
cell_embeddings <- meta %>% dplyr::select(UMAP_1,UMAP_2)
colnames(cell_embeddings) <- c('emb_1','emb_2')
CDplusHC_gut[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(cell_embeddings),assay = "RNA",key = "UMAP_")

## Test
final <- c(
  "Th1/Tem"="#1D77B4",
  "Th17"= '#CE1256', #"#B949C1",
  "Tcm"="#94C8DA",
  "Migrating"="#000000",
  "Treg"="#9170B9",
  "Cytotoxic"="#472D5A",
  "NK-like"="#5D6572",
  "Other"="grey",
  "LowQC"="floralwhite"
)

DimPlot(CDplusHC_gut,cols = final, group.by = 'celltype')

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------ 2) REMOVE LOW-QC ANNOTATED CELLS - RECALCULATE UMAP -------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
set.seed(1234)
## Setting cell identities for DEG-testing
Idents(CDplusHC_gut) <- 'celltype' 
table(Idents(CDplusHC_gut))

## Add more metadata
CDplusHC_gut$donor <- gsub('_.*','',colnames(CDplusHC_gut))
CDplusHC_gut$status <- ifelse(grepl('HC',CDplusHC_gut$donor),'HC','CD')


## Remove Low-QC cells
Idents(CDplusHC_gut) <- 'celltype'
CDplusHC_gut_LowQCrm <- subset(CDplusHC_gut,idents=c('LowQC'),invert=T)


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------- 3) GENERATE UMAP EMBEDDING SPLIT BY STATUS: COLON --------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓


## Set Colors
cols <- c(
  "Th1/Tem"="#1D77B4",
  "Th17"= '#CE1256', #"#B949C1",
  "Tcm"="#94C8DA",
  "Migrating"="#000000",
  "Treg"="#9170B9",
  "Cytotoxic"="#472D5A",
  "NK-like"="#5D6572",
  "Other"="grey70",
  "Tfh"="#0F8299"
)

## Get UMAP cell embeddings
emb_col <- Embeddings(CDplusHC_gut_LowQCrm,reduction = "umap") %>% 
           as.data.frame() %>% 
           tibble::rownames_to_column('cellbarcode') %>%
           dplyr::rename(umap_1=UMAP_1, umap_2=UMAP_2)

## Extract metadata for celltype annotation
met_col <- CDplusHC_gut_LowQCrm@meta.data %>% 
           dplyr::select(celltype,status) %>% 
           tibble::rownames_to_column('cellbarcode')

## Merge data for plotting
plt_df_col <- met_col %>% left_join(.,emb_col,by='cellbarcode')

## Identify centroids for plotting labels
centroids_col <- aggregate(cbind(umap_1, umap_2) ~ celltype, data = plt_df_col, FUN = mean)

## Set plotting order for celltypes (double-check with DimPlot()-function in Seurat)
cell_order <- c("Th17","Treg","Migrating","Cytotoxic","NK-like","Th1/Tem","Tcm","Tfh","Other")
plt_df_col$celltype <- factor(plt_df_col$celltype,levels=cell_order)

###### Generate UMAP Split-plot from integrated Gut Healthy & CD dataset
##  UMAP for Healthy
ggplot() +
  ## rasterize UMAP
  geom_point_rast(data = plt_df_col %>% filter(status=='HC'), 
                  aes(x = umap_1, y = umap_2, color = celltype), shape=16, size = 0.5, raster.dpi = 1000) +
  #geom_point(data = plt_df_col %>% filter(status=='HC'),aes(x = umap_1, y = umap_2, color = celltype), shape=16, size = 0.5) +
  scale_fill_manual(values = cols,name='Cell type') +  # Apply color scale
  scale_color_manual(values = cols,name='Cell type') +  # Apply color scale
  BuenColors::pretty_plot() + 
  theme_void()+
  labs(x = "UMAP-1",y = "UMAP-2",title='Healthy') +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,size = 8),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.key.size = unit(0.2, "cm")) +
  
  ##  UMAP for CD
  ggplot() +
  ## rasterize UMAP
  geom_point_rast(data = plt_df_col %>% filter(status=='CD'), 
                  aes(x = umap_1, y = umap_2, color = celltype), shape=16, size = 0.5, raster.dpi = 1000) +
  #geom_point(data = plt_df_col %>% filter(status=='CD'),aes(x = umap_1, y = umap_2, color = celltype), shape=16, size = 0.5) +
  scale_fill_manual(values = cols,name='Cell type') +  # Apply color scale
  scale_color_manual(values = cols,name='Cell type',
                     breaks = c("Th17","Treg","Migrating","Cytotoxic","NK-like","Th1/Tem","Tcm","Tfh","Other"),
                     guide = guide_legend(override.aes = list(size = 4))) +  # Apply color scale
  BuenColors::pretty_plot() +
  theme_void()+
  labs(x = "UMAP-1",y = "UMAP-2",title='CD') +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5,size = 8),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.key.size = unit(0.2, "cm")) -> p_splitUMAP


################## SAVING PLOTS
plt_dir <- './Suppl_Figure/Plots'
cowplot::ggsave2(p_splitUMAP, file = paste0(plt_dir,'supplementary_splitUMAP_Gut_HCvsCD_raster.pdf'), width = 5.5, height = 3)


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------- 4) GENERATE MARKER DOTPLOT INTEGRATED REFERENCE ----------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
library(BuenColors)

###### PLOT CD_CalbiTh17_TRATRB_2aaMis DATA
# Make order-vector for plotting
id_order <- c('Tcm','Migrating','Tfh','Th1/Tem','Th17','Cytotoxic','NK-like','Treg','Other')

CDplusHC_gut_LowQCrm$celltype <- factor(CDplusHC_gut_LowQCrm$celltype,levels = rev(id_order))
Idents(CDplusHC_gut_LowQCrm) <- 'celltype'

genes_gut <-c('TCF7','LEF1','CCR7','IL6R',                # Tcm
              'KLF2','SELL','S1PR1',                      # Migrating
              "TOX2","CXCR5","CXCL13","IL21",             # Tfh
              'CCL5','CXCR3','HOPX','TBX21',              #"Tem/Th1"
              'CCR6','CEBPD','IL4I1','IL17A','IL22',      #"Th17"
              'GZMK','GZMH','NKG7','PLEK',                #"Cytotoxic"
              'GZMA','CD160','CD8A','KLRK1',              #"NK"
              'CTLA4','FOXP3','TIGIT','IKZF2',            #"Treg"
              'ENTPD1','ZEB2','LAG3','IL2RA'              #"Other"
) 


## GENERATE MARKER DOTPLOT
DotPlot(object = CDplusHC_gut_LowQCrm, 
        features = genes_gut,
        #features = rev(genes_gut), # Use if you want to make plot long instead of wide
        col.min = -1,
        col.max = 1,
        dot.scale = 8,
        #col.min = 0, 
        #scale.by = "size", 
        dot.min = .01
        
) + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + # For Adding black outlines around dot
  RotatedAxis() +
  labs(x='',y='') +
  #scale_y_discrete(position = "right") + # Use if you want to add celltype labels on top instead of bottom axis
  #coord_flip() +                         # Use if you want to make plot long instead of wide
  BuenColors::pretty_plot(fontsize = 10) +
  scale_color_gradientn(colors=BuenColors::jdb_palette("ocean_earth")) +
  #scale_color_gradient2(low = "grey90", high = "black",
  #                      limits = c(0, 1.5),
  #                      breaks = c(0.0, 0.5, 1.0, 1.5)
  #)  +
  guides(color = guide_colorbar(frame.colour = "black",
                                #ticks.colour = "black",
                                ticks.colour = NA,
                                ticks = T, 
                                label = T,
                                label.theme =  element_text(size=7),
                                title.position = "left",
                                title = 'Scaled Expression',
                                title.hjust = .5,
                                title.theme = element_text(angle = 90,size=8)),
         size = guide_legend(title.position = "left",
                             title = 'Expression (%)',
                             label.theme =  element_text(size=7),
                             title.hjust = .5,
                             title.theme = element_text(angle = 90,size=8),
                             override.aes=list(shape=21, colour="black", fill="black") # For Adding black outlines around dot
         )
         
  ) +
  theme(legend.position = "right",
        legend.key.size = unit(0.4, 'cm'),
        legend.spacing.y = unit(.1, "cm"),
        legend.margin = margin(t = 0, r = 10, b = 0, l = 10),  # Adjust margins if needed
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0,size = 10), # Use if you want to add celltype labels on top instead of bottom axis
        axis.title = element_blank()
  ) -> p_dotplotMarker
p_dotplotMarker


################## SAVING PLOTS
plt_dir <- './Suppl_Figure/Plots'
cowplot::ggsave2(p_dotplotMarker, file = paste0(plt_dir,'/supplementary_markerDotPlot_Gut_HCvsCD.pdf'), width = 9.5, height = 2.8)









