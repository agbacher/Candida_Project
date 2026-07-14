#######################################################################################
# Author: Philipp Hofmann
# 
# This script will be used to reproduce Figure 5 - Panel G/H/J/K 
#
#######################################################################################

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(BuenColors)
library(ggpubr)
library(ggrastr)
library(tidyr)
library(ggrepel)

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ --------------------- (1) LOADING VELOCITY-PROCESSED DATA --------------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓


### LOAD METADATA WITH VELOCITY RESULTS
velo_data <- read.csv('./Figure5/metadata/CD_gut_CalbiRXT_metadata.csv')


### GENERATE SEURAT OBJECT
seu_h5 <- Read10X_h5('./Figure5/metadata/CD_GUT_CalbiRXT_h5/CD_GUT_CalbiRXT.h5')
seu_meta <- velo_data %>% tibble::column_to_rownames('cellID') %>% as.data.frame()
seu <- CreateSeuratObject(counts = seu_h5,meta.data = seu_meta)
seu <- NormalizeData(seu)

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ----------------- (2) GENERATE PAGA-VELOCITY - ARROW-GRAPH -------------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓


##### GROOUP COLORS FOR PAGA
group_colors <- c('Tcm'="#94c8da",'Tem'="#1d77b4",'Th17_1'="#d890a2",'Th17_2'="#9e0142")
## Calculate group positions from UMAP coordinates 
group_pos <- velo_data %>%
             dplyr::select(cluster,UMAP_1,UMAP_2) %>%
             group_by(cluster) %>%
             summarize(x = mean(UMAP_1), y = mean(UMAP_2)) %>%
             dplyr::rename(group = cluster)

### FROM VELOCITY ANALYSIS
## Extract the transition confidence matrix
#transitions_confidence <- py_to_r(adata$uns$paga$transitions_confidence$todense())
#colnames(transitions_confidence) <- group_pos$group
#rownames(transitions_confidence) <- group_pos$group
transitions_confidence <- matrix(c(0.0000000, 0.0000000, 0.0000000, 0.0000000,
                                    0.2120701, 0.0000000, 0.0000000, 0.0000000,
                                    0.0000000, 0.4728388, 0.0000000, 0.0000000,
                                    0.0000000, 0.0000000, 0.2793102, 0.0000000),
                                  nrow = 4,
                                  byrow = TRUE)
colnames(transitions_confidence) <- group_pos$group
rownames(transitions_confidence) <- group_pos$group

## Edges / connectivities
#conn <- py_to_r(adata$uns$paga$connectivities$todense())
#colnames(conn) <- group_pos$group
#rownames(conn) <- group_pos$group
conn <- matrix(c(0.0000000, 0.4355561, 0.1527472, 0.1710236,
                 0.4355561, 0.0000000, 0.4101528, 0.0264817,
                 0.1527472, 0.4101528, 0.0000000, 0.6375946,
                 0.1710236, 0.0264817, 0.6375946, 0.0000000),
               nrow = 4,
               byrow = TRUE)
colnames(conn) <- group_pos$group
rownames(conn) <- group_pos$group

## Edges / connectivities
conn <- py_to_r(adata$uns$paga$connectivities$todense())
colnames(conn) <- group_pos$group
rownames(conn) <- group_pos$group

## Format edge list
edges <- as.data.frame(as.table(conn)) %>%
            filter(Freq > 0.3) %>%  # adjust threshold to filter weak edges
            dplyr::rename(source = Var1, target = Var2, weight = Freq) %>%
            left_join(group_pos, by = c("source" = "group")) %>%
            dplyr::rename(x_source = x, y_source = y) %>%
            left_join(group_pos, by = c("target" = "group")) %>%
            dplyr::rename(x_target = x, y_target = y)

edges_confidence <- as.data.frame(as.table(transitions_confidence)) %>%
                    filter(Freq > 0) %>%  # Filter weak transitions
                    dplyr::rename(target = Var1, source = Var2, confidence = Freq) %>%   # << flipped direction of arrows since they point into the oposite direction
                    left_join(group_pos, by = c("source" = "group")) %>%
                    dplyr::rename(x_source = x, y_source = y) %>%
                    left_join(group_pos, by = c("target" = "group")) %>%
                    dplyr::rename(x_target = x, y_target = y)


### Offset arrowheads
# Function to shift target point back from its center
offset_arrow <- function(df, r = 0.3) {
  dx <- df$x_target - df$x_source
  dy <- df$y_target - df$y_source
  d <- sqrt(dx^2 + dy^2)
  
  df$x_target_offset <- df$x_target - r * dx / d
  df$y_target_offset <- df$y_target - r * dy / d
  return(df)
}

edges_confidence <- offset_arrow(edges_confidence)

ggplot() +
  ## UMAP background points
  geom_point_rast(data = velo_data, 
                 aes(x = UMAP_1, y = UMAP_2, 
                    color = cluster,
                    fill = cluster
                     ), 
                 size = 0.8, 
                 shape= 16,
                 alpha = 0.4, 
                 raster.dpi = 1000) +
  ## PAGA edges
  #geom_segment(data = edges,
  #             aes(x = x_source, y = y_source, xend = x_target, yend = y_target), #alpha = weight),
  #             color = "black", size = 1) +
  ## PAGA edges + Arrows
  geom_segment(data = edges_confidence,
               #aes(x = x_source, y = y_source, xend = x_target, yend = y_target),
               aes(x = x_source, y = y_source,xend = x_target_offset, yend = y_target_offset),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "black", size = 0.8)+
  ## PAGA nodes
  geom_point_rast(data = group_pos, 
                  aes(x = x, y = y, fill = group,color=group
                      ), 
                  size = 4, 
                  stroke=.1, 
                  raster.dpi = 1000) +
  ## Labels
  geom_text_repel(data = group_pos,
                  aes(x = x, y = y, label = group),
                  fontface = "bold",
                  size = 2.5, color = "black") +
  ## Color scale from scvelo
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  ## Theme tweaks
  scale_alpha(range = c(0.3, 1)) +
  theme_minimal() +
  theme(legend.position='none',plot.background = element_rect(fill = "transparent", colour = NA),) +
  BuenColors::pretty_plot(fontsize = 8) + labs(x = "UMAP-1",y = "UMAP-2") -> p_paga_Arrow
  Seurat::NoLegend() -> p_paga_Arrow
p_paga_Arrow

################## SAVING PLOTS
plt_dir <- './Figure5/Plots/'
cowplot::ggsave2(p_paga_Arrow, file = paste0(plt_dir,'UMAP_PAGAgraph_withArrow.pdf'), width = 2.5, height = 2.5)



# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------- (3) PLOT SPLICE VS. PSEUDOTIME PAGA-VELOCITY  ------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# See github post: https://github.com/theislab/scvelo/discussions/883
# NOTE: scv.pl.scatter() uses First order moments (Mu and Ms) values to plot 'spliced' expression vs.
#       pseudotime. [Mu: first order moment "unspliced", Ms: first order moment "spliced"]

## Load Velocity spliced expression data
plot_df <- read.csv('./Figure5/metadata/VeloSpliced_Expression_CD_gut_CalbiRXT_data.csv')


# List of genes that overlap in GSEA-analysis from Calbi-RXT TH17 in Gut and show interesting velocity
genes_of_interest <- c('TCF7','LEF1','IFNG','GZMA','GZMB','IL21','MX1','IRF4','CD38','GBP4','GBP5','STAT1')  # Replace with your actual genes
plot_df$gene <- factor(plot_df$gene,levels=genes_of_interest)

# Now, plot the data using ggplot2
ggplot(plot_df, aes(x = pseudotime, y = expression)) +
  facet_wrap(~gene, scales = "free_y") +
  #geom_point() + 
  geom_point(aes(color = cluster)) +
  geom_smooth(aes(group  = gene),color = scales::alpha("black", 0.7),se = F, method = "loess", span = 0.6, size = 1.1) +
  scale_color_manual(values=c( "Tcm"= "#94c8daff","Tem"="#1d77b4ff","Th17_1"="#d890a2ff","Th17_2"="#9e0142ff"))+
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Customize facet labels (gene names)
    strip.background = element_rect(fill = "lightgray")  # Optional: background color for facet labels
  ) +
  labs(title = '',
       x = "Velocity Pseudotime",
       y = 'Expression')+
  BuenColors::pretty_plot()+theme(legend.key.size = unit(0.4, 'cm')) -> p


###### SAVE PLOTS
plt_dir <- './Figure5/Plots/'
cowplot::ggsave2(p, file = paste0(plt_dir,'gois_splicedRNAexpr_vs_pseudotime_trendline.pdf'), width = 3.5, height = 3)


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------- (4) PLOT SPLICE VS. PSEUDOTIME PAGA-VELOCITY (SEPARATE) ------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

############################################################# PLOT INDIVIDUAL
# List of genes you want to plot
genes_of_interest <- unique(plot_df$gene)

# Create a list of plots for each gene
plots <- lapply(genes_of_interest, function(gene) {
  gene_data <- plot_df %>% filter(gene == !!gene)
  
  ggplot(gene_data, aes(x = pseudotime, y = expression, color = cluster)) +
    geom_point(size=0.1) + # for figure
    #geom_point(size=3) + 
    scale_color_manual(values=c( 
      "Tcm"= "#94c8daff",
      "Tem"="#1d77b4ff",
      "Th17_1"="#d890a2ff",
      "Th17_2"="#9e0142ff"
    )) +
    theme_minimal() + 
    labs(
      title = gene,
      x = "Pseudotime",
      y ="Spliced Expression"
    ) +
    BuenColors::pretty_plot() +
    theme(#plot.title = element_text(hjust = 0.5,size = 4),# For figure (Center the title)
          plot.title = element_text(hjust = 0.5,size = 8),# For figure (Center the title)
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.position='none'
          )
}
)

# Combine all plots into one grid using patchwork
combined_plot <- patchwork::wrap_plots(plots, ncol = 6)  # Adjust ncol for number of columns

# Print the combined plot
print(combined_plot)

###### SAVE PLOTS
plt_dir <- './Figure5/Plots/'
cowplot::ggsave2(combined_plot, file = paste0(plt_dir,'gois_splicedRNAexpr_vs_pseudotime_minimal.pdf'), width = 5, height = 2)


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------- (5) PLOT VELOCITY-PSEUDOTIME UMAP-EMBEDDING (SEPARATE) ------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

ggplot() +
  # UMAP background points
  geom_point_rast(data = velo_data, 
                  aes(x = UMAP_1, y = UMAP_2, color = pseudotime),
                  size = 0.8, 
                  raster.dpi = 1000) +
  scale_color_gradientn(colors = BuenColors::jdb_palette("solar_extra"),
                        breaks=c(0,1)
                        )+
  guides(colour = guide_colorbar(frame.colour = "black",
                                 ticks.colour = "black",
                                 ticks = T, 
                                 label = T,
                                 label.theme =  element_text(size=7),
                                 #direction = "horizontal",
                                 title.position = "left",
                                 title = 'Pseudotime',
                                 title.hjust = .5,
                                 title.theme = element_text(angle = 90,size=8))
         ) +
  BuenColors::pretty_plot(fontsize = 8) + 
  labs(x = "UMAP-1",y = "UMAP-2") +
  theme(#legend.position = "none",
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.key.size = unit(.3, 'cm')) -> umap_PAGAvelocity_pseudotime_legend

## Get legend
legend <- get_legend(umap_PAGAvelocity_pseudotime_legend)

ggplot() +
  # UMAP background points
  geom_point_rast(data = velo_data, 
                  aes(x = UMAP_1, y = UMAP_2, color = pseudotime),
                  size = 0.8, 
                  raster.dpi = 1000) +
  scale_color_gradientn(colors = BuenColors::jdb_palette("solar_extra"),
                        breaks=c(0,1)
  )+
  guides(colour = guide_colorbar(frame.colour = "black",
                                 ticks.colour = "black",
                                 ticks = T, 
                                 label = T,
                                 label.theme =  element_text(size=7),
                                 #direction = "horizontal",
                                 title.position = "left",
                                 title = 'Pseudotime',
                                 title.hjust = .5,
                                 title.theme = element_text(angle = 90,size=8))
  ) +
  BuenColors::pretty_plot(fontsize = 8) + 
  labs(x = "UMAP-1",y = "UMAP-2") +
  theme(legend.position = "none",
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.key.size = unit(.3, 'cm')) -> umap_PAGAvelocity_pseudotime


###### SAVE PLOTS
plt_dir <- './Figure5/Plots/'
cowplot::ggsave2(umap_PAGAvelocity_pseudotime, file = paste0(plt_dir,'UMAP_VelocityPseudotime.pdf'), width = 2.5, height = 2.5)
cowplot::ggsave2(legend, file = paste0(plt_dir,'velocityPseudotime_Legend.pdf'), width = 2.5, height = 2.5)


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------- (6) GENERATE CLONESIZE-UMAP VISUALIZATION ------------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

## CloneSize
velo_meta_mod_vdj_size <- velo_data %>%
                          group_by(cc_2aa_tcrdist) %>%
                          mutate(clone_size = n()) %>%  # Generate clone_size column
                          ungroup()

# Generate Plot [without Border]
ggplot() +
  # UMAP background points
  geom_point(data = velo_meta_mod_vdj_size %>% arrange(clone_size), 
             aes(x = UMAP_1, y = UMAP_2, color = clone_size,size = clone_size),alpha = 0.8) +
  # Color scale 
  scale_color_gradientn(#colors = c("grey90", "#FFD978", "#3A97FF", "#000000"), # original choice
    #colors = c("grey90", "#FFD978", "#3361A5", "#000000"),  # COLOR1 
    colors = c("grey90", "#FDB31A", "#3361A5", "#000000"),  # COLOR2
    breaks = c(1, 5, 10, 16))+     
  # Control the range of point sizes
  scale_size_continuous(range = c(0.8, 3)) +
  # Control the legend for the colorbar
  guides(colour = guide_colorbar(frame.colour = "black",
                                 ticks.colour = "black",
                                 ticks = T, 
                                 label = T,
                                 label.theme =  element_text(size=7),
                                 #direction = "horizontal",
                                 title.position = "left",
                                 title = 'Clone size',
                                 title.hjust = .5,
                                 title.theme = element_text(angle = 90,size=8)),
         # Remove size scaling color bar 
         size='none'
  ) +
  BuenColors::pretty_plot(fontsize = 8) + 
  labs(x = "UMAP-1",y = "UMAP-2") +
  theme(#legend.position = "none",
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.key.size = unit(.3, 'cm'))  -> umap_cloneSize_woBorder_legend

# Then extract the colorbar-only legend
legend_umapCloneSize <- get_legend(umap_cloneSize_woBorder_legend)


# Generate Plot [without Border]
ggplot() +
  # UMAP background points
  geom_point(data = velo_meta_mod_vdj_size %>% arrange(clone_size), 
             aes(x = UMAP_1, y = UMAP_2, color = clone_size,size = clone_size),alpha = 0.8) +
  # Color scale 
  scale_color_gradientn(#colors = c("grey90", "#FFD978", "#3A97FF", "#000000"), # original choice
    #colors = c("grey90", "#FFD978", "#3361A5", "#000000"),  # COLOR1 
    colors = c("grey90", "#FDB31A", "#3361A5", "#000000"),  # COLOR2
    breaks = c(1, 5, 10, 16))+     
  # Control the range of point sizes
  scale_size_continuous(range = c(0.8, 3)) +
  # Control the legend for the colorbar
  guides(colour = guide_colorbar(frame.colour = "black",
                                 ticks.colour = "black",
                                 ticks = T, 
                                 label = T,
                                 label.theme =  element_text(size=7),
                                 #direction = "horizontal",
                                 title.position = "left",
                                 title = 'Clone size',
                                 title.hjust = .5,
                                 title.theme = element_text(angle = 90,size=8)),
         # Remove size scaling color bar 
         size='none'
  ) +
  BuenColors::pretty_plot(fontsize = 8) + 
  labs(x = "UMAP-1",y = "UMAP-2") +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.key.size = unit(.3, 'cm'))  -> umap_cloneSize_woBorder


#################################
# -------- SAVE PLOT  --------- #
#################################

plt_dir <- './Figure5/Plots/'
cowplot::ggsave2(umap_cloneSize_woBorder, file = paste0(plt_dir,'/umap_CloneSize_CD_CalbiRXT_Th17cells.pdf'), width = 2.5, height = 2.5)
cowplot::ggsave2(legend_umapCloneSize, file = paste0(plt_dir,'/umap_CloneSize_LegendColorBar.pdf'), width = 2.5, height = 2.5)


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------- (7) GENERATE MARKER-DOTPLOT FOR VELOCITY-DATASET ------------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓


## Find Markers for Plotting
select_genes <- rownames(seu)
genes2remove <- c(grep("^TR[ABGD][VJ]|^MT-|^RP[SL][0-9]+|^HB[^P]|^HTOC[0-9]+$", select_genes, value = TRUE)) # remove hashtag genes "^HTOC[0-9]+$"
genes2use <- setdiff(select_genes,genes2remove)

Idents(seu) <- 'cluster'
cd_velo_marker <- FindAllMarkers(seu,features = genes2use, only.pos = TRUE) %>% scCustomize::Add_Pct_Diff()
top15_cd_velo_marker <- cd_velo_marker %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)


###### PLOT DATA
# Make order-vector for plotting
id_order <- c('Tcm','Tem','Th17_1','Th17_2')

seu$cluster <- factor(seu$cluster,levels = rev(id_order))
Idents(seu) <- 'cluster'

# Selected Markers
genes_gut <-c('TCF7','CCR7','LEF1','SELL',        
              "IL7R","S100A4","LGALS3",    
              'TGFB1','SLC4A10','IL22',
              'RORA','CCR6','CXCR6','IL23R','IL4I1','IL17A','IL12RB2','TNFRSF18',"IFNG",'CSF2','IL21'
) 

DotPlot(object = seu, 
        features = genes_gut,
        #idents = id_order,
        col.min = 0, 
        scale.by = "size", 
        dot.min = .01
) + 
  RotatedAxis() +
  labs(x='',y='') +
  #coord_flip() + 
  BuenColors::pretty_plot(fontsize = 10) +
  scale_color_gradient2(low = "grey90", high = "black",
                        limits = c(0, 1.5),
                        breaks = c(0.0, 0.5, 1.0, 1.5)
  )  +
  guides(color = guide_colorbar(frame.colour = "black",
                                #ticks.colour = "black",
                                ticks.colour = NA,
                                ticks = T, 
                                label = T,
                                label.theme =  element_text(size=7),
                                title.position = "left",
                                title = 'Avg. Expression',
                                title.hjust = .5,
                                title.theme = element_text(angle = 90,size=8)),
         size = guide_legend(title.position = "left",
                             title = 'Expression (%)',
                             label.theme =  element_text(size=7),
                             title.hjust = .5,
                             title.theme = element_text(angle = 90,size=8))) +
  theme(legend.position = "right",
        legend.key.size = unit(0.4, 'cm'),
        legend.spacing.y = unit(.1, "cm"),
        legend.margin = margin(t = 0, r = 10, b = 0, l = 10),  # Adjust margins if needed
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title = element_blank()
  ) -> p_dotplotMarker_gut_velo
p_dotplotMarker_gut_velo


#################################
# -------- SAVE PLOT  --------- #
#################################

plt_dir <- './Figure5/Plots/'
cowplot::ggsave2(p_dotplotMarker_gut_velo, file = paste0(plt_dir,'/DotPlot_cellMarker_CD_CalbiRXT_Th17cells.pdf'), width = 5.8, height = 2.2)
