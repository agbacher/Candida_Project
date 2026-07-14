#######################################################################################
# Author: Philipp Hofmann
# 
# This script will reproduce heatmaps from Figure 5-F and Supplementary Figure 4-C t
# 
#######################################################################################


library(dplyr)
library(scales)
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggrastr)


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ -------------- (1) LOAD DATA FOR HEATMAP (HC vs. CD) ------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

# set directories
save_dir <- "./Figure5/metadata/"

#### Read data
heatdat_status_short <-read.csv(paste0(save_dir,"data_MatrixPlot_GutTh17CalbiRXTcells_DonorLocation_short.csv"),row.names = 1)

# Transpose the matrix to cluster columns
heatdat_status_short<-heatdat_status_short[,!colnames(heatdat_status_short) %in% c('IRF7')]

oldrow <- rownames(heatdat_status_short)
newrow <- gsub('_.*','',oldrow)
rownames(heatdat_status_short) <- newrow


#################################################
## -- CLUSTER ROW & COLUMN DATA FOR HEATMAP -- ##
#################################################

# Column clustering
d_col <- dist(t(heatdat_status_short))
hc_col <- hclust(d_col)
col_order <- hc_col$labels[hc_col$order]

# Convert to long format
heatmap_df <- reshape2::melt(as.matrix(heatdat_status_short))
colnames(heatmap_df) <- c("Row", "Column", "Value")

# Reorder factor levels based on clustering
heatmap_df$Column <- factor(heatmap_df$Column, levels = col_order)

### Custom Row order
row_order <- rev(c("HC3", "HC4", "HC1", "HC2", "CD1", "CD5", "CD3", "CD2","CD8", "CD7", "CD4", "CD9", "CD6","CD10"))
heatmap_df$Row <- factor(heatmap_df$Row, levels = rev(row_order))


################################################
## --- GENERATE ROW/COL-CLUSTERED HEATMAP --- ##
################################################

# Define a color ramp with more steps
my_palette <- colorRampPalette(c("#ffffff",'#ffffff','grey90','#DBD9EA', "#2488F0", "#000080"))(100) # --> Final Color

# Set text axis
x_txt_size = 6
y_txt_size = 6

## Generate Heatmap
ggplot(heatmap_df, aes(x = Column, y = Row, fill = Value)) +
  geom_tile(color = "#f0f0f0") +
  scale_fill_gradientn(
    colours = my_palette,
    values = scales::rescale(c(-0.3, 0, 0.3)),  # keeps color focus around 0
    limits = c(-0.25, 0.25),  # or whatever your full data range is
    oob = squish,
    na.value = "#f0f0f0"
  ) +
  # modify legend
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black",
                               ticks = T, 
                               label = T,
                               label.theme =  element_text(size=y_txt_size),
                               #direction = "horizontal",
                               title.position = "top",
                               title = 'Scaled Expression',
                               title.hjust = .5,
                               title.theme = element_text(angle = 0,size=y_txt_size)))+
  theme_minimal() +
  theme(
    legend.position = 'right',
    legend.key.size = unit(0.4, 'cm'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    #axis.text.x = element_text(angle = 90, hjust =1,vjust=0.5,color = "black",size=x_txt_size),
    axis.text.x = element_text(angle = 45, hjust = 1,vjust=1.2, color = "black",size=x_txt_size),  # Set axis text color to black
    axis.text.y = element_text(color = "black",size=y_txt_size),  # Set y-axis text color to black
    legend.text = element_text(color = "black"),  # Set legend text color to black
    panel.grid = element_blank()
    
  ) +
  labs(x = '', y = '')+ #
  # Add frame when geom_tile(color='#f0f0f0') or geom_tile(color=NA)
  geom_rect(
    aes(xmin = min(as.numeric(Column)) - 0.5,
        xmax = max(as.numeric(Column)) + 0.5,
        ymin = min(as.numeric(Row)) - 0.5,
        ymax = max(as.numeric(Row)) + 0.5),
    inherit.aes = FALSE,
    color = "black",
    fill = NA,
    linewidth = 0.1
  ) -> p_heatmap # 

p_heatmap

### SAVE PLOTS
plt_dir <- './Figure5/Plots/'
cowplot::ggsave2(p_heatmap, file = paste0(plt_dir,'heatmap_small_Gut_HCvsCD.pdf'), width = 4.8, height = 1.9)
cowplot::ggsave2(p_heatmap, file = paste0(plt_dir,'NEW_heatmap_small_Gut_HCvsCD.pdf'), width = 4.8, height = 1.9)


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ---------------- (2) GENERATE LONG HEATMAP PER DONOR  [long CalbiRXT-TH17] ------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

# set directories
save_dir <- "./Figure5/metadata/"

#### Read data
heatdat_dnr_long <- read.csv(paste0(save_dir,"data_MatrixPlot_GutTh17CalbiRXTcells_Dnr_long.csv"),row.names = 1)

identical(colnames(heatdat_dnr_long),colnames(heatdat_dnr_long_new))

### Check
long_old <- reshape2::melt(as.matrix(heatdat_dnr_long))
colnames(long_old) <- c("Sample", "Gene", "Value_Old")
long_new <- reshape2::melt(as.matrix(heatdat_dnr_long_new))
colnames(long_new) <- c("Sample", "Gene", "Value_New")

merge <- long_old %>% left_join(.,long_new,by=c('Sample','Gene')) %>% mutate(dif=round(abs(Value_Old),5)-round(abs(Value_New),5))




#################################################
## -- CLUSTER ROW & COLUMN DATA FOR HEATMAP -- ##
#################################################

# Column clustering
d_col <- dist(t(heatdat_dnr_long))
hc_col <- hclust(d_col)
col_order <- hc_col$labels[hc_col$order]

# Convert to long format
heatmap_df <- reshape2::melt(as.matrix(heatdat_dnr_long))
colnames(heatmap_df) <- c("Row", "Column", "Value")

# Reorder factor levels based on clustering
col_order_long <- col_order
heatmap_df$Column <- factor(heatmap_df$Column, levels = col_order_long)

### Custom Row order
row_order <- rev(c("HC1", "HC2", "HC3", "HC4", "CD1", "CD2", "CD3", "CD4","CD5", "CD6", "CD7", "CD8", "CD9","CD10"))
heatmap_df$Row <- factor(heatmap_df$Row, levels = row_order)

################################################
## --- GENERATE ROW/COL-CLUSTERED HEATMAP --- ##
################################################

# Define a color ramp with more steps
my_palette <- colorRampPalette(c("#ffffff",'#ffffff','grey90','#DBD9EA', "#2488F0", "#000080"))(100) # --> Final Color


# Set text axis
x_txt_size = 6
y_txt_size = 7


## Generate Heatmap
ggplot(heatmap_df, aes(x = Column, y = Row, fill = Value)) +
  geom_tile(color = "#f0f0f0") +
  scale_fill_gradientn(
    colours = my_palette,
    values = scales::rescale(c(-0.3, 0, 0.3)),  # keeps color focus around 0
    limits = c(-0.25, 0.25),  # or whatever your full data range is
    oob = squish,
    na.value = "#f0f0f0"
  ) +
  # modify legend
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black",
                               ticks = T, 
                               label = T,
                               label.theme =  element_text(size=7),
                               #direction = "horizontal",
                               title.position = "left",
                               title = 'Scaled Expression',
                               title.hjust = .5,
                               title.theme = element_text(angle = 90,size=7)))+
  theme_minimal() +
  theme(
    legend.position = 'right',
    legend.key.size = unit(0.4, 'cm'),
    plot.margin = unit(rep(0, 4), 'lines'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    plot.title = element_text(hjust = 0.5),
    axis.ticks.x =  element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1,vjust=1.2, color = "black",size=x_txt_size),  # Set axis text color to black
    axis.text.y = element_text(color = "black",size=y_txt_size),  # Set y-axis text color to black
    legend.text = element_text(color = "black"),  # Set legend text color to black
    panel.grid = element_blank()
    
  ) +
  labs(x = '', y = '',title = 'Calbicans-specific gut-derived Th17 cells')+ #
  # Add frame when geom_tile(color='#f0f0f0') or geom_tile(color=NA)
  geom_rect(
    aes(xmin = min(as.numeric(Column)) - 0.5,
        xmax = max(as.numeric(Column)) + 0.5,
        ymin = min(as.numeric(Row)) - 0.5,
       ymax = max(as.numeric(Row)) + 0.5),
    inherit.aes = FALSE,
    color = "black",
    fill = NA,
    linewidth = 0.1
  ) -> p_heatmap # 

p_heatmap


### SAVE PLOTS
plt_dir <- './Figure5/Plots/'
cowplot::ggsave2(p_heatmap, file = paste0(plt_dir,'supplementary_heatmap_Gut_HCandCD_Th17CalbiRXT.pdf'), width = 12, height = 2.2)

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------------- (3) GENERATE LONG HEATMAP PER DONOR  [long TOTAL-TH17] ----------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

# set directories
save_dir <- "./Figure5/metadata/"

#### Read data
heatdat_dnr_long_totTH17 <- read.csv(paste0(save_dir,"data_MatrixPlot_Gut_TotalTH17_Dnr_long.csv"),row.names = 1)



#################################################
## -- CLUSTER ROW & COLUMN DATA FOR HEATMAP -- ##
#################################################

# Convert to long format
heatmap_df <- reshape2::melt(as.matrix(heatdat_dnr_long_totTH17))
colnames(heatmap_df) <- c("Row", "Column", "Value")

# Reorder factor levels based on clustering (keep same column order as before)
heatmap_df$Column <- factor(heatmap_df$Column, levels = col_order_long)

### Custom Row order
row_order <- rev(c("HC1", "HC2", "HC3", "HC4", "CD1", "CD2", "CD3", "CD4","CD5", "CD6", "CD7", "CD8", "CD9","CD10"))
heatmap_df$Row <- factor(heatmap_df$Row, levels = row_order)

################################################
## --- GENERATE ROW/COL-CLUSTERED HEATMAP --- ##
################################################

# Define a color ramp with more steps
my_palette <- colorRampPalette(c("#ffffff",'#ffffff','grey90','#DBD9EA', "#2488F0", "#000080"))(100) # --> Final Color

# Set text axis
x_txt_size = 6
y_txt_size = 7


## Generate Heatmap
ggplot(heatmap_df, aes(x = Column, y = Row, fill = Value)) +
  #geom_tile(color = "black",size = 0.2) +
  #geom_tile(color = "black") +
  geom_tile(color = "#f0f0f0") +
  #geom_tile(color=NA) +
  #geom_tile_rast() +
  scale_fill_gradientn(
    colours = my_palette,
    #colors=pals::brewer.greys(100),
    values = scales::rescale(c(-0.3, 0, 0.3)),  # keeps color focus around 0
    #values = scales::rescale(c(-0.25, 0, 0.31)),
    limits = c(-0.25, 0.25),  # or whatever your full data range is
    #limits = c(-0.15, 0.15),
    oob = squish,
    na.value = "#f0f0f0"
  ) +
  # modify legend
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black",
                               ticks = T, 
                               label = T,
                               label.theme =  element_text(size=7),
                               #direction = "horizontal",
                               title.position = "left",
                               title = 'Scaled Expression',
                               title.hjust = .5,
                               title.theme = element_text(angle = 90,size=7)))+
  theme(
    legend.position = 'right',
    legend.key.size = unit(0.4, 'cm'),
    plot.margin = unit(rep(0, 4), 'lines'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    plot.title = element_text(hjust = 0.5),
    axis.ticks.x =  element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1,vjust=1, color = "black",size=x_txt_size),  # Set axis text color to black
    axis.text.y = element_text(color = "black",size=y_txt_size),  # Set y-axis text color to black
    legend.text = element_text(color = "black"),  # Set legend text color to black
    panel.grid = element_blank()
    
  ) +
  labs(x = '', y = '',title = 'Total gut-derived Th17 cells')+ #
  # Add frame when geom_tile(color='#f0f0f0') or geom_tile(color=NA)
  geom_rect(
    aes(xmin = min(as.numeric(Column)) - 0.5,
        xmax = max(as.numeric(Column)) + 0.5,
        ymin = min(as.numeric(Row)) - 0.5,
        ymax = max(as.numeric(Row)) + 0.5),
    inherit.aes = FALSE,
    color = "black",
    fill = NA,
    linewidth = 0.1
  ) -> p_heatmap # 

p_heatmap


### SAVE PLOTS
plt_dir <- './Figure5/Plots/'
cowplot::ggsave2(p_heatmap, file = paste0(plt_dir,'supplementary_heatmap_Gut_HCandCD_totalTH17.pdf'), width = 12, height = 2.2)


