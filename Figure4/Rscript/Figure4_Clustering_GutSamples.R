#######################################################################################
# Date:   16.09.2024 
# Author: Philipp Hofmann
# 
# This Rscript will be used to analyse GEX data for GUT and ORAL tissue
# in Healthy Controls!
#
# We will generate UMAPs with cell annotation and densities
# 
# NOTE: inspired by Rpackages like 'Nebulosa' or 'schex'
#
#######################################################################################

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(BuenColors)
library(ggpubr)
library(ggrepel)
library(ggalt)
library(easybio)
#library(schex)
library(ggrastr)


###############################################################
# ---- 1) LOADING PROCESSED SINGLE-CELL DATA ---------------- #
###############################################################

## Create SeuratObject
seu_h5 <- Read10X_h5('./Figure4/metadata/HC_GUT_h5/HC_GUT.h5')
hc_gut <- CreateSeuratObject(counts = seu_h5)
hc_gut <- NormalizeData(hc_gut)


## Add metadata
meta <- read.csv('./Figure4/metadata/HC_gut_metadata.csv') %>%
        tibble::column_to_rownames('cellID')
hc_gut <- AddMetaData(hc_gut,meta)

## Add embedding
cell_embeddings <- meta %>% dplyr::select(UMAP_1,UMAP_2)
colnames(cell_embeddings) <- c('emb_1','emb_2')
hc_gut[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(cell_embeddings),assay = "RNA",key = "UMAP_")

## Test
final <- c(
  "Th1/Tem"="#1D77B4",
  "Th17"= '#CE1256', #"#B949C1",
  "Tcm"="#94C8DA",
  "Migrating"="#000000",
  "Treg"="#9170B9",
  "Cytotoxic"="#472D5A",
  "NK-like"="#5D6572"
)

DimPlot(hc_gut,cols = final, group.by = 'celltype')


######################################################
# ---- 3) LOAD GUT VDJ DATA (NEEEEW)) ---- #
######################################################

## LOAD VDJ_MASTER TABLE (TRA/TRB-aa sequences with  =<2aa mismatch in both, processed in scirpy using 'tcrdist' with cutoff=10)
airr_data_HC <- read.csv('./Figure4/metadata/Revision_BOG_filtered_VDJoutput_identity_TCRdist_cutoff10_sameV.csv') %>%
                dplyr::mutate(Barcode=stringr::str_split(X, "_") %>% purrr:::map_chr(., 1)) %>%
                dplyr::mutate(donorID=donor) %>%
                dplyr::select(-cc_aa_identity,-cc_aa_identity_size,-cc_aa_sameV_identity,-cc_aa_sameV_identity_size,-cc_2aa_sameV_tcrdist,-cc_2aa_sameV_tcrdist_size) %>%
                dplyr::mutate(cc_2aa_tcrdist=as.character(cc_2aa_tcrdist))

# Function to extract values safely
extract_value <- function(text, organ) {
  match <- stringr::str_extract(text, paste0(organ, ":[0-9]+"))  # Extract "organ:number"
  value <- as.numeric(stringr::str_extract(match, "[0-9]+"))  # Extract the number
  ifelse(is.na(value), 0, value)  # Replace NA with 0 if the organ is missing
}

############### GENERATE RESULTS FOR (WITH DATABASE)
### ADD DATASETOVERLAP COLUMN AND NEW CELL IDs
res_airr_data_HC <- airr_data_HC %>% 
                    # Summarize TCRs that come from different stimulations into one signature 'Blood'
                    mutate(origin=case_when(grepl('Th1|Th17|Calbi',stimu_tissue)~'blood',TRUE~stimu_tissue)) %>% 
                    na.omit() %>%
                    mutate(TRA_TRB_aa=paste0(cdr3_a_aa,'_',cdr3_b_aa)) %>%
                    dplyr::mutate(cell_id=paste0(donorID,'_',stimu_tissue,'_',Barcode))%>%
                    group_by(cc_2aa_tcrdist) %>%
                    mutate(dataset_overlap = paste(unique(origin), collapse='-')) %>%
                    ungroup() %>%
                    mutate(hyphen_count = sapply(strsplit(dataset_overlap, "-"), length)) %>%
                    dplyr::rename(n_overlap=hyphen_count) %>%   
                    tibble::column_to_rownames('X') 



## Calculate the number of cells from each dataset within clonotypes that overlap in Gut/Oral/Blood
cells_per_dataset_airr_data_HC <- res_airr_data_HC %>%
                                  na.omit() %>%
                                  # Compute cells per origin before collapsing
                                  group_by(cc_2aa_tcrdist, origin) %>%
                                  mutate(cells_per_origin_count = n()) %>%  # Temp column to store counts
                                  ungroup() %>%
                                  # Summarize into formatted string
                                  group_by(cc_2aa_tcrdist) %>%
                                  summarise(
                                    cells_per_origin = paste(unique(paste0(origin, ":", cells_per_origin_count)), collapse = "|"),
                                    cc_2aa_TRATRBA_tcrdist = dplyr::first(cc_2aa_tcrdist),  # Keep unique clone ID
                                    .groups = "drop"
                                  )  %>%
                                  dplyr::select(cc_2aa_tcrdist,cells_per_origin) 

final_airr_data_HC <- res_airr_data_HC %>% 
                      left_join(.,cells_per_dataset_airr_data_HC,by='cc_2aa_tcrdist') %>%
                      dplyr::mutate(dataset_overlap_mod=case_when(dataset_overlap=='gut-blood'~'blood-gut',
                                                                dataset_overlap=='oral-gut'~'gut-oral',
                                                                dataset_overlap=='gut-blood-oral' | dataset_overlap=='oral-gut-blood' | dataset_overlap=='gut-oral-blood' |
                                                                  dataset_overlap=='oral-blood-gut' | dataset_overlap=='blood-oral-gut' ~'blood-gut-oral',
                                                                dataset_overlap=='oral-blood'~'blood-oral',
                                                                TRUE~dataset_overlap)) %>%
                        mutate(
                          Ncell_blood = extract_value(cells_per_origin, "blood"),
                          Ncell_oral = extract_value(cells_per_origin, "oral"),
                          Ncell_gut = extract_value(cells_per_origin, "gut"),
                        )



## SUBSET TO CELLS FROM SINGLE-CELL DATA
datasets_calbiRXT <- c('blood','blood-gut','blood-gut-oral','blood-oral')

## ADD CALBI-REACTIVITY COLUMN
final_airr_data_HC <- final_airr_data_HC %>% mutate(Calbi_RXT_TRATRB_2aaMis=ifelse(dataset_overlap_mod %in% datasets_calbiRXT,'TRUE','FALSE'))

## Generate metadata for adding to 
seurat_CalbiRXT_HC <- hc_gut@meta.data %>% 
                      tibble::rownames_to_column('cell_id') %>% 
                      left_join(.,final_airr_data_HC[,!names(final_airr_data_HC) %in% c("donor","sample_map","stimu_tissue")],by='cell_id') %>% 
                      tibble::column_to_rownames('cell_id') %>% 
                      dplyr::select(names(final_airr_data_HC[,c("cc_2aa_tcrdist","cc_2aa_tcrdist_size","dataset_overlap_mod","cells_per_origin","Calbi_RXT_TRATRB_2aaMis")]))

table(seurat_CalbiRXT_HC$Calbi_RXT_TRATRB_2aaMis)             
table(seurat_CalbiRXT_HC$dataset_overlap_mod) 

gut_calbiRXT_DF_HC <- seurat_CalbiRXT_HC %>% tibble::rownames_to_column('cell_id') %>% filter(Calbi_RXT_TRATRB_2aaMis=='TRUE') %>% pull('cell_id')


##############################################
# ----- 3) UMAP WITH HIGHLIGHTED CELLS  ---- #
##############################################

table(hc_gut$celltype)


##### Calculate densities
# Gut-parameters
#cells_oi <- gut_calbiRXT_DF
cells_oi <- gut_calbiRXT_DF_HC
harmony_obj <- hc_gut
umap_red_key <- "umap"
umap1_key <- "UMAP_1"
umap2_key <- "UMAP_2"
anno <- "celltype"


# Extract UMAP coordinates
umap_data <- Embeddings(hc_gut, as.character(umap_red_key)) %>%
             as.data.frame() %>%
             dplyr::rename(UMAP_1=!!quo_name(umap1_key),UMAP_2=!!quo_name(umap2_key))  %>%
             tibble::rownames_to_column("cell")

# Add metadata (e.g., cells of interest)
umap_data <- umap_data %>% 
             mutate(group = ifelse(cell %in% cells_oi, 1, 0))

# Caclulate density values
density <- Nebulosa:::wkde2d(umap_data$UMAP_1, umap_data$UMAP_2, umap_data$group)

# Get density
get_dens <- Nebulosa:::get_dens(umap_data[,c("UMAP_1","UMAP_2","group")],density,"wkde")

# Add density vector
umap_data$density <- get_dens
umap_data <- umap_data %>% left_join(.,harmony_obj@meta.data%>%dplyr::select(!!quo_name(anno))%>%tibble::rownames_to_column('cell'),by='cell') %>%
             dplyr::rename(Annotation=!!quo_name(anno))

# Calculate cluster centroids for unique annotation labels
centroids <- aggregate(cbind(UMAP_1, UMAP_2) ~ Annotation, data = umap_data, FUN = mean)

# 4. Create the plot
DimPlot(harmony_obj, reduction = as.character(umap_red_key),cells.highlight = cells_oi,cols.highlight = c("#FA78FA"),group.by = as.character(anno),
        cols=c('gainsboro'),label = T,label.box = T,repel = T,label.size = 1.5,shuffle = T) +
  ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = density), size = 0.8, alpha = 0.7) +  # Color points by density
  scale_color_viridis_c(option = "inferno") +  # Apply color scale
  BuenColors::pretty_plot()+BuenColors::L_border()+
  theme(legend.position = "right")

## Manual Plot UMAP ORAL or Gut
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point_rast(data = subset(umap_data, as.character(group) == "0"),
                  aes(color = as.character(group)), 
                  size = 0.8, fill = "gainsboro", stroke = 0.4,
                  raster.dpi = 1000) +
  geom_point_rast(data = subset(umap_data, as.character(group) == "1"), 
             aes(color = as.character(group)), size = 0.8, shape = 21, 
             fill = "#797BA8", # Colon_color
             stroke = 0.4, 
             #color = "#AA1500", # Oral_color
             color = "#434579",
             raster.dpi = 1000) +  # Raster
  scale_color_manual(values = c("0" = "gainsboro", 
                                "1" = "#434579"
                                )) +  # Customize colors
  theme_minimal() +
  labs(#title = "Calbi-rxt cells (overlap blood and gut) in oral", 
       x = "UMAP-1", 
       y = "UMAP-2") +
  BuenColors::pretty_plot(fontsize = 8) +
  NoLegend() +
  theme(plot.background = element_rect(fill = "transparent", colour = NA)) -> p_gut_high

ggplot(umap_data[order(umap_data$density), ], aes(x = UMAP_1, y = UMAP_2)) + # Arrange UMAP-corrdinates by density to visibility better (especially important when small)
  #geom_point(aes(color = density), size = 0.8, shape=16,alpha=0.8) +  # Color points by density
  geom_point_rast(aes(color = density), size = 0.8, shape=16,alpha=0.8,raster.dpi = 1000) +
  #scale_color_viridis_c(option = "inferno") +  # Apply color scale
  scale_color_gradientn(colors = BuenColors::jdb_palette("solar_extra"),
                        breaks=c(min(umap_data$density),.00025,.0005,.00075,max(umap_data$density)),
                        labels=c('Lo','','','','Hi')
                        ) +
  BuenColors::pretty_plot(fontsize = 8) + labs(x = "UMAP-1",y = "UMAP-2") +
  guides(colour = guide_colorbar(frame.colour = "NA", ticks.colour = "NA",
                                  alpha = 1, 
                                  #barwidth = 0.5,  # Adjust bar width
                                  #barheight = 4    # Adjust bar height
  )) +
  theme(legend.position = "none",
        legend.title.position = "left",
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.title = element_text(angle = 90, hjust = 0.5),
        legend.key.size = unit(0.2, "cm")) -> p_gut_dense

ggplot(umap_data[order(umap_data$density), ], aes(x = UMAP_1, y = UMAP_2)) + # Arrange UMAP-corrdinates by density to visibility better (especially important when small)
  #geom_point(aes(color = density), size = 0.8, shape=16,alpha=0.8) +  # Color points by density
  geom_point_rast(aes(color = density), size = 0.8, shape=16,alpha=0.8,raster.dpi = 1000) +
  #scale_color_viridis_c(option = "inferno") +  # Apply color scale
  scale_color_gradientn(colors = BuenColors::jdb_palette("solar_extra"),
                        breaks=c(min(umap_data$density),.00025,.0005,.00075,max(umap_data$density)),
                        labels=c('Lo','','','','Hi')
  ) +
  BuenColors::pretty_plot(fontsize = 8) + labs(x = "UMAP-1",y = "UMAP-2") +
  guides(colour = guide_colorbar(frame.colour = "NA", ticks.colour = "NA",
                                 alpha = 1, 
                                 #barwidth = 0.5,  # Adjust bar width
                                 #barheight = 4    # Adjust bar height
  )) +
  theme(legend.position = "right",
        legend.title.position = "left",
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.title = element_text(angle = 90, hjust = 0.5),
        legend.key.size = unit(0.2, "cm")) -> p_gut_dense_legend


legend <- cowplot::get_legend(p_gut_dense_legend)


################## SAVING PLOTS
plt_dir <- './Figure4/Plots/'
cowplot::ggsave2(p_gut_high, file = paste0(plt_dir,'umap_Gut_Hightlight_HC_raster.pdf'), width = 2.5, height = 2.5)
cowplot::ggsave2(p_gut_dense, file = paste0(plt_dir,'umap_Gut_Density_HC_raster.pdf'), width = 2.5, height = 2.5)
cowplot::ggsave2(legend, file = paste0(plt_dir,'density_Gut_legend.pdf'), width = 2.5, height = 2.5)


###########################################
## CELL NUMBER HIST PER DONOR
###########################################

final <- c(
  "Th1/Tem"="#1D77B4",
  "Th17"= '#CE1256', #"#B949C1",
  "Tcm"="#94C8DA",
  "Migrating"="#000000",
  "Treg"="#9170B9",
  "Cytotoxic"="#472D5A",
  "NK-like"="#5D6572"
)

dnrID_meta <- hc_gut@meta.data %>% 
              dplyr::mutate(donorID=stringr::str_split(rownames(.), "_") %>% purrr:::map_chr(., 1)) %>%
              dplyr::select(donorID)
hc_gut <- AddMetaData(hc_gut,metadata=dnrID_meta)


SCpubr::do_BarPlot(hc_gut, 
                   group.by = "celltype",
                   split.by = "donorID",
                   plot.title = "",
                   position = "fill",
                   flip = FALSE,
                   return_data=T) -> p_GutBar

### Manual plot
dat_bar_gut <- p_GutBar$Data

### Manual order for plot
dat_bar_gut$celltype <- factor(dat_bar_gut$celltype, 
                                 levels = c("Th17","Th1/Tem", "Tcm", "Migrating","Cytotoxic", "Treg", "NK-like"))


## GUT
ggplot(dat_bar_gut, aes(x = donorID, y = freq*100, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "COLON",x='', y = "Cell type (%)", fill = "Cell Type") +
  #scale_fill_manual(values=v1) +
  #scale_fill_manual(values=v2) +
  #scale_fill_manual(values=v3) +
  #scale_fill_manual(values=v4) +
  #scale_fill_manual(values=orig_mod) +
  scale_fill_manual(values=final) +
  BuenColors::pretty_plot(fontsize = 10) +
  BuenColors::L_border() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = 0.5, size = 10, vjust = 1.5),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.text.x = element_text(angle = 45, hjust = 1)) -> p_bar_gut
p_bar_gut

plt_dir <- './Figure4/Plots/'
cowplot::ggsave2(p_bar_gut, file = paste0(plt_dir,'hist_Gut_CellPrct_HC.pdf'), width = 3, height = 2.5)





####################################################
# ------------ 4) Generate DotPlot --------------- #
####################################################

# Make order-vector for plotting
id_order <- c('Th17','Th1/Tem','Treg','Cytotoxic','NK-like','Tcm','Migrating')


hc_gut$celltype <- factor(hc_gut$celltype,levels = rev(id_order))
Idents(hc_gut) <- 'celltype'

genes_gut <-c('CCR6','IL4I1','IL17A','IL22',     # Th17
              "CXCR3","ITGA1","HOPX","EGR1",     # Th1/Tem
              'CTLA4','TIGIT','FOXP3','IKZF2', # Treg
               'GZMK','CCL4','GZMH','NKG7',       # Cytotoxic
              'KLRD1','XCL2','XCL1','KLRC1',     # NK
              'TCF7','LEF1','CCR7','IL6R',      # Tcm
              'LTB','KLF2','S1PR1','SELL'     # Migrating
              ) 

DotPlot(object = hc_gut, 
        features = genes_gut,
        #idents = id_order,
        col.min = 0, 
        scale.by = "size", 
        dot.min = .01) + 
  RotatedAxis() +
  labs(x='',y='') +
  #coord_flip() + 
  BuenColors::pretty_plot(fontsize = 10) +
  scale_color_gradient2(low = "grey90", high = "black")  +
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
  theme(legend.key.size = unit(0.4, 'cm'),
        axis.title = element_blank(),
        legend.position = "right",
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  ) -> p_dotplotMarker_gut
p_dotplotMarker_gut


#################################
# -------- SAVE PLOT  --------- #
#################################

plt_dir <- './Figure4/Plots/'
cowplot::ggsave2(p_dotplotMarker_gut, file = paste0(plt_dir,'/gut_cellMarker_DotPlot.pdf'), width = 7, height = 2.5)

