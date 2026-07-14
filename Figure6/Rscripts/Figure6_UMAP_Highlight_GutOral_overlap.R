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
# ▓ ----------------- (1) LOADING PROCESSED DATA CALBIRXT scDATA ------------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓


### LOAD METADATA WITH VELOCITY RESULTS
meta_data <- read.csv('./Figure5/metadata/CD_gut_CalbiRXT_metadata.csv')


### GENERATE SEURAT OBJECT
seu_h5 <- Read10X_h5('./Figure5/metadata/CD_GUT_CalbiRXT_h5/CD_GUT_CalbiRXT.h5')
seu_meta <- meta_data %>% tibble::column_to_rownames('cellID') %>% as.data.frame()
seu_meta$cc_2aa_tcrdist <- NULL
seu_meta$cc_2aa_tcrdist_size <- NULL
seu <- CreateSeuratObject(counts = seu_h5,meta.data = meta_data)
seu <- NormalizeData(seu)

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ----------------- (2) LOADING PROCESSED VDJ scDATA ---------------------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

## LOAD VDJ_MASTER TABLE (TRA/TRB-aa sequences with  =<2aa mismatch in both, processed in scirpy using 'tcrdist' with cutoff=10)
airr_data <- read.csv('./Figure4/metadata/Revision_BOG_filtered_VDJoutput_identity_TCRdist_cutoff10_sameV.csv') %>% # airr_data
                dplyr::mutate(Barcode=stringr::str_split(X, "_") %>% purrr:::map_chr(., 1)) %>%
                dplyr::mutate(donorID=donor) %>%
                dplyr::mutate(Status=ifelse(grepl('HC',donorID),'HC','CD')) %>%
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
res_airr_data <- airr_data %>% 
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
cells_per_dataset_airr_data <- res_airr_data %>%
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

final_airr_data <- res_airr_data %>% 
                      left_join(.,cells_per_dataset_airr_data,by='cc_2aa_tcrdist') %>%
                      dplyr::mutate(overlap_2aa=case_when(dataset_overlap=='gut-blood'~'blood-gut',
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


## ADD CALBI-REACTIVITY COLUMN
datasets_calbiRXT <- c('blood','blood-gut','blood-gut-oral','blood-oral')
final_airr_data_mod <- final_airr_data %>% 
                       mutate(Calbi_RXT_TRATRB_2aaDist=ifelse(overlap_2aa %in% datasets_calbiRXT, 'True', 'False')) 


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ----------------- (3) PLOTTIN HIGHLIGHT UMAP EMBEDDING ------------------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓


## Generate metadata for adding to 
head(final_airr_data_mod$cell_id)
head(colnames(seu))

seu_CalbiRXT_2aaDist <- seu@meta.data %>% 
                        tibble::rownames_to_column('cell_id') %>% 
                        left_join(.,final_airr_data_mod,by='cell_id') %>% 
                        tibble::column_to_rownames('cell_id') %>% 
                        dplyr::select(overlap_2aa,Calbi_RXT_TRATRB_2aaDist)


## ADDING metadata  
seu <- AddMetaData(seu,metadata = seu_CalbiRXT_2aaDist)

# Create a data frame for plotting
plot_df <- seu@meta.data %>%
           dplyr::select(UMAP_1,UMAP_2,cluster,overlap_2aa)

## Generate Plot [2-aa Distance TRA-TRB]
ggplot() +
   ## First layer: overlap (blood-gut)
   geom_point_rast(data = plot_df %>% filter(overlap_2aa=='blood-gut'),
                  aes(x = UMAP_1, y = UMAP_2),
                  size=0.5,
                  color = 'grey90',
                  raster.dpi = 1000) +
   ## Second layer: overlap (blood-gut-oral)
   geom_point_rast(data = plot_df %>% filter(overlap_2aa=='blood-gut-oral'),
                   aes(x = UMAP_1, y = UMAP_2),
                   size=2,
                  color = "black",
                  raster.dpi = 1000) +
   BuenColors::pretty_plot(fontsize = 8) +
   labs(x = "UMAP-1", y = "UMAP-2") +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.key.size = unit(0.3, 'cm')) -> umap_GutBloodOral_Highlighted_2aaDist 
umap_GutBloodOral_Highlighted_2aaDist


#################################
# -------- SAVE PLOT  --------- #
#################################

plt_dir <- './Figure6/Plots/'
cowplot::ggsave2(umap_GutBloodOral_Highlighted_2aaDist, file = paste0(plt_dir,'/umap_GutBloodOral_Highlighted_CD_CalbiRXT_Th17cells_2aa.pdf'), width = 2.5, height = 2.5)

