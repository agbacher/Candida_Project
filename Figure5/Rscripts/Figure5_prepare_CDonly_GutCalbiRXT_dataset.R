#######################################################################################
# Author: Philipp Hofmann
# 
# This script will be used to reproduce Figure 5 - Panel D (Volcano Plot - GSEA)
#
#######################################################################################


library(dplyr)
library(Seurat)
library(ggplot2)


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 
# ▓ ------------------------------- (1) LOADING PROCESSED SINGLE-CELL DATA ----------------------------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 

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


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 
# ▓ --------------------------------- (2) LOADING PROCESSED VDJ DATA ----------------------------------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 


## LOAD VDJ_MASTER TABLE (TRA/TRB-aa sequences with  =<2aa mismatch in both, processed in scirpy using 'tcrdist' with cutoff=10)
airr_data_new <- read.csv('./Figure4/metadata/Revision_BOG_filtered_VDJoutput_identity_TCRdist_cutoff10_sameV.csv') %>% # airr_data
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
res_airr_data_new <- airr_data_new %>% 
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
cells_per_dataset_airr_data_new <- res_airr_data_new %>%
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

final_airr_data_new <- res_airr_data_new %>% 
                        left_join(.,cells_per_dataset_airr_data_new,by='cc_2aa_tcrdist') %>%
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
final_airr_data_new <- final_airr_data_new %>% mutate(Calbi_RXT_TRATRB_2aaMis=ifelse(dataset_overlap_mod %in% datasets_calbiRXT,'TRUE','FALSE'))

## Generate metadata for adding to 
seurat_CalbiRXT_new <- CDplusHC_gut@meta.data %>% 
                       tibble::rownames_to_column('cell_id') %>% 
                       left_join(.,final_airr_data_new[,!names(final_airr_data_new) %in% c("donor","sample_map","stimu_tissue")],by='cell_id') %>% 
                       tibble::column_to_rownames('cell_id') %>% 
                       dplyr::select(names(final_airr_data_new[,c("cc_2aa_tcrdist","cc_2aa_tcrdist_size","dataset_overlap_mod","cells_per_origin","Calbi_RXT_TRATRB_2aaMis")]))

table(seurat_CalbiRXT_new$Calbi_RXT_TRATRB_2aaMis)             
table(seurat_CalbiRXT_new$dataset_overlap_mod) 

gut_calbiRXT_DF_new <- seurat_CalbiRXT_new %>% tibble::rownames_to_column('cell_id') %>% filter(Calbi_RXT_TRATRB_2aaMis=='TRUE') %>% pull('cell_id')


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 
# ▓ ------------------------------ (3) ADD CALBI-RXT TO SEURAT METADATA -------------------------------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 


## Modify 'Calbi_RXT_TRATRB_2aaMis' by converting NA to FALSE (NA mostly represent cells where no TCR-receptor was detected)
seurat_CalbiRXT_new <- seurat_CalbiRXT_new %>% 
                       mutate(Calbi_RXT_TRATRB_2aaMis=case_when(is.na(Calbi_RXT_TRATRB_2aaMis)~'FALSE',TRUE~Calbi_RXT_TRATRB_2aaMis)) %>%
                       dplyr::rename(Calbi_RXT_TRATRB_2aaMis_NEW=Calbi_RXT_TRATRB_2aaMis)
table(seurat_CalbiRXT_new$Calbi_RXT_TRATRB_2aaMis_NEW)    


## ADDING metadata  
CDplusHC_gut <- AddMetaData(CDplusHC_gut,metadata = seurat_CalbiRXT_new %>%dplyr::select('Calbi_RXT_TRATRB_2aaMis_NEW'))

CDplusHC_gut_CaRxt_2aa_new_meta <- CDplusHC_gut@meta.data %>% 
                                   dplyr::mutate(new_CaRxt_2aa=case_when(celltype=='Th17' & Calbi_RXT_TRATRB_2aaMis_NEW=='TRUE' ~ 'Th17_CalbiRXT',TRUE ~ celltype)) %>% 
                                   dplyr::select(new_CaRxt_2aa)

## Add to Seurat Object
CDplusHC_gut <- AddMetaData(CDplusHC_gut,metadata = CDplusHC_gut_CaRxt_2aa_new_meta)

## Check cell numbers
CDplusHC_gut$donor <- gsub('_.*','',colnames(CDplusHC_gut))
table(CDplusHC_gut[[c('donor','Calbi_RXT_TRATRB_2aaMis_NEW')]])
table(CDplusHC_gut[[c('donor','new_CaRxt_2aa')]])

### CHECK Calbicans-reactive CD vs HC (2aa Mismatch)
CDplusHC_gut$status <- ifelse(grepl('HC',CDplusHC_gut$donor),'HC','CD')
CDplusHC_gut$status_calbi_2aa_new <- paste0(CDplusHC_gut$status,'_',CDplusHC_gut$Calbi_RXT_TRATRB_2aaMis_NEW)

table(CDplusHC_gut[[c('donor','status_calbi_2aa_new')]])
table(CDplusHC_gut$status_calbi_2aa_new)

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 
# ▓ ------------------------- (4) SUBSET SEURAT TO CALBI-RXT CELLS (CD ONLY) --------------------------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 


## Subset To CD donors
Idents(CDplusHC_gut) <- 'status'
cd_gut <- subset(CDplusHC_gut,idents=c('CD'))

## Subset To Calbi-reactive cells
Idents(cd_gut) <- 'Calbi_RXT_TRATRB_2aaMis_NEW'
subset_gut <- subset(cd_gut,idents=c('TRUE'))
table(subset_gut$celltype)
sum(table(subset_gut$celltype))
subset_gut[["percent.mt"]] <- PercentageFeatureSet(object = subset_gut, pattern = "^MT-")
subset_gut <- NormalizeData(subset_gut, normalization.method = "LogNormalize", scale.factor = 10000)

## Load metadata
cdGut_calbirxt_meta <- read.csv('./Figure5/metadata/CD_gut_CalbiRXT_metadata.csv') %>%
                       tibble::column_to_rownames('cellID') %>%
                       dplyr::select(-donorID) 

subset_gut <- AddMetaData(subset_gut,cdGut_calbirxt_meta)

## Add embedding
cell_embeddings_Calbirxt <- cdGut_calbirxt_meta %>% dplyr::select(UMAP_1,UMAP_2)
subset_gut[["umap_CDcalbirxt"]] <- CreateDimReducObject(embeddings = as.matrix(cell_embeddings_Calbirxt),assay = "RNA",key = "UMAP_")

cols=c( 
  "Tcm"= "#94c8daff",
  "Tem"="#1d77b4ff",
  "Th17_1"="#d890a2ff",
  "Th17_2"="#9e0142ff"
)

## Visually inscpect
DimPlot(subset_gut, reduction = "umap_CDcalbirxt",group.by = 'cluster',cols = cols,label = T)

























