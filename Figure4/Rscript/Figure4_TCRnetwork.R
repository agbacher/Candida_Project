##################################################################################################################################
## Author: Philipp Hofmann
## 
## This Rscript will be used to generate network, showing C.albicans-reactive
## clonotypes in Oral and Colon Biopsies. Similar to Figure 5d --> https://doi.org/10.1038/s41591-023-02556-5
## 
## We will generate a nodes and edges file. The nodes file will have 3x columns (cc_2aa_TRATRBA_tcrdist,size,overlap) and the
## edges file will have 5x columns ("cc_2aa_TRATRBA_tcrdist", "tissue", "size", "donor","diagnosis").
##
## 
####################################################################################################################################

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓--------------------------- 1) LOADING R-LIBRARIES  --------------------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

library(magrittr)
library(dplyr)
library(networkD3)
library(pathfindR)
library(forcats)
library(ggplot2)
library(igraph)
library(webshot)

# Function to extract values safely
extract_value <- function(text, organ) {
  match <- stringr::str_extract(text, paste0(organ, ":[0-9]+"))  # Extract "organ:number"
  value <- as.numeric(stringr::str_extract(match, "[0-9]+"))  # Extract the number
  ifelse(is.na(value), 0, value)  # Replace NA with 0 if the organ is missing
}


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ -------------------- 2) LOAD VDJ DATA (WITH 2AA MISMATCH) -------------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

## LOAD VDJ_MASTER TABLE (TRA/TRB-aa sequences with  =<2aa mismatch in both, processed in scirpy using 'tcrdist' with cutoff=10)
airr_data_HC <- read.csv('./Figure4/metadata/Revision_BOG_filtered_VDJoutput_identity_TCRdist_cutoff10_sameV.csv') %>% # airr_data
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
res_airr_data_HC <- airr_data_HC %>% 
                    filter(Status == 'HC') %>%
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

                              
## check how many clonotypes per tissue
length(unique(res_airr_data_HC$cc_2aa_tcrdist[res_airr_data_HC$origin=='blood' & res_airr_data_HC$cc_2aa_tcrdist_size >= 2])) # -> BLOOD (clone_size > 1): N = 2331
length(unique(res_airr_data_HC$cc_2aa_tcrdist[res_airr_data_HC$origin=='gut' & res_airr_data_HC$cc_2aa_tcrdist_size >= 2]))   # -> GUT (clone_size > 1):   N = 3706
length(unique(res_airr_data_HC$cc_2aa_tcrdist[res_airr_data_HC$origin=='oral' & res_airr_data_HC$cc_2aa_tcrdist_size >= 1]))  # -> ORAL (clone_size > 1):  N = 731
                              

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------------ 3) GENERATE EDGES/NODES FILE FOR NETWORK ------------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

## This version below seems to be correct 
edges_healthy_new <- final_airr_data_HC %>%
                      filter(Status == 'HC') %>%
                      dplyr::select(cc_2aa_tcrdist, 
                                    cc_2aa_tcrdist_size, 
                                    dataset_overlap_mod, 
                                    n_overlap) %>%
                      dplyr::rename(dataset_overlap = dataset_overlap_mod) %>%
                      tidyr::separate_rows(dataset_overlap, sep = "-") %>%
                      filter(
                        (dataset_overlap %in% c("gut", "blood") & cc_2aa_tcrdist_size > 1) |
                          (dataset_overlap == "oral")  # keep all for tissue3
                      ) %>%
                      mutate(tissue = dataset_overlap) %>%
                      dplyr::select(-dataset_overlap) %>%
                      mutate(across(everything(), as.character))  %>%
                      mutate(n_overlap=as.numeric(n_overlap))

## Check counts per clone
actual_counts_new <- edges_healthy_new %>%
                     count(cc_2aa_tcrdist, name = "actual_n_overlap")

# Join actual counts back to the main data
edges_join_new <- edges_healthy_new %>%
                  left_join(actual_counts_new, by = "cc_2aa_tcrdist")

# Fix the counts
edges_trimmed_new <- edges_join_new %>%
                    group_by(cc_2aa_tcrdist) %>%
                    distinct(tissue, .keep_all = TRUE) %>%
                    filter(row_number() <= as.numeric(first(n_overlap))) %>%
                    ungroup() %>%
                    dplyr::select(-actual_n_overlap) %>%
                    group_by(cc_2aa_tcrdist) %>%
                    mutate(dataset_overlap = paste(unique(tissue), collapse='-')) %>%
                    ungroup()



# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------- 4) PLOTTING THE CLONAL-OVERLAP PER TISSUE NETWORK -------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# -------------------------------------------
# - Network Properties from Edge/Node files -
# -------------------------------------------
# Source variable - cc_2aa_TRATRBA_tcrdist
# Target variable - tissue
# NodeID/Name     - cc_2aa_TRATRBA_tcrdist
# Group           - Group for coloring ("blood", "oral", "gut","Xreactive")


## Network color, link color -> #DCDEDD
customCols <- JS(
  'd3.scaleOrdinal()
  .domain(["oral", "blood", "gut", "gut-oral", "blood-oral", "blood-gut", "blood-gut-oral"])
  .range(["#703E7E", "#3E4481", "#8D8E8F", "#AC2D77", "#AC2B76", "#8787A8", "#D0125D"])' # Row corresponds to groups
)


net_color <- c('oral'='#703E7E',           # group 1
               'gut'='#3E4481',            # group 2
               'blood'='#8D8E8F',          # group 3
               'gut-oral'='#AC2D77',       # group 4
               'blood-oral'='#AC2B76',     # group 5
               'blood-gut'='#8787A8',      # group 6
               'blood-gut-oral'='#D0125D'  # group 7
               )  # Generate simple network for whole dataset (Crohn's + Healthy Donor)


networkData_new <- data.frame(src=edges_trimmed_new$cc_2aa_tcrdist, target=edges_trimmed_new$tissue)
simpleNetwork(networkData_new,linkColour = '#DCDEDD',nodeColour="firebrick",charge = -10,linkDistance = 80)

### USE igraph workaround
# Make network input data: Source=cc_2aa_TRATRBA_tcrdist, Target=antigen_species
networkData_new <- data.frame(src=edges_trimmed_new$cc_2aa_tcrdist, target=edges_trimmed_new$tissue)

# Make net graph object
net_igraph_new <- graph_from_data_frame(networkData_new, directed = F)

# Get grouping (clusters) from network
wc <- cluster_walktrap(net_igraph_new,steps = 4)
members <- membership(wc)

# Convert igraph-network to network3D compatible network
net_3D_new <- igraph_to_networkD3(net_igraph_new, group = members)

# Merge net_3D nodes with manual curated nodes files to get correct grouping (colors)
net3D_nodes_new <- net_3D_new$nodes %>% dplyr::rename(grp_net3D=group)
nodes_final_new <- net3D_nodes_new %>% 
                  left_join(.,edges_trimmed_new %>%
                              distinct(cc_2aa_tcrdist, .keep_all = TRUE) %>% 
                              dplyr::rename(name=cc_2aa_tcrdist,
                                            group=tissue),
                            by='name') %>% 
                  dplyr::select(-grp_net3D) %>% 
                  mutate(group=case_when(is.na(group) & name=='oral' ~ '1',
                                         is.na(group) & name=='gut' ~ '2',
                                         is.na(group) & name=='blood' ~ '3',
                                         TRUE ~ group))  %>% 
                  mutate(cc_2aa_tcrdist_size=case_when(name=='oral' ~ '1',
                                                               name=='gut' ~ '1',
                                                               name=='blood' ~ '1',
                                                               TRUE ~ cc_2aa_tcrdist_size))   %>% 
                  mutate(group=case_when(dataset_overlap=='oral' ~ '1',
                                         dataset_overlap=='gut' ~ '2',
                                         dataset_overlap=='blood' ~ '3',
                                         dataset_overlap=='gut-oral' ~ '4',
                                         dataset_overlap=='blood-oral' ~ '5',
                                         dataset_overlap=='blood-gut' ~ '6',
                                         dataset_overlap=='blood-gut-oral' ~ '7',
                                         TRUE ~ group)) %>%
                  mutate(dataset_overlap = case_when(is.na(dataset_overlap) & name %in% c("oral", "gut", "blood") ~ name,
                                                     TRUE ~ dataset_overlap),
                         dataset_overlap = as.character(dataset_overlap))

# Set group names to overlap so that colors are matched correctly
unique(nodes_final_new$group)


### Generate Custom Network
forceNetwork(Links = net_3D_new$links, Nodes = nodes_final_new, 
             Source = 'source', Target = 'target', 
             NodeID = 'name', #Group = 'group', 
             Group = 'dataset_overlap', 
             Nodesize = 'cc_2aa_tcrdist_size',
             fontSize=0,
             legend=TRUE,
             opacity = 0.9,
             charge = -100,
             linkColour = '#DCDEDD',
             colourScale = JS(customCols),
             fontFamily = 'sans',
             opacityNoHover = 1,
             zoom=T)


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ---------------- 5) EXPORT EDGES & NODES FILES FOR CYTOSCAPE ----*------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

## Define Output-Directory
exp_dir <- './Figure4/metadata'

##########################
# ---- SAVE HEALTHY ---- #
##########################

# Write files to .csv
write.table(edges_trimmed_new,
            paste0(exp_dir,'/','edges_all_healthyDonor_threetissue_TCRnetwork_cytoscape_OralclonesEqual1_RestBigger1.csv'),
            quote = F,
            sep=",",
            col.names = T,
            row.names = F)


