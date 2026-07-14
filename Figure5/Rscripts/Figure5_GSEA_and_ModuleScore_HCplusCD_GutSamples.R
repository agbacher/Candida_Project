#######################################################################################
# Author: Philipp Hofmann
# 
# This script will be used to reproduce Figure 5 - Panel D (Volcano Plot - GSEA)
# and Panel E (Pathogenicity Score Plot)
#
#######################################################################################


library(dplyr)
library(dandelionR)
library(Seurat)
library(miloR)
library(scRepertoire)
library(SingleCellExperiment)
library(scater)
library(ggplot2)
library(Lamian)

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
# ▓ ------------------------------- (2) LOAD VDJ DATA (WITH 2AA MISMATCH) ------------------------------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 


## LOAD VDJ_MASTER TABLE (TRA/TRB-aa sequences with  =<2aa mismatch in both, processed in scirpy using 'tcrdist' with cutoff=10)
airr_data_new <-read.csv('./Figure4/metadata/Revision_BOG_filtered_VDJoutput_identity_TCRdist_cutoff10_sameV.csv') %>% # airr_data
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
# ▓ ------------------------------- (3) ADD CALBI-RXT TO SEURAT METADATA ------------------------------------- ▓
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

#table(CDplusHC_gut[[c('donor','status_calbi_2aa_old')]])
table(CDplusHC_gut[[c('donor','status_calbi_2aa_new')]])
sum(CDplusHC_gut[['status_calbi_2aa_new']]=='CD_TRUE')

#########################################################
# -------- EXPORT CANDIDA_REACTIVITY METADATA --------- #
#########################################################
calbi_reactivity_meta <- CDplusHC_gut@meta.data %>%
                         tibble::rownames_to_column('cellID') %>%
                         dplyr::select(cellID,Calbi_RXT_TRATRB_2aaMis_NEW,new_CaRxt_2aa)

table(calbi_reactivity_meta['new_CaRxt_2aa'])
sum(calbi_reactivity_meta$new_CaRxt_2aa=='Th17_CalbiRXT')

## Write METADATA
write.csv(calbi_reactivity_meta, './Figure5/metadata/HCplusCD_Gut_calbiReactivity_metadata.csv',row.names = F)


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 
# ▓ ------------------------------- DEG-TESTING GUT: Th17-cells (CD vs. HC) ---------------------------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 

## Setting cell identities for DEG-testing
CDplusHC_gut$status_celltype_CaRxt_2aa <- paste0(CDplusHC_gut$status,'_',CDplusHC_gut$new_CaRxt_2aa)
Idents(CDplusHC_gut) <- 'status_celltype_CaRxt_2aa' 
table(Idents(CDplusHC_gut))

## Select genes for FindMarkers() [remove VDJ,Mito,Ribo and HB genes]
select_genes <- rownames(CDplusHC_gut)
genes2remove <- genes2remove <- c(grep("^TR[ABGD][VJ]", select_genes, value = TRUE),
                                  grep("^MT-", select_genes, value = TRUE),
                                  grep("^RP[SL][0-9]+", select_genes, value = TRUE),
                                  grep("^HB[^P]", select_genes, value = TRUE),
                                  grep("^HTOC[0-9]+$", select_genes, value = TRUE)) # remove hashtag genes "^HTOC[0-9]+$"
genes2use <- setdiff(select_genes,genes2remove)

## Run differential test using 
set.seed(12)
CDplusHC_gut <- NormalizeData(CDplusHC_gut)
cd_vs_hc_th17<- FindMarkers(CDplusHC_gut,
                            ident.1=c('CD_Th17','CD_Th17_CalbiRXT'),
                            ident.2=c('HC_Th17','HC_Th17_CalbiRXT'),
                            features=genes2use, 
                            logfc.threshold = 0.1) %>% 
                 scCustomize::Add_Pct_Diff() %>% 
                 mutate(gene=rownames(.),direction=ifelse(avg_log2FC>0,'UP_CD','DWN_CD')) 


################### GUT: Th17-cells (CD vs. HC) ##################################
###### MSigDB Enrichment plot (from: https://github.com/clareaulab/perffseq_reproducibility)
library(msigdbr)
library(fgsea)
m_df<- msigdbr(species = "Homo sapiens", category = "H")
unique(m_df$gs_name)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

### LOAD PATHOGENIC GENES (Th17-pathogenicity genes curated by PB)
patho_th17_df <- readxl::read_xlsx('./Figure5/metadata/Th17_pathogenicity_genes.xlsx')


## ADD CUSTOM pTH17 GENESET TO MSigDBList
fgsea_sets$PATHOGENIC_TH17 <- c(fgsea_sets$PATHOGENIC_TH17,unique(patho_th17_df$Gene))

## Setup DGEs for GSE-Analsis (adj_pval < 0.05)
cd_vs_hc_th17.genes<- cd_vs_hc_th17 %>%
                      filter(p_val_adj<0.05) %>%    # filter for significant genes
                      arrange(desc(avg_log2FC)) %>% # Descending order, starting with upregulated in CD
                      dplyr::select(gene, avg_log2FC)

## Generate vector with log2FC and Genes from DGE
vec_th17gut_CDvsHC <- cd_vs_hc_th17.genes$avg_log2FC
names(vec_th17gut_CDvsHC) <- cd_vs_hc_th17.genes$gene

## Perform GSEA using DGEs
set.seed(63)
fgseaRes_th17gut_CDvsHC <- fgseaMultilevel(fgsea_sets, stats = vec_th17gut_CDvsHC)

## Set p-threshold for labels
padj_thresh <- 0.1

## Generate Plot
fgseaRes_th17gut_CDvsHC <- fgseaRes_th17gut_CDvsHC %>% 
                          mutate(estimate=ifelse(padj<padj_thresh,'significant','NA')) %>% 
                          mutate(log10_padj = -log10(padj))  %>%
                          mutate(pathway=gsub('HALLMARK_','',pathway)) %>%
                          mutate(pathway=gsub('_',' ',pathway))  %>%
                          mutate(pathway=case_when(pathway=='INTERFERON GAMMA RESPONSE'~'IFNg response',
                                                   pathway=='INTERFERON ALPHA RESPONSE'~'IFNa response',
                                                   pathway=='INFLAMMATORY RESPONSE'~'Inflammatory response',
                                                   pathway=='PATHOGENIC TH17'~'pTh17',
                                                   pathway=='HOMEOSTATIC TH17'~'npTh17',
                                                   pathway=='CHOLESTEROL HOMEOSTASIS'~'Cholesterol homeostasis',
                                                   TRUE~pathway))

## Select top enriched pathway
pathway_highlight <-  fgseaRes_th17gut_CDvsHC %>%
                      slice_max(order_by = log10_padj, n = 4, with_ties = FALSE) %>%
                      pull(pathway)
pathway_highlight


## Plot parameters
dot_size=2.5
repel_text=2

## Generate Plot 
fgseaRes_th17gut_CDvsHC %>%
  arrange((padj)) %>%
  ggplot(aes(x = ES, y = -log10(padj))) + 
  geom_hline(yintercept = -log10(padj_thresh),linetype='dashed', color = "grey") +
  annotate("text", x = min(fgseaRes_th17gut_CDvsHC$ES)-0.4, y = -log10(padj_thresh) - 0.1, 
           label = paste0("pad-threshold (<",padj_thresh,")"), color = "grey", size = 2, hjust = 0) +
  geom_point(color = "black",size=dot_size) +  # Default black points
  geom_point(data = subset(fgseaRes_th17gut_CDvsHC, pathway %in% pathway_highlight), aes(x = ES, y = log10_padj),color = "red",size=dot_size) +  # Red = Top enriched pathway
  ggrepel::geom_text_repel(
    data = fgseaRes_th17gut_CDvsHC %>% filter(pathway %in% pathway_highlight),
    aes(label = pathway),
    size = repel_text, 
    box.padding = 0.5,  # Space around labels 0.3
    point.padding = 0.2, # Space between label and point 0.2
    max.overlaps = 10   # Adjust this to control the number of displayed labels
  ) + 
  ylim(0,1.5)+
  xlim(-0.8,0.8)+
  labs(x = "Enrichment score [HC (neg)/CD (pos)]", y = "-log10 (Padj.)",title = 'GSEA: DGEs from Colon-Th17 cells (HC vs. CD)') +
  BuenColors::pretty_plot(fontsize = 6) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))-> p_gsea
p_gsea

### Checking Genes for Gene_set enrichment terms
fgseaRes_th17gut_CDvsHC$leadingEdge[fgseaRes_th17gut_CDvsHC$pathway=='IFNg response']
fgseaRes_th17gut_CDvsHC$leadingEdge[fgseaRes_th17gut_CDvsHC$pathway=='IFNa response']
fgseaRes_th17gut_CDvsHC$leadingEdge[fgseaRes_th17gut_CDvsHC$pathway=='Inflammatory response']
fgseaRes_th17gut_CDvsHC$leadingEdge[fgseaRes_th17gut_CDvsHC$pathway=='pTh17']


##### SAVE PLOT
plt_dir <- './Figure5/Plots/'
cowplot::ggsave2(p_gsea, file = paste0(plt_dir,'gsea_VolcanoPlt_Th17cells_Gut_CDvsHC.pdf'), width = 3, height = 2.5)


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 
# ▓ ---------------------  DEG-TESTING GUT: Calbicans-reactive Th17-cells (CD vs. HC)  ----------------------- ▓ 
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 


## Select genes for FindMarkers() [remove VDJ,Mito,Ribo and HB genes]
select_genes <- rownames(CDplusHC_gut)
genes2remove <- genes2remove <- c(grep("^TR[ABGD][VJ]", select_genes, value = TRUE),
                                  grep("^MT-", select_genes, value = TRUE),
                                  grep("^RP[SL][0-9]+", select_genes, value = TRUE),
                                  grep("^HB[^P]", select_genes, value = TRUE),
                                  grep("^HTOC[0-9]+$", select_genes, value = TRUE)) # remove hashtag genes "^HTOC[0-9]+$"
genes2use <- setdiff(select_genes,genes2remove)

## Run differential test using [CD_Th17_CalbiRXT=258,HC_Th17_CalbiRXT=326]
set.seed(56)
cd_vs_hc_th17calbi2aa <- FindMarkers(CDplusHC_gut,ident.1=c('CD_Th17_CalbiRXT'),ident.2=c('HC_Th17_CalbiRXT'),
                                     features=genes2use, 
                                     logfc.threshold = 0.1) %>% 
                          scCustomize::Add_Pct_Diff() %>% 
                          mutate(gene=rownames(.),direction=ifelse(avg_log2FC>0,'UP_CD','DWN_CD')) 


################### GUT: Calbicans-reactive Th17-cells (CD vs. HC) ##################################
###### MSigDB Enrichment plot (from: https://github.com/clareaulab/perffseq_reproducibility)
library(msigdbr)
library(fgsea)
m_df<- msigdbr(species = "Homo sapiens", category = "H")
unique(m_df$gs_name)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

### LOAD PATHOGENIC GENES (Th17-pathogenicity genes curated by PB)
patho_th17_df <- readxl::read_xlsx('./Figure5/metadata/Th17_pathogenicity_genes.xlsx')


## ADD CUSTOM pTH17 GENESET TO MSigDBList
fgsea_sets$PATHOGENIC_TH17 <- c(fgsea_sets$PATHOGENIC_TH17,unique(patho_th17_df$Gene))

## Setup DGEs for GSE-Analsis (NO adj_pval filtering for Calbicans-reactive cells)
cd_vs_hc_th17calbi2aa.genes<- cd_vs_hc_th17calbi2aa %>%
                              arrange(desc(avg_log2FC)) %>% # Descending order, starting with upregulated in CD
                              dplyr::select(gene, avg_log2FC)

## Generate vector with log2FC and Genes from DGE
vec_th17calbi2aagut_CDvsHC <- cd_vs_hc_th17calbi2aa.genes$avg_log2FC
names(vec_th17calbi2aagut_CDvsHC) <- cd_vs_hc_th17calbi2aa.genes$gene

## Perform GSEA using DGEs
set.seed(37)
fgseaRes_th17calbi2aagut_CDvsHC <- fgseaMultilevel(fgsea_sets, stats = vec_th17calbi2aagut_CDvsHC)

## Set p-threshold for labels
padj_thresh <- 0.1

## Generate Plot
fgseaRes_th17calbi2aagut_CDvsHC <- fgseaRes_th17calbi2aagut_CDvsHC %>% 
                                    mutate(estimate=ifelse(padj<padj_thresh,'significant','NA')) %>% 
                                    mutate(log10_padj = -log10(padj)) %>%
                                    mutate(pathway=gsub('HALLMARK_','',pathway)) %>%
                                    mutate(pathway=gsub('_',' ',pathway))  %>%
                                    mutate(pathway=case_when(pathway=='INTERFERON GAMMA RESPONSE'~'IFNg response',
                                                             pathway=='INTERFERON ALPHA RESPONSE'~'IFNa response',
                                                             pathway=='INFLAMMATORY RESPONSE'~'Inflammatory response',
                                                             pathway=='PATHOGENIC TH17'~'pTh17',
                                                             pathway=='CHOLESTEROL HOMEOSTASIS'~'Cholesterol homeostasis',
                                                             TRUE~pathway))
## Select top enriched pathway
pathway_highlight <-  fgseaRes_th17calbi2aagut_CDvsHC %>%
                      slice_max(order_by = log10_padj, n = 3, with_ties = FALSE) %>%
                      pull(pathway)
pathway_highlight


## Plot parameters

dot_size=2.5
repel_text=2
## Generate Plot 
fgseaRes_th17calbi2aagut_CDvsHC %>% 
  arrange((padj)) %>%
  #mutate(log_padj = -log10(padj)) %>%
  ggplot(aes(x = ES, y = -log10(padj))) + 
  geom_hline(yintercept = -log10(padj_thresh),linetype='dashed', color = "grey") +
  annotate("text", x = min(fgseaRes_th17calbi2aagut_CDvsHC$ES)-0.2, y = -log10(padj_thresh) - 0.1, 
           label = paste0("pad-threshold (<",padj_thresh,")"), color = "grey", size = 2, hjust = 0) +
  geom_point(color = "black",size=dot_size) +  # Default black points
  geom_point(data = subset(fgseaRes_th17calbi2aagut_CDvsHC, pathway %in% pathway_highlight), aes(x = ES, y = log10_padj),color = "red",size=dot_size) +  # Red = Top enriched pathway
  ggrepel::geom_text_repel(
    data = subset(fgseaRes_th17calbi2aagut_CDvsHC, pathway %in% pathway_highlight),
    aes(label = pathway),
    size = repel_text, 
    box.padding = 0.3,  # Space around labels
    point.padding = 0.2, # Space between label and point
    max.overlaps = 10   # Adjust this to control the number of displayed labels
  ) + 
  ylim(0,1.5)+
  xlim(-0.8,0.8)+
  labs(x = "Enrichment score [HC (neg)/CD (pos)]", y = "-log10 (Padj.)",title = 'GSEA: DEGs from Colon-CalbiRXT-Th17 cells (HCvs.CD)') +
  BuenColors::pretty_plot(fontsize = 6) +
  theme(plot.title = element_text(size=6),
        plot.background = element_rect(fill = "transparent", colour = NA)) -> p_gsea_calbi
p_gsea_calbi


### Checking Genes for Gene_set enrichment terms
fgseaRes_th17calbi2aagut_CDvsHC$leadingEdge[fgseaRes_th17calbi2aagut_CDvsHC$pathway=='IFNg response']
fgseaRes_th17calbi2aagut_CDvsHC$leadingEdge[fgseaRes_th17calbi2aagut_CDvsHC$pathway=='IFNa response']
fgseaRes_th17calbi2aagut_CDvsHC$leadingEdge[fgseaRes_th17calbi2aagut_CDvsHC$pathway=='pTh17']

intersect(fgseaRes_th17calbi2aagut_CDvsHC$leadingEdge[fgseaRes_th17calbi2aagut_CDvsHC$pathway=='IFNg response'][[1]],fgseaRes_th17calbi2aagut_CDvsHC$leadingEdge[fgseaRes_th17calbi2aagut_CDvsHC$pathway=='IFNa response'][[1]])

##### SAVE PLOT
plt_dir <- './Figure5/Plots/'
cowplot::ggsave2(p_gsea_calbi, file = paste0(plt_dir,'gsea_VolcanoPlt_Th17CalbiRXTcells_Gut_CDvsHC.pdf'), width = 3, height = 2.5)


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 
# ▓ ---------------------------- 2) MODULE SCORE WITH GSEA-ENRICHED GENES ------------------------------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 

## subset to Th17 and Th17-C.albi-reactive cells
Idents(CDplusHC_gut) <- 'new_CaRxt_2aa'
sub <- subset(CDplusHC_gut,ident = c('Th17','Th17_CalbiRXT'))
table(sub$new_CaRxt_2aa)

## Get Geneset for Module-Scoring
df_Th17Calbi <- data.frame(gene=c(unname(fgseaRes_th17calbi2aagut_CDvsHC$leadingEdge[fgseaRes_th17calbi2aagut_CDvsHC$pathway=='IFNg response'])[[1]],
                                  unname(fgseaRes_th17calbi2aagut_CDvsHC$leadingEdge[fgseaRes_th17calbi2aagut_CDvsHC$pathway=='IFNa response'])[[1]],
                                  unname(fgseaRes_th17calbi2aagut_CDvsHC$leadingEdge[fgseaRes_th17calbi2aagut_CDvsHC$pathway=='pTh17'])[[1]]
                                  ),
                           geneset=c(
                             rep("IFNg response",length(unname(fgseaRes_th17calbi2aagut_CDvsHC$leadingEdge[fgseaRes_th17calbi2aagut_CDvsHC$pathway=='IFNg response'])[[1]])),
                             rep("IFNa response",length(unname(fgseaRes_th17calbi2aagut_CDvsHC$leadingEdge[fgseaRes_th17calbi2aagut_CDvsHC$pathway=='IFNa response'])[[1]])),
                             rep("pTh17",length(unname(fgseaRes_th17calbi2aagut_CDvsHC$leadingEdge[fgseaRes_th17calbi2aagut_CDvsHC$pathway=='pTh17'])[[1]]))
                             )
                           )
unique(df_Th17Calbi$gene)
length(unique(df_Th17Calbi$gene))

## Add Module
sub <- NormalizeData(sub)
sub <- AddModuleScore(sub,features = list(unique(df_Th17Calbi$gene)),name='patho_score_')


######################################
# --- GENERATE MODULE-SCORE PLOT --- #
######################################
# NOTE: The plot in the publication was generate using GraphPad Prism. However,
#       the following code shows an alternative version using ggplot2.


## Generate Seurat ModuleScore-ViolinPlot to get plotting data 
VlnPlot(sub,feature='patho_score_1',
        split.by = 'status',
        split.plot = T,
        pt.size = 0,
        idents = c('Th17','Th17_CalbiRXT')) -> p

# Get dataframe from ggplot-object generated using Seurat and make custom plot
data <- p[1]$data
data <- data %>% mutate(split = factor(split, levels = c("HC", "CD")))

## Generate Plot
library(ggplot2)
library(ggpubr)
ggplot(data, aes(x = split, y = patho_score_1, fill = split)) +
  geom_violin(width = 0.75, trim = FALSE,color = "black",linewidth = 0.4) +
  geom_boxplot( width = 0.25, outlier.shape = NA,color = "black",linewidth = 0.4,alpha = 0.5) +
  stat_compare_means( comparisons = list(c("HC", "CD")),
                      method = "wilcox.test",
                      label = "p.format",
                      size = 2
                      ) +
  scale_fill_manual(values = c('HC'="#9b9d9d",'CD'="#797BA8")) +
  scale_color_manual(values = c('HC'="black",'CD'="#434579")) +
  labs(x = NULL,y = "Pathogenicity Score") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),
        plot.background = element_rect(fill = "transparent", colour = NA)
        ) + 
  facet_grid(. ~ ident, scales = 'free') -> p_ModuleScore
p_ModuleScore

##### SAVE PLOT
plt_dir <- './Figure5/Plots/'
cowplot::ggsave2(p_ModuleScore, file = paste0(plt_dir,'moduleScore_GSAE_Gut_CDvsHC.pdf'), width = 3.5, height = 3.8)













