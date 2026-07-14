import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import matplotlib
matplotlib.use('TkAgg')
from scipy.spatial.distance import cdist
warnings.filterwarnings("ignore", category=UserWarning, message="Support for Awkward Arrays is currently experimental")
warnings.filterwarnings("ignore", category=UserWarning, message="Trying to modify attribute `.obs`")
warnings.filterwarnings("ignore", category=UserWarning, message="Variable names are not unique")
#warnings.filterwarnings("ignore", category=ImplicitModificationWarning, message="Trying to modify attribute `._uns` of view, initializing view as actual.")

import random
random.seed(42)

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ --------------------------- (1) LOAD GUT-DATASET ----------------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

## Genes of Interest
all_genes = ['HLA-DRB1','IL7','IFITM3','P2RY14','C1S','MT2A','GZMA','CD38','HLA-DQA1','IRF8','HLA-DMA','CD274','IRF4','NOD1','PELI1','GPR18','RSAD2','TRIM25','HERC6','GBP4',
             'FGL2','STAT1','IFIT3', 'IFIH1','CMTR1','IFI44L','TRAFD1','IFI30','IRF1','STAT2','MX1','SLAMF7','LAP3','SLC25A28','OAS2','JAK2','TXNIP','EPSTI1','HLA-B','HELZ2',
             'OGFR','ICAM1','OAS3','CD74','NFKB1','ST3GAL5','IRF9','TAPBP','TAP1','IFIT1','IRF2','ISG20','DHX58','RBCK1','HLA-G','PSMB8','CCL5','CASP3','APOL6',
             'RNF213', 'ADAR','SP110','IFIT2','PARP12','DDX60','MX2','EIF2AK2','MYD88','RNF31','PTPN1','MTHFD2','PIM1','PTPN2','IL2RB','PML','UPP1','NLRC5','IL7','IFITM3','C1S',
             'PARP9','RSAD2','CCRL2','TRIM25','HERC6','GBP4','CSF1','IFIT3','IFIH1','CMTR1','IFI44L','TRAFD1','IFI30','IRF1','STAT2','MX1','LAP3','SLC25A28','TRIM5','UBA7',
             'TXNIP','EPSTI1','HELZ2','OGFR','CD74','IRF9','TAP1','IRF2','ISG20','DHX58','PSMB8','OAS1','TENT5A','ADAR','SP110','IFIT2','PARP12','DDX60','EIF2AK2','RNF31','IL21','EOMES',
              'GZMB','IFNG','CCL4','GNLY','IL17F','LAG3'
             ]

genes=["GZMA","GZMB","IL21","IFNG","CD38","IRF4","CD274","GBP4","GPR18","STAT1","CMTR1","TRIM25","IFI44L","STAT2","SLAMF7","JAK2","CCL5",
           "ICAM1","NLRC5","NFKB1","CD74","IRF7","PARP9","IRF1","IFI30","IFIH1","MX1","IRF9","TRIM5","ISG20","IRF2","EOMES","CCL4"]

## Load Seurat pre-processed Object
adata = sc.read_10x_h5("./Figure5/metadata/HCplusCD_GUT_h5/HCplusCD_GUT.h5")
adata # --> n_obs × n_vars = 80481 × 36601


## Add metadata
metadata = pd.read_csv("./Figure5/metadata/HCplusCD_Gut_metadata.csv")
reactivity = pd.read_csv("./Figure5/metadata/HCplusCD_Gut_calbiReactivity_metadata.csv")

# Merge on cellID
combined_metadata = metadata.merge(reactivity,on="cellID",how="left",validate="one_to_one")
# Use cellID as the index
combined_metadata = combined_metadata.set_index("cellID")
# Check how well the IDs match adata
print("Cells in adata:", adata.n_obs)
print("Cells in metadata:", combined_metadata.shape[0])
print("Matching cells:",adata.obs_names.isin(combined_metadata.index).sum())

# Align metadata to the exact order of adata.obs_names
combined_metadata = combined_metadata.reindex(adata.obs_names)

# Add columns to adata.obs
adata.obs = adata.obs.join(combined_metadata)


## Subset to Th17 cells
ad_gut_th17calbi = adata[adata.obs["new_CaRxt_2aa"].isin(["Th17_CalbiRXT"])]

## Re-normalize after subsetting
ad_gut_th17calbi.layers['counts'] = ad_gut_th17calbi.X.copy()
sc.pp.normalize_total(ad_gut_th17calbi, target_sum=10000)
sc.pp.log1p(ad_gut_th17calbi)
ad_gut_th17calbi.layers['scaled'] = ad_gut_th17calbi.X.copy()
sc.pp.scale(ad_gut_th17calbi, max_value=10, layer='scaled')



### Add additional metadata
location_map = {"HC1":"colon",
                "HC2":"colon",
                "HC3":"colon",
                "HC4":"colon",
                "CD1":"colon",
                "CD2":"colon",
                "CD3":"colon",
                "CD4":"colon",
                "CD5":"colon",
                "CD6":"ileum",
                "CD7":"ileum",
                "CD8":"ileum",
                "CD9":"ileum",
                "CD10":"ileum"}

ad_gut_th17calbi.obs["donor_id"] = [name.split('_')[0] for name in ad_gut_th17calbi.obs_names.to_list()]
ad_gut_th17calbi.obs['biopsy_location'] = ad_gut_th17calbi.obs['donor_id'].map(location_map)
ad_gut_th17calbi.obs["status"] = (ad_gut_th17calbi.obs["donor_id"].str.extract(r"^(HC|CD)", expand=False))

## Generating Heatmap
sc.pl.matrixplot(ad_gut_th17calbi,genes,groupby=['status'],layer='scaled',cmap='RdBu_r',vmax=0.25, vmin=-0.25)


# Subsetting to whole Total th17 cells
# Make sure names match format (e.g., prefixed)
# If adata_filter.obs_names are already formatted correctly, you can proceed:
adata_filter_totalTH17 = adata[adata.obs["new_CaRxt_2aa"].isin(["Th17", "Th17_CalbiRXT"])]
print(adata_filter_totalTH17)

## Add donor_ID
adata_filter_totalTH17.obs["donor_id"] = [name.split('_')[0] for name in adata_filter_totalTH17.obs_names.to_list()]

## Renormalizing
adata_filter_totalTH17.layers['counts'] = adata_filter_totalTH17.X.copy() # .X -> contains cellbender corrected counts
sc.pp.normalize_total(adata_filter_totalTH17, target_sum=10000)
sc.pp.log1p(adata_filter_totalTH17)
adata_filter_totalTH17.layers['scaled'] = adata_filter_totalTH17.X.copy()
sc.pp.scale(adata_filter_totalTH17, max_value=10, layer='scaled')


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ --------------- (3) EXPORT MATRIXPLOT-DATA FOR GGPLOT2 ----------------- ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

###############################################
## -- COLON: EXPORT MATRIXPLOT DATA SHORT -- ##
###############################################

# Create the MatrixPlot object
mp_dnr_location_short = sc.pl.matrixplot(ad_gut_th17calbi,genes,groupby=['donor_id','biopsy_location'],layer='scaled',vmax=0.25, vmin=-0.25,return_fig=True)
mp_dnr_location_short = sc.pl.matrixplot(ad_gut_th17calbi,genes,groupby=['donor_id'],layer='scaled',vmax=0.25, vmin=-0.25,return_fig=True)
# Access the DataFrame with the plotted values
df_dnr_location_status_short = mp_dnr_location_short.values_df
# Save to CSV so you can import it into R
df_dnr_location_status_short.to_csv("./Figure5/metadata/data_MatrixPlot_GutTh17CalbiRXTcells_DonorLocation_short.csv")

##############################################
## -- COLON: EXPORT MATRIXPLOT DATA LONG -- ##
##############################################
## Create the MatrixPlot object
mp_donor_long = sc.pl.matrixplot(ad_gut_th17calbi,list(dict.fromkeys(all_genes)),groupby=['donor_id'],layer='scaled',vmax=0.25, vmin=-0.25,return_fig=True)
## Access the DataFrame with the plotted values
df_donor_long = mp_donor_long.values_df
## Save to CSV so you can import it into R
df_donor_long.to_csv("./Figure5/data_MatrixPlot_GutTh17CalbiRXTcells_Dnr_long.csv")


## Create the MatrixPlot object
mp_donor_long_totTH17 = sc.pl.matrixplot(adata_filter_totalTH17,list(dict.fromkeys(all_genes)),groupby=['donor_id'],layer='scaled',vmax=0.25, vmin=-0.25,return_fig=True)
## Access the DataFrame with the plotted values
df_donor_long_totTH17 = mp_donor_long_totTH17.values_df
## Save to CSV so you can import it into R
df_donor_long_totTH17.to_csv("./Figure5/data_MatrixPlot_Gut_TotalTH17_Dnr_long.csv")

