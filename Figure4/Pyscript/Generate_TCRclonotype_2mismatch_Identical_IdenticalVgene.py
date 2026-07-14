import anndata as ad
import scanpy as sc
import muon as mu
import scirpy as ir
import numpy as np
import pandas as pd
import re
import matplotlib
matplotlib.use('TkAgg')
import warnings
warnings.filterwarnings("ignore", category=UserWarning, message="Support for Awkward Arrays is currently experimental")
warnings.filterwarnings("ignore", category=UserWarning, message="Trying to modify attribute `.obs`")
warnings.filterwarnings("ignore", category=UserWarning, message="Variable names are not unique")
#warnings.filterwarnings("ignore", category=ImplicitModificationWarning, message="Trying to modify attribute `._uns` of view, initializing view as actual.")

import random
random.seed(42)

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ --------------------------- (1) LOAD BOG-PROJECT DATA ------------------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

### Step1: Import all BOG-Project Data and generate multimodal MuData object (GEX,VDJ)
samples = ['HC1','HC2','HC3','HC4',
           'CD1','CD2','CD3','CD4','CD5','CD6','CD7','CD8','CD9','CD10'
           ]


# Create a list of AnnData objects (one for each sample)
samples_list = []
for sample in samples:
    #print(sample)

    if sample in ['CD1','CD3','CD4','CD5']:
        datasets = ['Th1_blood','Th17_blood','Calbi_blood','gut','oral','oral_exp2']
    else:
        datasets = ['Th1_blood','Th17_blood','Calbi_blood','gut','oral']

    for location in datasets:
        sample_name = f'{sample}_{location}_cellranger_out'
        samples_list.append(sample_name)



# Create a list of AnnData objects (one for each sample)
adatas_tcr_bog = {}
adatas_gex_bog = {}
for sample in samples_list:
    print(f'Sample: {sample}')

    # Skip Oral samples for Donors that were not collected
    if any(sub in sample for sub in ['CD10_oral','CD7_oral','CD6_oral']):
        continue
        print(f'Sample {sample} skipped - lacking data')

    adata_tcr_bog = ir.io.read_10x_vdj(f'./cellranger_out/{sample}/vdj_t/filtered_contig_annotations.csv')

    ## Load Anndata object with cellbender-corrected UMI-counts in ad.X and cellranger raw UMI-counts in ad.raw
    cr_gene_filtered_h5 = f'./cellranger_out/{sample}/count/sample_filtered_feature_bc_matrix.h5'
    print("...Loading CellRanger MTX-file done!")
    adata_gex_bog = sc.read_10x_h5(cr_gene_filtered_h5, gex_only=True)
    adata_gex_bog.var_names_make_unique()

    adatas_tcr_bog[sample] = adata_tcr_bog
    adatas_gex_bog[sample] = adata_gex_bog

# Merge anndata objects and make multimodal MuData Object
import anndata
ad_gex_bog = anndata.concat(adatas_gex_bog, index_unique="_")
ad_tcr_bog = anndata.concat(adatas_tcr_bog, index_unique="_")

# Generating Multiome (GEX+TCR) Muon-Data Object
mdata = mu.MuData({"gex": ad_gex_bog,
                   "airr": ad_tcr_bog
                   })

# Adding metadata
# Define the mapping between samples and desired names
sample_name_mapping = {
    '1_HC1_Th1_blood_cellranger_out': 'HC1_Th1',
    '1_HC1_Th17_blood_cellranger_out': 'HC1_Th17',
    '1_HC1_Calbi_blood_cellranger_out': 'HC1_Calbi',
    '1_HC2_Th1_blood_cellranger_out': 'HC2_Th1',
    '1_HC2_Th17_blood_cellranger_out': 'HC2_Th17',
    '1_HC2_Calbi_blood_cellranger_out': 'HC2_Calbi',
    '1_HC3_Th1_blood_cellranger_out': 'HC3_Th1',
    '1_HC3_Th17_blood_cellranger_out': 'HC3_Th17',
    '1_HC3_Calbi_blood_cellranger_out': 'HC3_Calbi',
    '1_HC4_Th1_blood_cellranger_out': 'HC4_Th1',
    '1_HC4_Th17_blood_cellranger_out': 'HC4_Th17',
    '1_HC4_Calbi_blood_cellranger_out': 'HC4_Calbi',
    '1_CD1_Th1_blood_cellranger_out': 'CD1_Th1',
    '1_CD1_Th17_blood_cellranger_out': 'CD1_Th17',
    '1_CD1_Calbi_blood_cellranger_out': 'CD1_Calbi',
    '1_CD2_Th1_blood_cellranger_out': 'CD2_Th1',
    '1_CD2_Th17_blood_cellranger_out': 'CD2_Th17',
    '1_CD2_Calbi_blood_cellranger_out': 'CD2_Calbi',
    '1_CD5_Th1_blood_cellranger_out': 'CD5_Th1',
    '1_CD5_Th17_blood_cellranger_out': 'CD5_Th17',
    '1_CD5_Calbi_blood_cellranger_out': 'CD5_Calbi',
    '1_CD3_Th1_blood_cellranger_out': 'CD3_Th1',
    '1_CD3_Th17_blood_cellranger_out': 'CD3_Th17',
    '1_CD3_Calbi_blood_cellranger_out': 'CD3_Calbi',
    '1_CD4_Th1_blood_cellranger_out': 'CD4_Th1',
    '1_CD4_Th17_blood_cellranger_out': 'CD4_Th17',
    '1_CD4_Calbi_blood_cellranger_out': 'CD4_Calbi',
    '1_CD6_Th1_blood_cellranger_out': 'CD6_Th1',
    '1_CD6_Th17_blood_cellranger_out': 'CD6_Th17',
    '1_CD6_Calbi_blood_cellranger_out': 'CD6_Calbi',
    '1_CD7_Th1_blood_cellranger_out': 'CD7_Th1',
    '1_CD7_Th17_blood_cellranger_out': 'CD7_Th17',
    '1_CD7_Calbi_blood_cellranger_out': 'CD7_Calbi',
    '1_CD8_Th1_blood_cellranger_out': 'CD8_Th1',
    '1_CD8_Th17_blood_cellranger_out': 'CD8_Th17',
    '1_CD8_Calbi_blood_cellranger_out': 'CD8_Calbi',
    '1_CD9_Th1_blood_cellranger_out': 'CD9_Th1',
    '1_CD9_Th17_blood_cellranger_out': 'CD9_Th17',
    '1_CD9_Calbi_blood_cellranger_out': 'CD9_Calbi',
    '1_CD10_Th1_blood_cellranger_out': 'CD10_Th1',
    '1_CD10_Th17_blood_cellranger_out': 'CD10_Th17',
    '1_CD10_Calbi_blood_cellranger_out': 'CD10_Calbi',
    '1_HC1_gut_cellranger_out': 'HC1_gut',
    '1_HC4_gut_cellranger_out': 'HC4_gut',
    '1_HC2_gut_cellranger_out': 'HC2_gut',
    '1_HC3_gut_cellranger_out': 'HC3_gut',
    '1_CD1_gut_cellranger_out': 'CD1_gut',
    '1_CD2_gut_cellranger_out': 'CD2_gut',
    '1_CD5_gut_cellranger_out': 'CD5_gut',
    '1_CD3_gut_cellranger_out': 'CD3_gut',
    '1_CD4_gut_cellranger_out': 'CD4_gut',
    '1_CD6_gut_cellranger_out': 'CD6_gut',
    '1_CD7_gut_cellranger_out': 'CD7_gut',
    '1_CD8_gut_cellranger_out': 'CD8_gut',
    '1_CD9_gut_cellranger_out': 'CD9_gut',
    '1_CD10_gut_cellranger_out': 'CD10_gut',
    '1_HC1_oral_cellranger_out': 'HC1_oral',
    '1_HC4_oral_cellranger_out': 'HC4_oral',
    '1_HC2_oral_cellranger_out': 'HC2_oral',
    '1_HC3_oral_cellranger_out': 'HC3_oral',
    '1_CD1_oral_cellranger_out': 'CD1_oral',
    '1_CD1_oral_exp2_cellranger_out': 'CD1_oral',
    '1_CD2_oral_cellranger_out': 'CD2_oral',
    '1_CD5_oral_cellranger_out': 'CD5_oral',
    '1_CD5_oral_exp2_cellranger_out': 'CD5_oral',
    '1_CD3_oral_cellranger_out': 'CD3_oral',
    '1_CD3_oral_exp2_cellranger_out': 'CD3_oral',
    '1_CD4_oral_cellranger_out': 'CD4_oral',
    '1_CD4_oral_exp2_cellranger_out': 'CD4_oral',
    '1_CD8_oral_cellranger_out': 'CD8_oral',
    '1_CD9_oral_cellranger_out': 'CD9_oral',
}


# Set global metadata on `mdata.obs`
mdata.obs["sample"] = [name.split('_')[1] for name in mdata.obs_names.to_list()]
# Apply the mapping to the 'sample' column in mdata.obs
mdata.obs['sample_map'] = mdata.obs['sample'].map(sample_name_mapping)
mdata.obs['donor'] = [name.split('_')[0] for name in mdata.obs.sample_map.to_list()]
mdata.obs['stimu_tissue'] = [name.split('_')[1] for name in mdata.obs.sample_map.to_list()]
mdata.obs.groupby("donor", dropna=False).size()
mdata.obs.groupby("stimu_tissue", dropna=False).size()

## Sanity check if all samples are included
sample_list_1 = pd.Series(list(sample_name_mapping.keys())).unique().tolist()
sample_list_2 = mdata.obs['sample'].unique().tolist()
set1 = set(sample_list_1)
set2 = set(sample_list_2)
# Items in mapping but not in obs
diff_1_not_2 = set1 - set2
# Items in obs but not in mapping
diff_2_not_1 = set2 - set1
## If print is empty then all samples are included
print("In mapping but not in obs:", diff_1_not_2)
print("In obs but not in mapping:", diff_2_not_1)

##### SAVING MuDATA OBJECT
#mdata.write("./Analysis/h5ad_files/ProjBOG_GEXandVDJ_MuData_BLOODandGUTandORAL.h5mu")
#sc.write('./Analysis/h5ad_files/ProjBOG_GEXandVDJ_MuData_BLOODandGUTandORAL.h5ad', mdata['gex'], compression="gzip")

##### LOAD MuDATA OBJECT
mdata =mu.read_h5mu("./Analysis/h5ad_files/ProjBOG_GEXandVDJ_MuData_BLOODandGUTandORAL.h5mu")

# Apply the mapping to the 'sample' column in mdata.obs
mdata.obs.groupby("donor", dropna=False).size()
mdata.obs.groupby("stimu_tissue", dropna=False).size()


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------- (2) PRE-PROCESSING OF AIRR-DATA  ------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

adata_ir = mdata.mod["airr"]
adata_ir.obs["stimu_tissue"] = mdata.obs["stimu_tissue"].astype("category")


## AIRR-data pre-processing
ir.pp.index_chains(adata_ir)
ir.tl.chain_qc(adata_ir)

## Plot Distribution for AIRR-status
ir.pl.group_abundance(adata_ir, groupby="chain_pairing", target_col="stimu_tissue")

## Filter cells based TCR chain-pairing (remove Multichain=likely doublets, Orphan chain=single alpha or beta, no IR)
mu.pp.filter_obs(adata_ir, "chain_pairing", lambda x: ~np.isin(x, ["no IR", "orphan VDJ", "orphan VJ", "multichain"]))

## Check How many cells surfived fitlering
adata_ir 

## Sanity check after filtering
ir.pl.group_abundance(adata_ir, groupby="chain_pairing", target_col="stimu_tissue")

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------- (3) CALCULATE IDENTICAL TRA/TRB  ------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

# using default parameters, `ir_dist` will compute nucleotide sequence identity
ir.pp.ir_dist(adata_ir,
              sequence="aa",
              metric='identity',
              key_added='ir_aa_identity')
ir.tl.define_clonotypes(adata_ir,
                        receptor_arms="all",
                        dual_ir="primary_only",
                        distance_key='ir_aa_identity',
                        key_added='cc_aa_identity')


# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------- (4) CALCULATE IDENTICAL TRA/TRB SAME V-GENE ------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

# using default parameters, `ir_dist` will compute nucleotide sequence identity
#ir.pp.ir_dist(mdata,
#              sequence="aa",
#              metric='identity',
#              key_added='ir_aa_identity')

ir.tl.define_clonotypes(adata_ir,
                        receptor_arms="all",
                        dual_ir="primary_only",
                        same_v_gene = True,
                        distance_key='ir_aa_identity',
                        key_added='cc_aa_sameV_identity')



# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------- (5) CALCULATE TCR-DISTANCE ON TRA/TRB WITH 2-AA MISTMATCH  ------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

### USING IMPLEMENTATION IN SCIRPY
# computing CD3 neighborhood based on clonotype clusters using TCRdist-algorithm and allowing N=2 amino acid changes (cutoff=10)
ir.pp.ir_dist(adata_ir,
              metric="tcrdist",
              sequence="aa",
              cutoff=10,
              key_added='ir_dist_2aa_tcrdist')

ir.tl.define_clonotype_clusters(adata_ir,
                                sequence="aa",
                                metric="tcrdist",
                                receptor_arms="all",
                                dual_ir="primary_only",
                                distance_key='ir_dist_2aa_tcrdist',
                                key_added='cc_2aa_tcrdist'
                                )

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------- (6) CALCULATE TCR-DISTANCE ON TRA/TRB WITH 2-AA MISTMATCH  (same Vgene) ------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

### Defining clonotypes based TCRdist with 2aa mismatch, respecting V-gene
ir.tl.define_clonotype_clusters(adata_ir,
                                sequence="aa",
                                metric="tcrdist",
                                receptor_arms="all",
                                dual_ir="primary_only",
                                same_v_gene = True,
                                distance_key='ir_dist_2aa_tcrdist',
                                key_added='cc_2aa_sameV_tcrdist'
                                )

# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
# ▓ ------------- (7) EXPORTING RESULTS FOR AIRR-ANALYIS FOR IMPORT IN Rsoftware ------------ ▓
# ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

## Get VDJ data
dict_rename_tcrdist = {
    "VJ_1_junction_aa": "cdr3_a_aa",
    "VJ_1_v_call": "v_a_gene",
    "VJ_1_j_call": "j_a_gene",
    "VDJ_1_junction_aa": "cdr3_b_aa",
    "VDJ_1_v_call": "v_b_gene",
    "VDJ_1_j_call": "j_b_gene"
}
airr_data = ir.get.airr(adata_ir, ['v_call','j_call','junction_aa'], ['VJ_1', 'VDJ_1'])
airr_data = airr_data.rename(columns=dict_rename_tcrdist)

# Check cell order
(mdata.obs_names == adata_ir.obs_names).all()

# Join into adata_ir.obs
adata_ir.obs = adata_ir.obs.join(airr_data)

# Adding the AIRR-data to general observation of the Object for ease of export
cols = [
    "v_a_gene","j_a_gene","cdr3_a_aa",
    "v_b_gene","j_b_gene","cdr3_b_aa",
    "cc_aa_identity","cc_aa_identity_size",
    "cc_aa_sameV_identity","cc_aa_sameV_identity_size",
    "cc_2aa_tcrdist","cc_2aa_tcrdist_size",
    "cc_2aa_sameV_tcrdist","cc_2aa_sameV_tcrdist_size"
]
## adding results columns to muData
for c in cols:
    mdata.obs[c] = adata_ir.obs[c]

# Export TCR-clone data
cols_to_export=['donor','sample_map','stimu_tissue',
                'v_a_gene', 'j_a_gene', 'cdr3_a_aa', 'v_b_gene', 'j_b_gene', 'cdr3_b_aa',
                'cc_aa_identity','cc_aa_identity_size',
                'cc_aa_sameV_identity','cc_aa_sameV_identity_size',
                'cc_2aa_tcrdist', 'cc_2aa_tcrdist_size',
                'cc_2aa_sameV_tcrdist', 'cc_2aa_sameV_tcrdist_size'
                ]
mdata.obs[cols_to_export].\
        to_csv(f'./Figure4/metadata/Revision_BOG_filtered_VDJoutput_identity_TCRdist_cutoff10_sameV.csv')
