# 31.10.2023 - Charlotte Boys

# %% Import libraries
import scanpy as sc
import platform as platform
import os as os
import yaml as yaml
import numpy as np
import pandas as pd

# load muon and mofax
import muon as mu
import mofax as mofa

from helper_functions import adata_to_views, preprocess_data

# %% Parse snakemake object TODO: update for new repo
if 'snakemake' in locals():
    h5ad_input = snakemake.input['h5ad_input']
    coldata_path = snakemake.input['coldata_path']
    h5ad_output = snakemake.output['h5ad_output']
    h5mu_output = snakemake.output['h5mu_output']
    sample_key = snakemake.params['sample_key']
    condition_key = snakemake.params['condition_key']
    groupby = snakemake.params['groupby']
    drop_views = snakemake.params['drop_views']
    min_cells = snakemake.params['min_cells']
    batch_key = snakemake.params['batch_key']
    if snakemake.wildcards['subset'] == 'all':
        sample_subset = 'all'
        cell_subset = 'all'
    else:
        sample_subset = snakemake.params['sample_subset']
        cell_subset = snakemake.params['cell_subset']
else:
    sample_key = 'EdgarID'
    condition_key = 'Cohort'
    groupby = 'celltype'
    sample_subset = 'all'
    cell_subset = 'all'
    min_cells = 50
    batch_key = 'Project1'
    drop_views = ['NA'] # ['DTL-2', 'PT-3']

    h5ad_input = 'results/preprocessing/JulioNereid.h5ad'
    coldata_path = 'data/Metadata_For_Charlotte.csv'

    if platform.system() == 'Linux':
        h5ad_input = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', h5ad_input)
        coldata_path = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', coldata_path)

    elif platform.system() == 'Darwin':
        h5ad_input = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', h5ad_input)
        coldata_path = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', coldata_path)

# %% Read in AnnData object
adata = sc.read_h5ad(h5ad_input)

# %%
# Keep only samples in subset
if sample_subset != 'all':
    adata = adata[adata.obs[sample_key].isin(sample_subset)]

# %% Keep only cell types we're interested in
if cell_subset != 'all':
    adata = adata[adata.obs[groupby].isin(cell_subset)]

# %% Add coldata to adata.obs
coldata = pd.read_csv(coldata_path)
#coldata = coldata.rename(columns= {"sample_ID": "Sample"})
coldata_rna = coldata[coldata[sample_key].isin(adata.obs[sample_key])]

# %%
intersecting_columns = list(set(adata.obs) & set(coldata_rna))
print("INFO: Intersecting columns between metadata and adata.obs")
print(intersecting_columns)
# %%
if len(intersecting_columns) > 1:
    new_obs = adata.obs.reset_index().merge(coldata_rna, how = "left", on = intersecting_columns)
else:
    new_obs = adata.obs.reset_index().merge(coldata_rna, how = "left", on = [sample_key])
# %%
new_obs = new_obs.drop_duplicates(subset='index').set_index('index')
# for all columns of type object, convert to string
object_columns = new_obs.select_dtypes(include=['object']).columns
new_obs[object_columns] = new_obs[object_columns].astype(str)
# for all columns of type float64, convert to float32 
float_columns = new_obs.select_dtypes(include=['float64']).columns
new_obs[float_columns] = new_obs[float_columns].astype('float32')

# %%
adata.obs = new_obs.reindex(adata.obs.index)

print("INFO: adata.obs columns")
print(list(adata.obs))
# %% 
# Require minimum genes per cell and minimum cells per gene
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
# %%
# Based on low number of samples covered (15-17)
# adata = adata[~adata.obs[groupby].isin(['Plasma_cells', 'Proliferating'])]

# Crucial for being able to save and read the mudata object later on
# Otherwise you'll get an error message like TypeError: AnnData.__init__() got an unexpected keyword argument 'TAL'
# If one of your cell type views is called EC/TAL
# Or TypeError: AnnData.__init__() got an unexpected keyword argument 'P'
# if one of your cell type  views is called VSM/P
adata.obs[groupby] = adata.obs[groupby].str.replace('/', '_')
# Remove based on low number of samples covered
adata = adata[~adata.obs[groupby].isin(drop_views)]
#adata.X = np.round(adata.X)
adata.layers['counts'] = adata.X
# %% Create a multi-view structure
mdata_original = adata_to_views(adata,
                        groupby=groupby,
                        sample_key=sample_key,
                        obs_keys=['ID', 'Project1', 'Cohort', 'Prep', 'Method', 'Age', 'Race',
                        'eGFRatBx_NEPTUNE', 'UPCRatBx_NEPTUNE', 'Sex',
                        'RAASBlockPreBx', 'InterstitialFibrosis', 'TubularAtrophy'], # add those to mdata.obs
                        min_prop=0.05, # min nnz values (filter features)
                        min_smpls=3, # min samples per view (filter features)
                        min_cells=min_cells, # min cells per view (filter samples)
                        min_counts=100, # min counts per view (filter samples)
                        mode='sum', # mode of aggregation
                        verbose=True,
                        large_n=5, # edgeR-like filtering
                        min_total_count=15,
                        min_count=10,
                        keep_psbulk_stats = True,
                        skip_checks = True # Ignore that counts are not integers (because of SoupX correction)
                        )

# Enforce minimum number of views per sample:
# lambda function to check if sample is in a view
in_view = lambda view: mdata_original.obs_names.isin([*mdata_original.mod[view].obs_names])*1
# Check for two or more views per sample
mask = sum(list(map(in_view, [*mdata_original.mod.keys()]))) >= 2
pass_samples = mdata_original.obs_names[mask]
for current_mod in mdata_original.mod.keys():
    s = mdata_original.mod[current_mod].obs_names.isin(pass_samples)
    mdata_original.mod[current_mod] = mdata_original.mod[current_mod][s, :].copy()
mdata_original.update()
mdata_original.uns['psbulk_stats'] = mdata_original.uns['psbulk_stats'][mask]

print("INFO: debugging mdata")
print(list(mdata_original.obs))
print("INFO: debugging views")
print(list(mdata_original.mod['PT'].obs))

# %% Pre-process the pseudobulk profiles
mdata = preprocess_data(mdata_original,
                        target_sum = 10000,
                        remove_mitochondrial = False,
                        exclude_highly_expressed = False,
                        batch_key = batch_key, # first find HVGs within each batch then merge
                        n_top_genes = None) # default number of top genes

if 'snakemake' in locals():
    print('INFO: Writing AnnData to {0}'.format(h5ad_output))
    adata.write(h5ad_output)
    print('INFO: Writing MuData to {0}'.format(h5mu_output))
    mdata.write(h5mu_output)
# %%
