# 31.10.2023 - Charlotte Boys

# %% Import libraries
import scanpy as sc
import platform as platform
import os as os
import yaml as yaml
import numpy as np
import pandas as pd
import muon as mu
import mofax as mofa
from helper_functions import adata_to_views, preprocess_data

# %% Parse snakemake object TODO: update for new repo
if 'snakemake' in locals():
    # Input
    h5ad_input = snakemake.input['h5ad_input']
    coldata_path = snakemake.input['coldata_path']
    annotation_dict = snakemake.input['annotation_dict']
    # Output
    h5ad_output = snakemake.output['h5ad_output']
    h5mu_output = snakemake.output['h5mu_output']
    # Params
    sample_key = snakemake.params['sample_key']
    groupby = snakemake.params['groupby']
    drop_views = snakemake.params['drop_views']
    min_cells = snakemake.params['min_cells']
    batch_key = snakemake.params['batch_key']
    condition_key = snakemake.params['condition_key']
    sample_subset = snakemake.params['sample_subset']
else:
    sample_key = 'OCEAN_MS_Exp_ID'
    groupby = 'LR_Cluster2'
    sample_subset = ['Exp18', 'Exp19', 'Exp154', 'Exp155', 'Exp2', 'Exp21', 'Exp175', 'Exp169',
      'Exp167', 'Exp84', 'Exp85', 'Exp152', 'Exp131', 'Exp109', 'Exp214', 'Exp216', 'Exp219',
      'Exp221', 'Exp226', 'Exp228', 'Exp150', 'Exp148', 'Exp146', 'Exp145', 'Exp144', 'Exp141',
      'Exp138', 'Exp135', 'Exp133', 'Exp16', 'Exp129', 'Exp128', 'Exp126', 'Exp121', 'Exp123',
      'Exp119', 'Exp117', 'Exp115', 'Exp113', 'Exp111', 'Exp1', 'Exp197', 'Exp199', 'Exp201',
        'Exp203', 'Exp206', 'Exp208', 'Exp210', 'Exp212', 'Exp185', 'Exp186', 'Exp184', 'Exp56',
        'Exp93', 'Exp46', 'Exp36', 'Exp26', 'Exp158', 'Exp139', 'Exp92', 'Exp54', 'Exp183', 'Exp45',
          'Exp25', 'Exp35', 'Exp91', 'Exp182', 'Exp15', 'Exp58', 'Exp42', 'Exp31', 'Exp30', 'Exp29',
          'Exp28', 'Exp204', 'Exp94', 'Exp12', 'Exp173', 'Exp74', 'Exp64', 'Exp60', 'Exp61', 'Exp59',
            'Exp57', 'Exp55', 'Exp8', 'Exp9', 'Exp53', 'Exp52', 'Exp51', 'Exp50', 'Exp49', 'Exp48', 'Exp47',
            'Exp44', 'Exp6', 'Exp7', 'Exp43', 'Exp41', 'Exp40', 'Exp39', 'Exp38', 'Exp37', 'Exp34', 'Exp4',
              'Exp5', 'Exp33', 'Exp32']
    min_cells = 50
    batch_key = 'Project1'
    condition_key = 'Cohort_Charlotte'
    drop_views = ['tPC-IC']
    h5ad_input = 'data/OCEAN_v3_Nu_102025a_CZI.h5ad'
    coldata_path = 'data/metadata.csv'
    annotation_dict = 'annotations/OCEAN_LR.csv'

    if platform.system() == 'Linux':
        h5ad_input = os.path.join('/mnt/sds-hd/sd22b002/projects/OCEAN_MCFA', h5ad_input)
        coldata_path = os.path.join('/mnt/sds-hd/sd22b002/projects/OCEAN_MCFA', coldata_path)
    elif platform.system() == 'Darwin':
        h5ad_input = os.path.join('/Users/charlotteboys/GitHub/OCEAN_MCFA', h5ad_input)
        coldata_path = os.path.join('/Users/charlotteboys/GitHub/OCEAN_MCFA', coldata_path)

# %% Read in AnnData object
adata = sc.read_h5ad(h5ad_input)

# %% Add annotations from annotation dict
mapping_df = pd.read_csv(annotation_dict, sep=';')
print("INFO: mapping")
celltype_mapping2 = dict(zip(mapping_df['LR_Cluster'],
                                       mapping_df['LR_Cluster2']))
print(celltype_mapping2)
# %%
adata.obs['LR_Cluster2'] = adata.obs['LR_Cluster'].map(celltype_mapping2)
print("INFO: adata.obs")
print(adata.obs.head())

# %%
# Keep only samples in subset
if sample_subset != 'all':
    adata = adata[adata.obs[sample_key].isin(sample_subset)]

# %% Add coldata to adata.obs
coldata = pd.read_csv(coldata_path)
coldata = coldata[coldata[sample_key].isin(adata.obs[sample_key])]

# %%
intersecting_columns = list(set(adata.obs) & set(coldata))
print("INFO: Intersecting columns between metadata and adata.obs")
print(intersecting_columns)
# %%
if len(intersecting_columns) > 1:
    new_obs = adata.obs.reset_index().merge(coldata, how = "left", on = intersecting_columns)
else:
    new_obs = adata.obs.reset_index().merge(coldata, how = "left", on = [sample_key])
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
# Crucial for being able to save and read the mudata object later on
# Otherwise you'll get an error message like TypeError: AnnData.__init__() got an unexpected keyword argument 'TAL'
# If one of your cell type views is called EC/TAL
# Or TypeError: AnnData.__init__() got an unexpected keyword argument 'P'
# if one of your cell type  views is called VSM/P
adata.obs[groupby] = adata.obs[groupby].str.replace('/', '_')
# Remove based on low number of samples covered
adata = adata[~adata.obs[groupby].isin(drop_views)]
adata.layers['counts'] = adata.X
# %% Create a multi-view structure
mdata_original = adata_to_views(adata,
                        groupby=groupby,
                        sample_key=sample_key,
                        obs_keys=[batch_key, condition_key],
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
