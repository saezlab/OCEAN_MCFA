# %% Import libraries
import platform as platform
import os as os
import scanpy as sc

# %%
os.chdir('/mnt/sds-hd/sd22b002/projects/OCEAN_MCFA/workflow/scripts/mofacell')
from helper_functions import adata_to_views, preprocess_data

# %% Parse snakemake object
if 'snakemake' in locals():
    h5ad_input = snakemake.input['h5ad_input']
    sample_key = snakemake.params['sample_key']
    condition_key = snakemake.params['condition_key']
    batch_key = snakemake.params['batch_key']
    groupby = snakemake.params['groupby']
    min_cells = snakemake.params['min_cells']
    pod_csv = snakemake.output['pod_csv']
    if snakemake.wildcards['subset'] == 'all':
        sample_subset = 'all'
        cell_subset = 'all'
    else:
        sample_subset = snakemake.params['sample_subset']
        cell_subset = snakemake.params['cell_subset']

else:
    h5ad_input = 'OCEAN_MCFA/results/mofacell/sn10xLRCluster2.h5ad'#data/OCEAN_v3_Nu_102025a_CZI.h5ad'
    sample_key = 'OCEAN_MS_Exp_ID'
    condition_key = 'Cohort_Manuscript'
    batch_key = 'Project1'
    groupby = 'LR_Cluster'
    min_cells = 20
    pod_csv = 'OCEAN_MCFA/results/mofacell/pod_psbulk.csv'
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

    if platform.system() == 'Linux':
        h5ad_input = os.path.join('/mnt/sds-hd/sd22b002/projects', h5ad_input)
        pod_csv = os.path.join('/mnt/sds-hd/sd22b002/projects', pod_csv)
# %%
adata = sc.read_h5ad(h5ad_input)
# # %%
# adata = adata[adata.obs[sample_key].isin(sample_subset)]

# # %%
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# # %%
# adata.obs[groupby] = adata.obs[groupby].str.replace('/', '_')
# adata.layers['counts'] = adata.X
# %% Store pseudobulk counts in mdata object
mdata = adata_to_views(adata,
                        groupby=groupby,
                        sample_key=sample_key,
                        obs_keys=[condition_key, batch_key],
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

# %% Enforce minimum number of views per sample:
# lambda function to check if sample is in a view
in_view = lambda view: mdata.obs_names.isin([*mdata.mod[view].obs_names])*1
# Check for two or more views per sample
mask = sum(list(map(in_view, [*mdata.mod.keys()]))) >= 2
pass_samples = mdata.obs_names[mask]
for current_mod in mdata.mod.keys():
    s = mdata.mod[current_mod].obs_names.isin(pass_samples)
    mdata.mod[current_mod] = mdata.mod[current_mod][s, :].copy()
mdata.update()
mdata.uns['psbulk_stats'] = mdata.uns['psbulk_stats'][mask]
# %%
mdata_log1p = preprocess_data(mdata.copy(),
                        target_sum = 10000,
                        remove_mitochondrial = False,
                        exclude_highly_expressed = False,
                        batch_key = batch_key,
                        n_top_genes = None) # default number of top genes

# %%
df = mdata_log1p.mod['POD'].to_df().rename(columns=lambda x: x.replace('POD:', ''))
df.to_csv(pod_csv)
# %%
# for key in mdata.mod.keys():
#     df = mdata.mod[key].to_df().rename(columns=lambda x: x.replace(key + ':', ''))
#     df.to_csv('/Users/charlotteboys/GitHub/MichiganIgAN/results/downstream/psbulk_expression_' + key + '.csv')