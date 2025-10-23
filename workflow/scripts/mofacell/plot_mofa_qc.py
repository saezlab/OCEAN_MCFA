# (C) Robin Fallegger, adapted by Charlotte Boys 16.02.2024
# %%
import scanpy as sc
import platform as platform
import os as os
import yaml as yaml

import decoupler as dc
import pandas as pd

# load muon and mofax
import muon as mu

import numpy as np
import matplotlib.pyplot as plt
import PyComplexHeatmap as pch


# %%
from helpers_multiome_pp import join_mdata, get_psbulk_stats

# %%
# if snakemake object is in local
if 'snakemake' in locals():
     # load original data
    reference_file = snakemake.input[0]
    #projected_file = snakemake.input[1]
    colors = snakemake.params['colors']
    sample_key = snakemake.params['sample_key']
    condition_key = snakemake.params['condition_key']
    omic = snakemake.params['omic']
    files = [reference_file]#, projected_file]

else:
    wildcards = {'dataset': 'Julio_OCEAN_Nereid', 'subset': 'all'}
    omic = 'snRNA'
    sample_key = 'ID'
    condition_key = 'Cohort'
    reference_file = 'results/mofacell/{dataset}~{subset}.h5mu'.format(**wildcards)
    #projected_file = 'results/mofacell/{omic}/{source}/{psbulk}/{proj_omic}_projected.h5mu'.format(**wildcards)
    
    #colors = {'omic': {'snRNA': "#e4ccf1", 'scRNA': "#6a2786"}, condition_key: {'CKD': "#032ca6", 'AKI': '#f25e3d', 'Ref': "#f2dbae"}}

    target_sum = 1e4

    files = []
    for file in [reference_file]:
        if platform.system() == 'Linux':
            files.append( os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', file))
        elif platform.system() == 'Darwin':
            files.append(os.path.join('/Volumes/sd22b002/projects/MichiganIgAN', file))

    for file in files:
        #check if input file exists
        if not os.path.exists(file):
            raise FileNotFoundError(file + ' does not exist')

omics = [omic]#[wildcards['omic'], wildcards['proj_omic']]

# %%
mdatas = []
#check if input file exists
for input_file, omic in zip(files, omics):

    ### Load RNA modalities
    print('INFO: Loading {0} pseudobulked MuData'.format(omic))
    mdata = mu.read_h5mu(input_file)
    mdata.obs['omic'] = omic

    # reindex obs to have unique index (patient id + omic)
    mdata.obs.index = mdata.obs.index + '_' + mdata.obs['omic']
    mdata.obs[sample_key] = mdata.obs.index.str.split('_').str[0]
    
    mdata.uns['psbulk_stats'].index = mdata.uns['psbulk_stats'].index + '_' + omic

    # do the same for each view
    for view in mdata.mod.keys():
        mdata.mod[view].obs.index = mdata.mod[view].obs.index + '_' + omic

    mdatas.append(mdata)

# %%
# We don't need multiome right now
#mdata, mdata_obs = join_mdata(mdatas[0], mdatas[1], genes='left', views='inner', obsm_keys= ['X_mofa'], uns_keys=['psbulk_stats', 'mofa'])
print(mdata)

# %%
#mu.pp.filter_obs(mdata, condition_key, lambda x: x != 'AKI')

# %%
cells = get_psbulk_stats(mdata, 'cells')
cells = np.log10(cells)
counts = get_psbulk_stats(mdata, 'counts')
counts = np.log10(counts)

# %%
n_genes = []
for view in mdata.mod.keys():
    # Counts the number of non-nans - i.e. all the genes in the view
    # Check with Robin that this is what he was aiming for
    # This doesn't change depending on the sample
    df = pd.DataFrame(np.count_nonzero(~np.isnan(mdata.mod[view].X), axis = 1),
                                 index = mdata.mod[view].obs_names,
                                 columns=[view])
    # df = df / mdata.mod[view].shape[1]
    n_genes.append(df)
    
n_genes = pd.concat(n_genes, axis = 1)
n_genes = np.log10(n_genes)

# %%
#condition_color = {cond: colors[condition_key].get(cond, None) for cond in colors[condition_key].keys() if cond in mdata.obs[condition_key].unique()}
#omic_color = {cond: colors['omic'].get(cond, None) for cond in colors['omic'].keys() if cond in mdata.obs['omic'].unique()}

row_ha = pch.HeatmapAnnotation(Condition = pch.anno_simple(mdata.obs.filter([condition_key], axis=1).astype('object'), label = True),#, colors=condition_color),
                                Omic = pch.anno_simple(mdata.obs.filter(['omic'], axis=1), label = False),#, colors=omic_color),
                                vgap = 3,
                                legend_width=10,
                                label_side= 'bottom',
                                legend = False,
                                axis = 0, verbose=0)

cm_genes = pch.ClusterMapPlotter(data=n_genes,
                            right_annotation=row_ha,
                            col_cluster=False,row_cluster=True,
                            label='log10(genes)',row_dendrogram=True,
                            col_dendrogram = False,
                            show_rownames=False,show_colnames=True,
                            verbose=0,legend_gap=5,
                            cmap = 'Reds', plot=False, na_col = 'grey')

cm_cells = pch.ClusterMapPlotter(data=cells,
                            right_annotation=row_ha,
                            col_cluster=False,row_cluster=False,
                            label='log10(cells)',row_dendrogram=False,
                            col_dendrogram = False,
                            show_rownames=False,show_colnames=True,
                            verbose=0,legend_gap=5,
                            cmap = 'Greens', plot=False, na_col = 'grey')

cm_counts = pch.ClusterMapPlotter(data=counts,
                            right_annotation=row_ha,
                            col_cluster=False,row_cluster=False,
                            label='log10(counts)',row_dendrogram=False,
                            col_dendrogram = False,
                            show_rownames=False,show_colnames=True,
                            verbose=0,legend_gap=5,
                            cmap = 'Blues',  plot=False, na_col = 'grey')

# %%
plt.figure(figsize=(20, 6))
ax,legend_axes = pch.composite(cmlist = [cm_genes, cm_cells, cm_counts], main = 0, axis = 1, width_ratios=[1.1,1,1], col_gap=7, legend_gap=15)
cm_genes.ax.set_title('Genes')
cm_cells.ax.set_title('Cells')
cm_counts.ax.set_title('Counts')
plt.suptitle('Pseudo-bulk stats', fontsize = 14)
if 'snakemake' in locals():
    plt.savefig(snakemake.output[0], dpi = 300)