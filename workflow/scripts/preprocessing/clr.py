# %% 
import pandas as pd
import scanpy as sc
import os
import platform
from skbio.stats.composition import clr, multiplicative_replacement

# %% Parse snakemake
if 'snakemake' in locals():
    h5ad_input = snakemake.input['h5ad_input']
    clr_csv = snakemake.output['clr_csv']
    prop_csv = snakemake.output['prop_csv']
    groupby = snakemake.params['groupby']
    sample_key = snakemake.params['sample_key']
else:
    h5ad_input = 'results/preprocessing/Julio_OCEAN_Nereid.h5ad'
    groupby ='celltype' 
    sample_key = 'EdgarID'
    clr_csv = 'results/preprocessing/clr~Julio_OCEAN_Nereid~sn10x.csv'
    prop_csv = 'results/preprocessing/prop~Julio_OCEAN_Nereid~sn10x.csv'

    if platform.system() == 'Linux':
        h5ad_input = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', h5ad_input)
        clr_csv = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', clr_csv)
        prop_csv = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', prop_csv)

    elif platform.system() == 'Darwin':
        h5ad_input = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', h5ad_input)
        clr_csv = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', clr_csv)
        prop_csv = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', prop_csv)

# %%
adata = sc.read_h5ad(h5ad_input)

# %% Regularise groupbys
# Assuming adata is your AnnData object and 'groupby' is the key for the groupby column
adata.obs[groupby] = adata.obs[groupby].str.replace('+', 'pos')
adata.obs[groupby] = adata.obs[groupby].str.replace('-', '_')
adata.obs[groupby] = adata.obs[groupby].str.replace('/', '_')
adata.obs[groupby] = adata.obs[groupby].str.replace(' ', '')

# %%
# Group the DataFrame by sample and count the occurrences of each cell type
grouped = adata.obs.groupby(sample_key)[groupby].value_counts().unstack(fill_value=0)
# Remove rows which are all zero 
grouped = grouped.loc[(grouped != 0).any(axis=1)]
# Calculate the proportions within each group
proportions = grouped.div(grouped.sum(axis=1), axis=0).values

# %%
# Get row and column names
sample_names = grouped.index
groupby_names = grouped.columns

# Create a new DataFrame with proportions and appropriate row and column names
proportions_df = pd.DataFrame(proportions, index=sample_names, columns=groupby_names)

# %% Replace zeros with small positive delta making sure all proportions still add to 1
# Then perform clr transform
clr_mat = clr(multiplicative_replacement(proportions))
clr_df = pd.DataFrame(clr_mat, index=sample_names, columns=groupby_names)

# %%
clr_df.to_csv(clr_csv)
proportions_df.to_csv(prop_csv)

# %%
