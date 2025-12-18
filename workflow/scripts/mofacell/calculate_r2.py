# 31.10.2023 Charlotte Boys
# %% Import libraries
import platform as platform
import os as os
import mofax as mofa
import pandas as pd
import warnings

# %% Parse snakemake object
if 'snakemake' in locals():
    model_input = snakemake.input['model_input']
    variance_explained_by_sample_output = snakemake.output['variance_explained_by_sample_output']
    variance_explained_total_output = snakemake.output['variance_explained_total_output']
    condition_key = snakemake.params['condition_key']
    sample_key = snakemake.params['sample_key']
else:
    model_input = 'results/mofacell/model~combined_sn10xLRCluster2.h5ad'
    variance_explained_by_sample_output = 'results/mofacell/r2_by_sample~combined_sn10xLRCluster2.csv'
    variance_explained_total_output = 'results/mofacell/r2_total~combined_sn10xLRCluster2.csv'
    condition_key = 'Cohort_Charlotte'
    sample_key = 'OCEAN_MS_Exp_ID'
    if platform.system() == 'Linux':
        model_input = os.path.join('/mnt/sds-hd/sd22b002/projects', model_input)
        variance_explained_by_sample_output = os.path.join('/mnt/sds-hd/sd22b002/projects/OCEAN_MCFA', variance_explained_by_sample_output)
        variance_explained_total_output = os.path.join('/mnt/sds-hd/sd22b002/projects/OCEAN_MCFA', variance_explained_total_output)

    elif platform.system() == 'Darwin':
        model_input = os.path.join('/Users/charlotteboys/GitHub/OCEAN_MCFA', model_input)
        variance_explained_by_sample_output = os.path.join('/Users/charlotteboys/GitHub/OCEAN_MCFA', variance_explained_by_sample_output)
        variance_explained_total_output = os.path.join('/Users/charlotteboys/GitHub/OCEAN_MCFA', variance_explained_total_output)

# %% Load mofa model
model = mofa.mofa_model(model_input)

# %% Calculate R2 values
warnings.filterwarnings('ignore')
model.metadata[sample_key] = model.metadata.index.values
variance_explained_sample = model.calculate_variance_explained(
        factors=None, # all factors
        groups=None,
        group_label= sample_key,
        views=None,
        per_factor = True
    )
variance_explained_sample_all = model.calculate_variance_explained(
        factors=None, # all factors
        groups=None,
        group_label= sample_key,
        views=None,
        per_factor = False # uses all factors to reconstruct
    )
variance_explained_sample = variance_explained_sample.join(model.metadata[condition_key], on="Group")
variance_explained_sample_all = variance_explained_sample_all.join(model.metadata[condition_key], on="Group")
variance_explained_sample_all['Factor'] = 'All'
variance_explained_by_sample = pd.concat([variance_explained_sample, variance_explained_sample_all], axis = 0, ignore_index = True)

# %% Calculate R2 values without splitting by sample
variance_explained_total = model.calculate_variance_explained(
        factors=None, # all factors
        groups=None,
        group_label= None,
        views=None,
        per_factor = True
    )

variance_explained_total_allfactors = model.calculate_variance_explained(
        factors=None, # all factors
        groups=None,
        group_label= None,
        views=None,
        per_factor = False # uses all factors to reconstruct
    )

variance_explained_total_allfactors['Factor'] = 'All'
variance_explained_total = pd.concat([variance_explained_total, variance_explained_total_allfactors], axis = 0, ignore_index = True)

# %% Save R2 values
variance_explained_by_sample.to_csv(variance_explained_by_sample_output, index=False)
# %%
variance_explained_total.to_csv(variance_explained_total_output, index=False)
# %%
