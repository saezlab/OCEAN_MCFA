# 31.10.2023 Charlotte Boys
# %% Import libraries
import platform as platform
import os as os
import yaml as yaml

import pandas as pd

import liana as li

# load muon and mofax
import muon as mu
import mofax as mofa


# %% Parse snakemake object TODO: update for new repo
if 'snakemake' in locals():
    h5mu_input = snakemake.input['h5mu_input']
    model_output = snakemake.output['model_output']
    scores_output = snakemake.output['scores_output']
    scores_wide_output = snakemake.output['scores_wide_output']
    loadings_output = snakemake.output['loadings_output']
    nfactors = snakemake.params['nfactors']
    hvg = snakemake.params['hvg']
    sample_key = snakemake.params['sample_key']

else:
    h5mu_input = 'MichiganIgAN/results/mofacell/Julio_OCEAN_Nereid~all.h5mu'
    model_output = 'MichiganIgAN/results/mofacell/model~Julio_OCEAN_Nereid~all~hvg.h5ad'
    scores_output = 'MichiganIgAN/results/mofacell/scores~Julio_OCEAN_Nereid~all~hvg.csv'
    scores_wide_output = 'MichiganIgAN/results/mofacell/widescores~Julio_OCEAN_Nereid~all~hvg.csv'
    loadings_output = 'MichiganIgAN/results/mofacell/loadings~Julio_OCEAN_Nereid~all~hvg.csv'
    nfactors = 20
    hvg = 'hvg'
    sample_key = 'EdgarID'

    if platform.system() == 'Linux':
        h5mu_input = os.path.join('/mnt/sds-hd/sd22b002/projects', h5mu_input)
        model_output = os.path.join('/mnt/sds-hd/sd22b002/projects', model_output)
        scores_output = os.path.join('/mnt/sds-hd/sd22b002/projects', scores_output)
        scores_wide_output = os.path.join('/mnt/sds-hd/sd22b002/projects', scores_wide_output)
        loadings_output = os.path.join('/mnt/sds-hd/sd22b002/projects', loadings_output)

    elif platform.system() == 'Darwin':
        h5mu_input = os.path.join('/Users/charlotteboys/GitHub', h5mu_input)
        model_output = os.path.join('/Users/charlotteboys/GitHub', model_output)
        scores_output = os.path.join('/Users/charlotteboys/GitHub', scores_output)
        scores_wide_output = os.path.join('/Users/charlotteboys/GitHub', scores_wide_output)
        loadings_output = os.path.join('/Users/charlotteboys/GitHub', loadings_output)

# %% Read mdata
print('INFO: Loading assay')
mdata = mu.read_h5mu(h5mu_input)

# %% Set up to run with HVGs
# First of all, if we're running with LIANA, the hvg selection has already been done
if any("&" in key for key in mdata.mod.keys()):
    use_var = None
elif hvg == 'all':
    use_var = None
elif hvg == 'hvg':
    use_var = 'highly_variable' # I think there's a bug in muon that this doesn't actually make a difference
    # So we filter to only keep HVGs:
    for view in mdata.mod.keys():
        mu.pp.filter_var(mdata.mod[view], var = 'highly_variable')
    
# %% Run mofa
print('INFO: Running MOFA')
mu.tl.mofa(mdata,
           use_obs='union',
           convergence_mode='medium',
           n_factors=nfactors,
           seed=1337,
           outfile=model_output,
           use_var=use_var
           )
print(mdata)
print(mdata.varm)
# %% Extract factor scores and loadings
print('INFO: Extracting MOFA factor scores and loadings')
factor_scores_obs = li.utils.get_factor_scores(mdata, obsm_key='X_mofa')

print(factor_scores_obs.head())
factor_scores_wide = factor_scores_obs.set_index(sample_key).filter(like='Factor', axis=1)
factor_scores_long = factor_scores_wide.rename_axis('Sample').reset_index().melt(id_vars=['Sample'],
                                                           var_name='Factor',
                                                           value_name='Score')

variable_loadings =  li.utils.get_variable_loadings(mdata, varm_key = 'LFs', view_sep=':')
print(variable_loadings.head())
variable_loadings_long = pd.melt(variable_loadings, 
                                 id_vars=['view', 'variable'],
                                 var_name='Factor',
                                 value_name='Loading')
# %% Save factor scores and loadings
print('INFO: Saving MOFA factor scores and loadings')
factor_scores_long.to_csv(scores_output)
factor_scores_wide.to_csv(scores_wide_output)
variable_loadings_long.to_csv(loadings_output)
# %%
