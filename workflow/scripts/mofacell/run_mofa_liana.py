# 26.04.2024 Charlotte Boys
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
