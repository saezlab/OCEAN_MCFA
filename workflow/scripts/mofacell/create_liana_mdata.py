# 26.04.2023 - Charlotte Boys

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
# 
import liana as li
import decoupler as dc

from helper_functions import adata_to_views, preprocess_data

# %% Parse snakemake object 
if 'snakemake' in locals():
    h5ad_input = snakemake.input['h5ad_input']
    h5mu_output = snakemake.output['h5mu_output']
    sample_key = snakemake.params['sample_key']
    groupby = snakemake.params['groupby']
else:
    sample_key = 'EdgarID'
    groupby = 'celltype'
    h5ad_input = 'results/mofacell/LIANA~JulioNereid~sn10xcanon2.h5ad'

    if platform.system() == 'Linux':
        h5ad_input = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', h5ad_input)

    elif platform.system() == 'Darwin':
        h5ad_input = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', h5ad_input)

# %% Read in AnnData object
print('INFO: Loading data')
adata = sc.read_h5ad(h5ad_input)

# %% 
print('INFO: Creating multi-view structure')
mdata = li.multi.lrs_to_views(adata,
                              sample_key = sample_key,
                              score_key='magnitude_rank',
                              obs_keys=['ID', 'Project1', 'Cohort', 'Prep', 'Method', 'Age', 'Race',
                                        'eGFRatBx_NEPTUNE', 'UPCRatBx_NEPTUNE', 'Sex',
                                        'RAASBlockPreBx', 'InterstitialFibrosis', 'TubularAtrophy'], # add those to mdata.obs
                              lr_prop = 0.1, # minimum required proportion of samples to keep an LR
                              lrs_per_sample = 20, # minimum number of interactions to keep a sample in a specific view
                              lrs_per_view = 20, # minimum number of interactions to keep a view
                              samples_per_view = 10, # minimum number of samples to keep a view
                              min_variance = 0, # minimum variance to keep an interaction
                              lr_fill = 0, # fill missing LR values across samples with this
                              verbose=True
                              )


mdata

if 'snakemake' in locals():
    print('INFO: Writing AnnData to {0}'.format(h5mu_output))
    mdata.write(h5mu_output)