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


# %% Parse snakemake object
if 'snakemake' in locals():
    h5mu_input = snakemake.input['h5mu_input']
    h5mu_liana = snakemake.input['h5mu_liana']
    h5mu_output = snakemake.output['h5mu_output']
    cell_subset = snakemake.params['cell_subset']
    hvg = snakemake.params['hvg']
else:
    h5mu_input = 'MichiganIgAN/results/mofacell/Julio_OCEAN_Nereid~sn10xcanon2.h5mu'
    h5mu_liana = 'MichiganIgAN/results/mofacell/LIANA~Julio_OCEAN_Nereid~sn10xcanon2.h5mu'
    cell_subset = 'all' #['CNT', 'DCT', 'IC']
    hvg = 'hvg' #'all'

    if platform.system() == 'Linux':
        h5mu_input = os.path.join('/mnt/sds-hd/sd22b002/projects', h5mu_input)
        h5mu_liana = os.path.join('/mnt/sds-hd/sd22b002/projects', h5mu_liana)

    elif platform.system() == 'Darwin':
        h5mu_input = os.path.join('/Users/charlotteboys/GitHub', h5mu_input)
        h5mu_liana = os.path.join('/Users/charlotteboys/GitHub', h5mu_liana)

# %% Read mdata
print('INFO: Loading assay')
mdata = mu.read_h5mu(h5mu_input)
mdata_liana = mu.read_h5mu(h5mu_liana)

# %%
print('INFO: Subsetting mudata objects')
# Make sure we're looking at the same samples
mu.pp.filter_obs(mdata, var = mdata_liana.obs_names)
mu.pp.filter_obs(mdata_liana, var = mdata.obs_names)

if cell_subset != 'all':
    to_exclude = [view for view in mdata_liana.mod.keys() if not any(string in view.split('&') for string in cell_subset)]
    for view in to_exclude:
        if view in list(mdata_liana.mod.keys()):
            print('Removing {0} due to exclusion list'.format(view))
            mdata_liana.mod[view].var['to_remove'] =  mdata_liana.mod[view].var_names.str.startswith(view)
            mu.pp.filter_var(mdata_liana.mod[view], 'to_remove', lambda x: ~x)
            del mdata_liana.mod[view]
# %%
# to_exclude = [view for view in mdata.mod.keys() if not view in cell_subset]
# for view in to_exclude:
#     if view in list(mdata.mod.keys()):
#         print('Removing {0} due to exclusion list'.format(view))
#         mdata.mod[view].var['to_remove'] =  mdata.mod[view].var_names.str.startswith(view)
#         mu.pp.filter_var(mdata.mod[view], 'to_remove', lambda x: ~x)
#         del mdata.mod[view]
# %% Set up to run with HVGs
if hvg == 'hvg':
    for view in mdata.mod.keys():
        mu.pp.filter_var(mdata.mod[view], var = 'highly_variable')

# %% Filter interactions for those where one or both of the ligand and receptor are HVGs in their respective views
if hvg == 'hvg':
    for view in mdata_liana.mod.keys():
        sender = view.split('&')[0]
        receiver = view.split('&')[1]
        hvg_ligands = [i[1] for i in mdata.mod[sender].var_names.str.split(':')]
        hvg_receptors = [i[1] for i in mdata.mod[receiver].var_names.str.split(':')]
        interactions = [i[1] for i in mdata_liana.mod[view].var_names.str.split(':')]
        #ligands = [i.split('^')[0] for i in interactions]
        #receptors = [i.split('^')[1] for i in interactions]
        mask = [
                    interaction.split('^')[0] in hvg_ligands or interaction.split('^')[1] in hvg_receptors
                    for interaction in interactions
                ]
        mdata_liana.mod[view].var['highly_variable'] = mask
        mu.pp.filter_var(mdata_liana.mod[view], var = 'highly_variable')
# %%
print('INFO: Combining mudata')
# Initialize an empty dictionary to store AnnData objects
adict = {}

# Iterate through the keys of mdata_liana.mod
for key in mdata_liana.mod.keys():
    # Retrieve the AnnData object corresponding to the key
    adict[key] = mdata_liana.mod[key]
for key in mdata.mod.keys():
    adict[key] = mdata.mod[key]

mdata_combined = mu.MuData(adict)
mdata_combined.obs = mdata.obs

mdata_combined
# %%
if 'snakemake' in locals():
    print('INFO: Writing AnnData to {0}'.format(h5mu_output))
    mdata_combined.write(h5mu_output)