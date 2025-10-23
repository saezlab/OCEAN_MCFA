# TODO finish this script
# - [ ] write the rule

# %% Import libraries
import liana as li
import decoupler as dc
import pandas as pd
import platform
import os
import scanpy as sc

# %% Parse snakemake object
if 'snakemake' in locals():
    scores_wide = snakemake.input['scores_wide']
    metadata_csv = snakemake.input['metadata_csv']
    clr = snakemake.input['clr']
    r2_csv = snakemake.input['r2_csv']
    r2_total_csv = snakemake.input['r2_total_csv']
    associations_csv = snakemake.output['associations_csv']
    factors_obs_csv = snakemake.output['factors_obs_csv']

else:
    scores_wide = 'results/mofacell/widescores~Julio_OCEAN_Nereid~sn10xcanon2~hvg.csv'
    metadata_csv = 'results/preprocessing/updated_metadata_clean.csv'
    clr = 'results/preprocessing/clr~Julio_OCEAN_Nereid~sn10x.csv'
    r2_csv = 'results/mofacell/r2_by_sample~Julio_OCEAN_Nereid~sn10xcanon2~hvg.csv'
    r2_total_csv = 'results/mofacell/r2_total~Julio_OCEAN_Nereid~sn10xcanon2~hvg.csv'
    associations_csv = 'results/mofacell/associations~Julio_OCEAN_Nereid~sn10xcanon2~hvg.csv'
    factors_obs_csv = 'results/mofacell/factors_obs~Julio_OCEAN_Nereid~sn10xcanon2~hvg.csv'

    if platform.system() == 'Linux':
        scores_wide = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', scores_wide)
        metadata_csv = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', metadata_csv)
        clr = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', clr)
        r2_csv = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', r2_csv)
        r2_total_csv = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', r2_total_csv)
        associations_csv = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', associations_csv)
        factors_obs_csv = os.path.join('/mnt/sds-hd/sd22b002/projects/MichiganIgAN', factors_obs_csv)

    elif platform.system() == 'Darwin':
        scores_wide = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', scores_wide)
        metadata_csv = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', metadata_csv)
        clr = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', clr)
        r2_csv = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', r2_csv)
        r2_total_csv = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', r2_total_csv)
        associations_csv = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', associations_csv)
        factors_obs_csv = os.path.join('/Users/charlotteboys/GitHub/MichiganIgAN', factors_obs_csv)


# %% Read in factor scores
factor_scores = pd.read_csv(scores_wide, index_col=0).rename_axis('Sample').filter(like='Factor', axis=1)

# %% Read in proportions
clr_props = pd.read_csv(clr, index_col=0).reindex(factor_scores.index)
clr_props.columns = [col + '_clr' for col in clr_props.columns]

# %% Read in metadata
metadata = pd.read_csv(metadata_csv, index_col=0).reindex(factor_scores.index)

# %% Save data for plotting
obs = pd.concat([metadata, clr_props], axis='columns')
to_plot = pd.concat([factor_scores, obs], axis = 1)
to_plot.to_csv(factors_obs_csv)

# %% Read in r2
r2 = pd.read_csv(r2_csv).rename(columns={'Group': 'Sample'})
r2 = r2[~(r2['Factor'] == 'All')]
r2_total = pd.read_csv(r2_total_csv)

# %% Find best-explained views for each Factor
bev = r2_total.loc[r2_total.groupby('Factor')['R2'].idxmax()]
# %% Find the set of samples where R2 is NA for the best explained view for each factor in r2
na_samples = r2.merge(bev[['Factor', 'View']], on=['Factor', 'View'], how='inner')
na_samples = na_samples[na_samples['R2'].isna()]


# %% Remove samples with NA R2 for the best explained view for each factor
for row in zip([*na_samples['Sample']], [*na_samples['Factor']]):
    factor_scores.loc[row] = pd.NA

# %% Run associations
associations = dc.get_metadata_associations([factor_scores, obs], obs_keys = [*obs.columns])

# %% Save data
associations.to_csv(associations_csv)



# %%
