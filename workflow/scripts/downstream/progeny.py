# 07.02.2024 Charlotte Boys
# %% Import libraries
import platform as platform
import os as os
import decoupler as dc
import pandas as pd
import numpy as np

# %% Parse snakemake object
if 'snakemake' in locals():
    loadings = snakemake.input['loadings']
    acts_csv = snakemake.output['acts_csv']
    pvals_csv = snakemake.output['pvals_csv']
else:
    loadings = "results/mofacell/loadings~combined_sn10xLRCluster2.csv"# 'results/mofacell/loadings~Julio_OCEAN_Nereid~LIANA_all~all.csv'
    acts_csv = 'results/downstream/progeny_acts~combined_sn10xLRCluster2.csv'
    pvals_csv = 'results/downstream/progeny_pvals~combined_sn10xLRCluster2.csv'
    if platform.system() == 'Linux':
        loadings = os.path.join('/mnt/sds-hd/sd22b002/projects/OCEAN_MCFA', loadings)
        acts_csv = os.path.join('/mnt/sds-hd/sd22b002/projects/OCEAN_MCFA', acts_csv)
        pvals_csv = os.path.join('/mnt/sds-hd/sd22b002/projects/OCEAN_MCFA', pvals_csv)
    elif platform.system() == 'Darwin':
        loadings = os.path.join('/Users/charlotteboys/GitHub/OCEAN_MCFA', loadings)
        acts_csv = os.path.join('/Users/charlotteboys/GitHub/OCEAN_MCFA', acts_csv)
        pvals_csv = os.path.join('/Users/charlotteboys/GitHub/OCEAN_MCFA', pvals_csv)

# %%
loadings_all = pd.read_csv(loadings, index_col=0)
print(loadings_all.head())
# %%
progeny = dc.get_progeny(organism='human', top = 500)
print(progeny.head())
# %% We run by cell type because each cell type has a different feature set
def pathway_activity(cell_type):
    loadings_df = loadings_all[loadings_all['view'] == cell_type]
    mat = loadings_df.pivot(index='Factor', columns='variable', values='Loading')
    pathway_acts, pathway_pvals = dc.run_ulm(mat=mat, net=progeny, verbose=True)
    return (pathway_acts.melt(ignore_index = False).assign(View = cell_type),
            pathway_pvals.melt(ignore_index = False).assign(View = cell_type))
# %%
views = np.unique(loadings_all['view'])
cell_types = [x for x in views if '&' not in x]
cell_types
# %%
results = [pathway_activity(cell_type) for cell_type in cell_types]
acts = pd.concat([result[0] for result in results], axis=0)
pvals = pd.concat([result[1] for result in results], axis=0)

# %%
acts.to_csv(acts_csv)
pvals.to_csv(pvals_csv)
# %%
